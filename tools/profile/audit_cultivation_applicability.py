"""Reproducible coverage audit separating adapters, algebra and fast execution."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import stim
from constructed_folded import make_model
from parity_envelope_audit import audit as audit_parity
from parity_envelope_audit import axes, certify, layers, operations
from sahay_gadget import SahayGadget
from study_sahay_gadget import entry_validation
from study_wan_corpus import fresh_reset_wires
from study_wan_scaling import (
    REFERENCE,
    family,
    load_reference,
    masks_for,
    permute_masks,
    scope_plan,
)
from wan_corpus import CASES, REVISION, Exporter, digest, export_corpus, pauli_text, verify_source


def rank(rows: list[int]) -> int:
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)


def css_summary(checks: list[dict]) -> dict:
    data = sorted({q for check in checks for q in check["support"]})
    rows = {
        axis: [sum(1 << q for q in check["support"]) for check in checks if check["axis"] == axis]
        for axis in "XZ"
    }
    ranks = {axis: rank(rows[axis]) for axis in "XZ"}
    commutes = all((x & z).bit_count() % 2 == 0 for x in rows["X"] for z in rows["Z"])
    return dict(
        data_qubits=len(data),
        physical_data=data,
        check_counts={axis: len(rows[axis]) for axis in "XZ"},
        independent_ranks=ranks,
        redundant_checks=len(checks) - sum(ranks.values()),
        commutes=commutes,
        logical_qubits=len(data) - sum(ranks.values()) if commutes else None,
    )


def css_from_text(text: str) -> list[dict]:
    # Pure-X/Z products select the code checks; the supplied logical magic
    # readouts are mixed/Y products or are inside the final inverse-T pair.
    result: list[dict] = []
    lines = text.splitlines()
    last_phase = max(
        (i for i, line in enumerate(lines) if line.startswith(("T ", "T_DAG "))), default=-1
    )
    for i, line in enumerate(lines):
        if not line.startswith("MPP "):
            continue
        for group in stim.Circuit(line)[0].target_groups():
            if all(t.is_x_target for t in group):
                axis = "X"
            elif all(t.is_z_target for t in group):
                axis = "Z"
            else:
                continue
            result.append(dict(axis=axis, support=[t.value for t in group], line=i + 1))
    # Chan and SOFT d3 put an all-X logical measurement before the final T
    # inverse. Their CSS checks come after that inverse. SOFT d5 instead has
    # its CSS checks before a final phase layer and mixed-Y logical readout.
    later = [check for check in result if check["line"] > last_phase + 1]
    return later if later else result


def canonical_fold_d3(model: dict) -> tuple[dict, dict]:
    """Pull back a noiseless decoder; preserve all original measurement records."""
    original = model["events"]
    first_native = next(i for i, e in enumerate(original) if e["name"] in {"CSX", "CCZ"})
    last_native = max(i for i, e in enumerate(original) if e["name"] in {"CSX", "CCZ"})
    root = model["ghz_root"]
    cat_readout = next(
        i for i in range(last_native + 1, len(original)) if original[i]["name"] == "POST_Z"
    )
    decode_start = cat_readout
    while original[decode_start]["name"] == "POST_Z":
        decode_start += 1
    projection = next(
        i for i in range(decode_start, len(original)) if original[i]["name"] == "POST_Z"
    )
    decoder = stim.Circuit()
    for e in original[decode_start:projection]:
        if e["name"] not in {"H", "CX"}:
            raise ValueError("decoder is not a noiseless Clifford circuit")
        decoder.append(e["name"], e["support"])
    width = decoder.num_qubits
    inverse = decoder.to_tableau().inverse()

    def pulled(q: int, axis: str) -> Any:
        p = stim.PauliString(width)
        p[q] = axis
        value = inverse(p)
        if value.sign != 1:
            raise ValueError("decoder introduces a negative Pauli")
        return value

    logical = {axis: axes(pulled(model["logical"], axis)) for axis in "XZ"}
    if any(set(logical[axis].values()) != {axis} for axis in "XZ"):
        raise ValueError("decoder does not expose pure CSS logical axes")
    result = copy.deepcopy(model)
    events = []
    record = 0
    block = 0
    for i, e in enumerate(original):
        if decode_start <= i < projection:
            continue
        event = copy.deepcopy(e)
        event["original_event"] = i
        if event["name"] == "POST_Z":
            q = event["support"][0]
            if i >= projection:
                p = axes(pulled(q, "Z"))
                if len(set(p.values())) != 1 or next(iter(p.values())) not in "XZ":
                    raise ValueError("decoder measurement is not a pure CSS check")
                event.update(name="MEASURE_PRODUCT", axis=next(iter(p.values())), support=list(p))
            else:
                event.update(name="MEASURE", axis="Z", inverted=False)
            event["record"] = record
            events.append(event)
            events.append(dict(name="DETECTOR", records=[record], reference_parity=0))
            record += 1
            if first_native < i < last_native and q == root:
                block += 1
        else:
            if event["name"] in {"CSX", "CS_DAG_X", "CCZ"}:
                event["native_block"] = block
            events.append(event)
    if block != 1:
        raise ValueError("did not find exactly two folded checks")
    result["events"] = events
    del result["logical"]
    result["logical_x"] = {str(q): a for q, a in logical["X"].items()}
    result["logical_z"] = {str(q): a for q, a in logical["Z"].items()}
    return result, dict(
        decoder_start_event=decode_start,
        terminal_projection_event=projection,
        decoder_elementary_gates=projection - decode_start,
        original_records=record,
        description="noiseless decoder pulled back to CSS products and logical axes; "
        "all records retained",
    )


def check_rewrite(clifft: Any, original: dict, rewritten: dict) -> dict:
    rng = np.random.default_rng(4950928)
    error, comparisons = 0.0, 0
    cases = []
    for index in range(4):
        faults = {
            channel["location"]: int(rng.integers(1, 4 ** len(channel["support"])))
            if channel["kind"] != "flip"
            else 1
            for channel in original["channels"]
            if index and rng.random() < 0.02
        }
        texts = [Exporter(m, 0, faults).export()[0] for m in (original, rewritten)]
        samples = clifft.sample(
            clifft.compile(texts[0]), 8, seed=270 + index, threads=1, batch_size=1
        )
        values = [
            clifft.record_probabilities(clifft.compile(fresh_reset_wires(t)), samples.measurements)
            for t in texts
        ]
        difference = float(np.max(np.abs(values[0] - values[1])))
        error = max(error, difference)
        comparisons += len(values[0])
        cases.append(
            dict(
                faults=faults,
                rejected_records=int(np.any(samples.detectors, axis=1).sum()),
                maximum_absolute_error=difference,
            )
        )
    if error > 1e-11:
        raise AssertionError("Noiseless decoder rewrite changed record probabilities")
    return dict(cases=cases, comparisons=comparisons, maximum_absolute_error=error)


def folded_check(model: dict, reference: Any) -> tuple[dict, Any]:
    if not all(isinstance(e, dict) for e in model["events"]):
        return dict(selected=False, layer="input representation", reason="list-event schema"), None
    try:
        gadget = SahayGadget(model, reference)
    except (ValueError, KeyError) as error:
        return dict(selected=False, layer="model recognizer", reason=str(error)), None
    plans = family(gadget.masks, gadget.rank)
    return dict(
        selected=True,
        certificate=gadget.certificate,
        peak_entries=plans["peak_entries"],
        numeric_payload_bytes=plans["numeric_payload_bytes"],
        plans_over_budget=plans["plans_exceeding_existing_limit"],
    ), gadget


def load_pattern(reference: Path) -> Any:
    folder = reference / "src/python/clifft/_folded"
    for path in folder.glob("*.py"):
        relative = "src/python/clifft/_folded/" + path.name
        if path.read_bytes() != subprocess.check_output(
            ["git", "show", REFERENCE + ":" + relative]
        ):
            raise ValueError("Unpinned folded reference " + path.name)
    spec = importlib.util.spec_from_file_location(
        "applicability_folded_reference",
        folder / "__init__.py",
        submodule_search_locations=[str(folder)],
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parity_entry(clifft: Any, text: str, checks: list[dict], envelopes: dict) -> list[dict]:
    data = {q for check in checks for q in check["support"]}
    results = []
    for candidate in envelopes["candidate_envelopes"]:
        if not candidate["accepted"]:
            continue
        certificate = candidate["certificate"]
        row: dict = dict(first_line=candidate["first_line"])
        if set(certificate["phase_qubits"]) != data:
            row.update(
                tested=False,
                reason="phase support differs from terminal code; "
                "needs its intermediate code boundary",
            )
            results.append(row)
            continue
        spectators = sorted(set(certificate["middle_qubits"]) - data)
        probes = [dict.fromkeys(check["support"], check["axis"]) for check in checks]
        probes += [{q: axis} for q in spectators for axis in "XYZ"]
        # These odd-length colour codes have all-data X/Z as logical axes.
        # The all-Y sign does not affect this purity test.
        if len(data) % 2 != 1 or any(len(check["support"]) % 2 for check in checks):
            raise ValueError("All-data logical probe requires this odd even-check CSS code")
        probes += [dict.fromkeys(sorted(data), axis) for axis in "XYZ"]
        prefix = "\n".join(text.splitlines()[: candidate["first_line"] - 1]) + "\n"
        prefix += "\n".join("EXP_VAL " + pauli_text(p) for p in probes) + "\n"
        program = clifft.compile(prefix)
        writes = [line for line in program.inspect().splitlines() if "WRITE_EXPECTATION" in line]
        if len(writes) != len(probes):
            raise AssertionError("Probe output accounting changed")
        classical = [" I sign=" in line for line in writes]
        fixed = {
            q: [axis for j, axis in enumerate("XYZ") if classical[len(checks) + 3 * i + j]]
            for i, q in enumerate(spectators)
        }
        selected = list(range(len(checks))) + [
            len(checks) + 3 * i + "XYZ".index(fixed[q][0])
            for i, q in enumerate(spectators)
            if fixed[q]
        ]
        sample = clifft.sample(program, 128, seed=8027, threads=1, batch_size=1)
        deviation = float(np.max(np.abs(np.abs(sample.exp_vals[:, selected]) - 1)))
        purity_error = float(np.max(np.abs((sample.exp_vals[:, -3:] ** 2).sum(axis=1) - 1)))
        row.update(
            tested=True,
            code_checks_are_tracked_signs=all(classical[: len(checks)]),
            spectator_axes=fixed,
            all_spectators_have_tracked_axis=all(fixed.values()),
            all_selected_probes_at_width_one=all(" w1 " in writes[i] for i in selected),
            selected_probe_lines=[writes[i] for i in selected],
            prefix_peak_active_width=program.peak_active_width,
            sampled_shots=128,
            maximum_signed_probe_deviation=deviation,
            maximum_logical_purity_error=purity_error,
            includes_entry_flags=True,
            scope="entry code and single-qubit spectator probes; "
            "no full gadget/noise binding or continuation",
        )
        indices = {q: i for i, q in enumerate(sorted(data))}
        x_rows = [
            sum(1 << indices[q] for q in check["support"])
            for check in checks
            if check["axis"] == "X"
        ]
        masks = masks_for(x_rows, len(data), [])
        order = scope_plan(masks, len(x_rows))["order"]
        plans = family(permute_masks(masks, order), len(x_rows))
        row["unary_code_plan_screen"] = dict(
            peak_entries=plans["peak_entries"],
            numeric_payload_bytes=plans["numeric_payload_bytes"],
            plans_over_budget=plans["plans_exceeding_existing_limit"],
            scope="single-bit factors on this code; excludes actual fault binding and continuation",
        )
        results.append(row)
    return results


def representation_checks(model: dict, reference: Any) -> list[dict]:
    from wan_corpus import gate_lines

    variants = {}
    no_annotations = copy.deepcopy(model)
    for event in no_annotations["events"]:
        event.pop("native_block", None)
    variants["native block annotations removed"] = no_annotations
    decomposed = copy.deepcopy(model)
    expanded = []
    for event in decomposed["events"]:
        if event["name"] in {"CSX", "CS_DAG_X"}:
            for line in gate_lines(event["name"], event["support"]):
                name, *qubits = line.split()
                expanded.append(dict(event, name=name, support=list(map(int, qubits))))
        else:
            expanded.append(event)
    decomposed["events"] = expanded
    variants["exact native-gate decomposition with original noise sites"] = decomposed
    renamed = copy.deepcopy(model)
    used = sorted(
        {q for event in renamed["events"] for q in event.get("support", [])}
        | {q for channel in renamed["channels"] for q in channel["support"]}
    )
    mapping = {q: 3 * i + 11 for i, q in enumerate(reversed(used))}
    for event in renamed["events"]:
        if "support" in event:
            event["support"] = [mapping[q] for q in event["support"]]
    for channel in renamed["channels"]:
        channel["support"] = [mapping[q] for q in channel["support"]]
    for axis in "xyz":
        renamed["logical_" + axis] = {
            str(mapping[int(q)]): value for q, value in renamed["logical_" + axis].items()
        }
    variants["qubit identifiers permuted"] = renamed
    results = []
    for name, variant in variants.items():
        outcome, _ = folded_check(variant, reference)
        expected = name != "native block annotations removed"
        if outcome["selected"] != expected:
            raise AssertionError("Representation control changed: " + name)
        results.append(dict(change=name, outcome=outcome))
    return results


def negative_checks(clifft: Any, text: str) -> list[dict]:
    ops = operations(clifft, text)
    phases = layers(ops)
    forward, middle, reverse = [], [], []
    for (a, b), (c, d) in zip(phases, phases[1:]):
        if b - a < 2:
            continue
        try:
            certify(ops[a:b], ops[b:c], ops[c:d])
        except ValueError:
            continue
        forward, middle, reverse = ops[a:b], ops[b:c], ops[c:d]
        break
    if not middle:
        raise AssertionError("No positive flagged control")
    from dataclasses import replace

    variants = {}
    missing = middle.copy()
    del missing[next(i for i, op in enumerate(missing) if op.name == "CX")]
    variants["noninverse CNOT network"] = (forward, missing, reverse)
    wrong = middle.copy()
    index = next(i for i, op in enumerate(wrong) if op.name == "M")
    wrong[index] = replace(wrong[index], name="MX")
    variants["wrong flag readout basis"] = (forward, wrong, reverse)
    variants["unclosed flag"] = (forward, middle[:index] + middle[index + 1 :], reverse)
    variants["noninverse phase layer"] = (forward, middle, reverse[:-1])
    results = []
    for name, arguments in variants.items():
        try:
            certify(*arguments)
        except ValueError as error:
            results.append(dict(mutation=name, reason=str(error)))
        else:
            raise AssertionError("Malformed parity envelope accepted: " + name)
    return results


def main() -> None:
    import clifft

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    verify_source(args.source)
    args.work.mkdir(parents=True, exist_ok=False)
    corpus = args.work / "corpus"
    manifest = export_corpus(args.source, corpus, 0.001)
    reference = load_reference(args.reference)
    pattern = load_pattern(args.reference)
    result: dict = dict(
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        input_revision=REVISION,
        reference_revision=REFERENCE,
        extension_sha256=digest(Path(clifft.__file__).parent / "_clifft_core.abi3.so"),
        source_sha256={
            name: digest(Path(__file__).with_name(name))
            for name in (
                "audit_cultivation_applicability.py",
                "parity_envelope_audit.py",
                "sahay_gadget.py",
                "study_sahay_gadget.py",
                "study_wan_scaling.py",
                "wan_corpus.py",
            )
        },
        external={},
        controls={},
    )
    names = [*CASES, "soft-d3", "soft-d5"]
    for name in names:
        model = None
        if name in CASES:
            path = args.source / "data/alternatives/inputs" / name / "model.json"
            model = json.loads(path.read_text())
            text = (corpus / (name + ".stim")).read_text()
        else:
            path = (
                args.source
                / "data/inputs"
                / ("clifft_d3_p001.stim" if name == "soft-d3" else "soft_cultivation_d5_p0005.stim")
            )
            text = (
                path.read_text() if name == "soft-d3" else (corpus / (name + ".stim")).read_text()
            )
        recognized, reason = pattern.recognize(text, clifft.parse)
        row: dict = dict(
            source_path=str(path.relative_to(args.source)),
            source_sha256=digest(path),
            pattern_recognizer=dict(selected=recognized is not None, reason=reason),
            inverse_t_envelopes=audit_parity(clifft, text) if not name.startswith("fold") else None,
        )
        if model is not None:
            row["literal_folded_model"], _ = folded_check(model, reference)
        checks = (
            [e for e in model["events"] if e["name"] == "MEASURE_PRODUCT"]
            if model is not None and "logical_x" in model
            else css_from_text(text)
        )
        row["terminal_css"] = css_summary(checks) if checks else None
        if row["inverse_t_envelopes"] is not None:
            row["parity_entry"] = parity_entry(clifft, text, checks, row["inverse_t_envelopes"])
        if name == "fold-d3":
            assert model is not None
            rewritten, change = canonical_fold_d3(model)
            row["canonicalization"] = change
            row["canonicalization"]["validation"] = check_rewrite(clifft, model, rewritten)
            row["canonical_folded_model"], gadget = folded_check(rewritten, reference)
            if gadget is None:
                raise AssertionError("Fold d3 canonicalization did not satisfy the certificate")
            row["canonical_entry"] = entry_validation(clifft, gadget)
            row["terminal_css"] = css_summary(gadget.checks)
        elif name == "fold-d5":
            assert model is not None
            _, gadget = folded_check(model, reference)
            row["canonical_entry"] = entry_validation(clifft, gadget)
            result["representation_controls"] = representation_checks(model, reference)
        if name == "chan-v1-d3":
            result["negative_parity_checks"] = negative_checks(clifft, text)
        result["external"][name] = row
        print(name, row["literal_folded_model"] if model else "no model schema", flush=True)
    for distance in (3, 5, 7, 9):
        model = make_model(reference.Surface(distance))
        result["controls"][f"constructed-d{distance}"], _ = folded_check(model, reference)
    for distance in (3, 5, 7):
        relative = f"tests/fixtures/folded/f{distance}.stim"
        text = subprocess.check_output(["git", "show", REFERENCE + ":" + relative], text=True)
        recognized, reason = pattern.recognize(text, clifft.parse)
        result["controls"][f"reconstructed-f{distance}"] = dict(
            selected=recognized is not None,
            reason=reason,
            scope="template and region certificate only; "
            "entry/growth validation is in prior studies",
        )
    result["export_manifest_sha256"] = digest(corpus / "manifest.json")
    result["exported_cases"] = list(manifest["cases"])
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
