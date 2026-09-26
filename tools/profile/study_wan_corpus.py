"""Bounded baseline, fault-witness and eligibility study for the Wan corpus."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import os
import platform
import re
import resource
import statistics
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np
from wan_corpus import REVISION, Exporter, digest, pauli_text, verify_source


def fresh_reset_wires(text: str) -> str:
    """Purify resets with fresh zero wires for the visible-record probability API.

    Old wires are never touched again and are traced by the probability query.
    No measurement is added or removed. This is only a noiseless oracle input;
    timed sampling always uses the original resets and physical wire count.
    """
    mapping: dict[int, int] = {}
    next_qubit = 0
    lines = []

    def allocate(q):
        nonlocal next_qubit
        mapping[q] = next_qubit
        next_qubit += 1
        return mapping[q]

    def target(match):
        q = int(match[0])
        return str(mapping[q] if q in mapping else allocate(q))

    for line in text.splitlines():
        gate, targets = line.split(" ", 1)
        if gate.startswith(("DETECTOR", "OBSERVABLE_INCLUDE")):
            continue
        if gate in {"R", "RX"}:
            for q in map(int, targets.split()):
                fresh = allocate(q)
                if gate == "RX":
                    lines.append(f"H {fresh}")
        else:
            tokens = [
                t if t.startswith("rec[") else re.sub(r"\d+", target, t) for t in targets.split()
            ]
            lines.append(gate + " " + " ".join(tokens))
    return "\n".join(lines) + "\n"


def xor_mask(indices: list[int]) -> int:
    result = 0
    for i in indices:
        result ^= 1 << i
    return result


def nullspace(rows: list[int], n: int) -> list[int]:
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    basis = []
    for free in range(n):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot, row in sorted(pivots.items()):
            if (vector & row).bit_count() % 2:
                vector ^= 1 << pivot
        basis.append(vector)
    return basis


def witness(clifft: Any, root: Path, case: str) -> dict:
    model = json.loads((root / "data/alternatives/inputs" / case / "model.json").read_text())
    candidate = json.loads(
        (root / "data/alternatives/certificates" / case / "fault_pattern.json").read_text()
    )
    locations = candidate["growth_faults"] + [candidate["final_fault"]]
    for item in locations:
        channel = model["channels"][item["source_location"]]
        if channel["event_index"] != item["event_index"] or channel["support"] != item["support"]:
            raise ValueError("Fault witness does not match its model")
    faults = {v["source_location"]: v["outcome"]["code"] for v in locations}
    report = dict(faults=faults)
    for label, selected in (("reference", {}), ("four_faults", faults)):
        text, _ = Exporter(model, 0, selected).export()
        program = clifft.compile(fresh_reset_wires(text))
        basis = nullspace([xor_mask(d) for d in model["detectors"]], program.num_measurements)
        if len(basis) > 12:
            raise ValueError("Witness record space exceeds the bounded oracle")
        words = [0]
        for vector in basis:
            words += [word ^ vector for word in words]
        records = np.array(
            [[(word >> i) & 1 for i in range(program.num_measurements)] for word in words],
            dtype=np.uint8,
        )
        start = time.perf_counter()
        probabilities = clifft.record_probabilities(program, records)
        bad = records[:, model["observables"][0]].sum(axis=1) % 2
        acceptance = float(probabilities.sum())
        error = float(probabilities[bad == 1].sum())
        expected = 1 if label == "reference" else (0.25 if case == "chan-v1-d5" else 0)
        expected_error = 0 if label == "reference" else expected
        if abs(acceptance - expected) > 1e-12 or abs(error - expected_error) > 1e-12:
            raise AssertionError((case, label, acceptance, error))
        report[label] = dict(
            A=acceptance,
            B=error,
            expected_A=expected,
            expected_B=expected_error,
            accepted_record_dimension=len(basis),
            records=len(records),
            peak_active_width=program.peak_active_width,
            seconds=time.perf_counter() - start,
            probabilities=probabilities.tolist(),
        )
    return report


def summarize(result: Any, *, survivors: bool = False) -> dict:
    shots = int(result.total_shots) if survivors else len(result.measurements)
    if survivors:
        accepted = int(result.passed_shots)
        errors = int(np.any(result.observables, axis=1).sum())
    else:
        keep = ~np.any(result.detectors, axis=1)
        accepted = int(keep.sum())
        errors = int(np.any(result.observables[keep], axis=1).sum())
    return dict(
        attempts=shots,
        accepted=accepted,
        errors=errors,
        A=accepted / shots,
        B=errors / shots,
        LER=errors / accepted if accepted else None,
    )


def one_fault_acceptance(clifft: Any, corpus: Path, manifest: dict) -> dict:
    """Check the first acceptance coefficient against physical one-fault sampling."""
    report = {}
    for case, meta in manifest["cases"].items():
        if "reference" not in meta:
            continue
        text = (corpus / (case + ".stim")).read_text()
        ref = meta["reference"]
        physical = len(meta["noise_locations"])
        retained = ref["locations_in_coefficient_convention"]
        # The references sometimes sum identity-response sites out of the series.
        # Restore their degree-one mass before conditioning on physical fault count.
        expected = float((Fraction(ref["a"][1]) + physical - retained) / physical)
        program = clifft.compile(
            text,
            expected_detectors=meta["expected_detectors"],
            expected_observables=meta["expected_observables"],
            postselection_mask=[1] * len(meta["expected_detectors"]),
        )
        if len(program.noise_site_probabilities) != physical:
            raise ValueError("Compiler and model count different physical noise sites")
        np.testing.assert_allclose(program.noise_site_probabilities, manifest["p"], rtol=1e-12)
        shots = 128 if case == "fold-d5" else (2048 if case == "rp2-d5" else 65536)
        start = time.perf_counter()
        result = clifft.sample_k_survivors(
            program, shots, k=1, seed=49501, threads=1, batch_size=1, keep_records=True
        )
        counts = summarize(result, survivors=True)
        sigma = math.sqrt(expected * (1 - expected) / shots)
        if abs(counts["A"] - expected) > 6 * sigma or counts["errors"]:
            raise AssertionError((case, counts, expected, sigma))
        report[case] = dict(
            **counts,
            expected_acceptance=expected,
            standard_error=sigma,
            z=(counts["A"] - expected) / sigma,
            seconds=time.perf_counter() - start,
            physical_locations=physical,
            reference_locations=retained,
        )
        print(case, "one-fault acceptance", counts["A"], "expected", expected, flush=True)
    return report


def folded_boundary(clifft: Any, root: Path) -> dict:
    """Inspect the data-code invariant before the first large Sahay check."""
    model = json.loads((root / "data/alternatives/inputs/fold-d5/model.json").read_text())
    end = next(i for i, event in enumerate(model["events"]) if event.get("native_block") == 2)
    exporter = Exporter(model, 0.001)
    for i, event in enumerate(model["events"][:end]):
        exporter.dictionary_event(event, i)
    checks = [event for event in model["events"] if event["name"] == "MEASURE_PRODUCT"]
    text = "\n".join(exporter.lines) + "\n"
    text += "".join(
        "EXP_VAL " + pauli_text(dict.fromkeys(e["support"], e["axis"])) + "\n" for e in checks
    )
    program = clifft.compile(text)
    probes = [line for line in program.inspect().splitlines() if "WRITE_EXPECTATION" in line]
    if len(probes) != 40 or any(" w1 I sign=" not in line for line in probes):
        raise AssertionError("The compiled prefix no longer certifies the signed data code")
    result = clifft.sample(program, 1024, seed=495, threads=1, batch_size=1)
    deviation = float(np.max(np.abs(np.abs(result.exp_vals) - 1)))
    if deviation > 1e-12:
        raise AssertionError("Sampled prefix leaves a code-check eigenspace")
    return dict(
        prefix_event_stop_exclusive=end,
        prefix_peak_active_width=program.peak_active_width,
        code_check_count=len(probes),
        code_probes=probes,
        sampled_shots=1024,
        max_abs_expectation_deviation_from_one=deviation,
        native_block_gate_counts={
            gate: sum(e.get("native_block") == 2 and e["name"] == gate for e in model["events"])
            for gate in ("CSX", "CS_DAG_X", "CCZ")
        },
        scope="Prefix includes cat preparation; certifies the signed data code only, "
        "not classical ancillas or a complete specialized region.",
    )


def benchmark(
    clifft: Any, text: str, meta: dict, *, shots: int, p: float, scheduled: bool = False
) -> dict:
    rows = {}
    for mode in ("full", "postselected"):
        postselect = mode == "postselected"
        start = time.perf_counter()
        circuit = clifft.parse(text)
        detector_count = sum(node.gate.name == "DETECTOR" for node in circuit.nodes)
        options = {}
        scheduler = None
        if scheduled:
            passes = clifft.default_hir_pass_manager()
            scheduler = clifft.ActiveWidthSchedulePass(search_budget=None)
            passes.add(scheduler)
            options["hir_passes"] = passes
        program = clifft.compile(
            text,
            expected_detectors=meta["expected_detectors"],
            expected_observables=meta["expected_observables"],
            postselection_mask=[1] * detector_count if postselect else None,
            **options,
        )
        compile_seconds = time.perf_counter() - start
        if program.peak_active_width > 24:
            raise ValueError("Refusing a dense worker above width 24")
        sampler = clifft.sample_survivors if postselect else clifft.sample
        kwargs = dict(threads=1, batch_size=1)
        if postselect:
            kwargs["keep_records"] = True
        start = time.perf_counter()
        sampler(program, 8, seed=495, **kwargs)
        probe_seconds = time.perf_counter() - start
        # Keep every case bounded while collecting repeated timings of useful size.
        batch = max(8, min(shots, int(0.3 * 8 / max(probe_seconds, 1e-6))))
        timings, counts = [], []
        for i in range(3):
            start = time.perf_counter()
            result = sampler(program, batch, seed=49500 + i, **kwargs)
            timings.append(time.perf_counter() - start)
            counts.append(summarize(result, survivors=postselect))
        rows[mode] = dict(
            p=p,
            compile_seconds=compile_seconds,
            peak_active_width=program.peak_active_width,
            qubits=program.num_qubits,
            measurements=program.num_measurements,
            detectors=program.num_detectors,
            noise_sites=len(program.noise_site_probabilities),
            batch_shots=batch,
            sample_seconds=timings,
            counts=counts,
            median_us_per_attempt=statistics.median(timings) * 1e6 / batch,
        )
        if scheduler is not None:
            rows[mode]["scheduler"] = dict(
                applied=scheduler.applied,
                incumbent_peak=scheduler.incumbent_peak,
                result_peak=scheduler.result_peak,
                incumbent_dense_work=scheduler.incumbent_dense_work,
                result_dense_work=scheduler.result_dense_work,
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--shots", type=int, default=32768)
    parser.add_argument("--eligibility-only", action="store_true")
    parser.add_argument("--strata-only", action="store_true")
    parser.add_argument("--schedule", action="store_true")
    parser.add_argument("--boundary-only", action="store_true")
    parser.add_argument("--control-circuit", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.package.resolve()))
    import clifft

    if Path(clifft.__file__).resolve().parent != args.package.resolve() / "clifft":
        raise ValueError("A different Clifft package shadowed the requested build")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    verify_source(args.source)
    manifest = json.loads((args.corpus / "manifest.json").read_text())
    report = dict(
        upstream_revision=REVISION,
        corpus_p=manifest["p"],
        python=sys.version,
        platform=platform.platform(),
        package=str(args.package),
        extension_sha256=digest(next((args.package / "clifft").glob("_clifft_core*.so"))),
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        active_width_schedule=args.schedule,
        cases={},
    )

    def save():
        report["peak_process_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    if args.strata_only:
        report["one_fault_acceptance"] = one_fault_acceptance(clifft, args.corpus, manifest)
        save()
        return
    if args.boundary_only:
        report["folded_boundary"] = folded_boundary(clifft, args.source)
        save()
        return
    if args.eligibility_only and args.control_circuit:
        control = getattr(clifft, "compile")(
            args.control_circuit.read_text(), specialize_folded=True
        )
        if not control.has_folded_regions:
            raise AssertionError("The positive control did not select specialization")
        report["positive_control"] = dict(
            circuit_sha256=digest(args.control_circuit),
            selected=True,
            peak_active_width=control.peak_active_width,
        )

    for name, meta in manifest["cases"].items():
        text = (args.corpus / (name + ".stim")).read_text()
        if hashlib.sha256(text.encode()).hexdigest() != meta["circuit_sha256"]:
            raise ValueError("Generated circuit digest mismatch")
        print(name, flush=True)
        row = dict(circuit_sha256=meta["circuit_sha256"], model_sha256=meta.get("model_sha256"))
        report["cases"][name] = row
        if args.eligibility_only:
            recognize = importlib.import_module("clifft._folded").recognize
            matched, reason = recognize(text, clifft.parse)
            row.update(selected=matched is not None, reason=reason)
            fallback = getattr(clifft, "compile")(text, specialize_folded=True)
            row.update(
                public_api_selected=fallback.has_folded_regions,
                fallback_peak_active_width=fallback.peak_active_width,
            )
            save()
            continue
        row["noisy"] = benchmark(
            clifft, text, meta, shots=args.shots, p=manifest["p"], scheduled=args.schedule
        )
        if "model_path" in meta:
            model = json.loads((args.source / meta["model_path"]).read_text())
            zero_text, zero_meta = Exporter(model, 0).export()
        else:
            zero_text = re.sub(
                r"((?:DEPOLARIZE[123]|[XYZ]_ERROR|M[XYZ]?|MR[XYZ]?)\()[^)]*(\))", r"\g<1>0\2", text
            )
            zero_meta = meta
        row["noiseless"] = benchmark(
            clifft, zero_text, zero_meta, shots=min(1024, args.shots), p=0, scheduled=args.schedule
        )
        for mode in row["noiseless"].values():
            if any(c["accepted"] != c["attempts"] or c["errors"] for c in mode["counts"]):
                raise AssertionError("Nontrivial fault-free output: " + name)
        save()
    if not args.eligibility_only:
        report["fixed_faults"] = {
            case: witness(clifft, args.source, case)
            for case in ("chan-v1-d5", "chan-four-round-d5")
        }
        # This is the sole matching all-orders acceptance reference in this corpus.
        if manifest["p"] == 0.001:
            soft = report["cases"]["soft-d5"]["noisy"]["postselected"]
            counts = soft["counts"]
            n = sum(c["attempts"] for c in counts)
            a = sum(c["accepted"] for c in counts) / n
            expected = 0.143955584096
            sigma = math.sqrt(expected * (1 - expected) / n)
            report["soft_acceptance"] = dict(
                A=a, expected_A=expected, attempts=n, standard_error=sigma, z=(a - expected) / sigma
            )
            if abs(a - expected) > 6 * sigma:
                raise AssertionError("SOFT acceptance differs from the exact reference")
    save()


if __name__ == "__main__":
    main()
