"""Validate and measure an external folded-code instrument with full records."""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import stim
from sahay_gadget import Bound, Boundary, SahayGadget, parity
from sahay_kernel import SahayKernel
from study_wan_corpus import fresh_reset_wires
from study_wan_scaling import REFERENCE, load_reference
from wan_corpus import Exporter, digest, pauli_text, verify_source


def physical_pauli(gadget: SahayGadget, mask: int, axis: str) -> str:
    return pauli_text({q: axis for i, q in enumerate(gadget.data) if mask >> i & 1})


def preparation(gadget: SahayGadget, boundary: Boundary) -> str:
    stabilizers = []
    for rows, axis in ((gadget.x_rows, "X"), (gadget.z_rows, "Z")):
        for i, row in enumerate(rows):
            op = stim.PauliString(gadget.n)
            for q in range(gadget.n):
                if row >> q & 1:
                    op[q] = axis
            sign = (boundary.signs >> i) & 1 if axis == "X" else parity(boundary.offset & row)
            op.sign = (-1) ** sign
            stabilizers.append(op)
    z = stim.PauliString(gadget.n)
    for q in range(gadget.n):
        if gadget.logical_z >> q & 1:
            z[q] = "Z"
    if parity(boundary.offset & gadget.logical_z):
        raise ValueError("Oracle preparation requires a positive logical-Z base")
    stabilizers.append(z)
    encoding = stim.Tableau.from_stabilizers(stabilizers).to_circuit()
    lines = [
        op.name + " " + " ".join(str(gadget.data[t.value]) for t in op.targets_copy())
        for op in encoding
    ]
    # Physical logical rotations avoid depending on the encoder's choice of X
    # destabilizer and therefore exercise relative phases independently.
    x = {q: "X" for i, q in enumerate(gadget.data) if gadget.logical_x >> i & 1}
    z_axes = {q: "Z" for i, q in enumerate(gadget.data) if gadget.logical_z >> i & 1}
    y = {
        q: "Y" if q in x and q in z_axes else "X" if q in x else "Z"
        for q in x.keys() | z_axes.keys()
    }
    alpha, beta = boundary.logical
    theta = 2 * math.atan2(abs(beta), abs(alpha))
    phi = float(np.angle(beta) - np.angle(alpha)) if abs(alpha * beta) > 1e-12 else 0.0
    lines += [
        f"R_PAULI({theta / np.pi:.17g}) " + pauli_text(y),
        f"R_PAULI({phi / np.pi:.17g}) " + pauli_text(z_axes),
    ]
    lines += [f"X {q}" for i, q in enumerate(gadget.ancillas) if boundary.ancillas >> i & 1]
    return "\n".join(lines) + "\n"


def oracle_circuit(
    gadget: SahayGadget, boundary: Boundary, faults: dict[int, int], axis: str = "magic"
) -> str:
    exporter = Exporter(gadget.model, 0, faults)
    for i in range(gadget.start, len(gadget.model["events"])):
        exporter.dictionary_event(gadget.model["events"][i], i)
    lines = [line for line in exporter.lines if not line.startswith("DETECTOR")]
    if axis == "magic":
        lines.append("R_PAULI(-0.25) " + physical_pauli(gadget, gadget.logical_z, "Z"))
        lines.append("MPP " + physical_pauli(gadget, gadget.logical_x, "X"))
    else:
        x = {q: "X" for i, q in enumerate(gadget.data) if gadget.logical_x >> i & 1}
        z = {q: "Z" for i, q in enumerate(gadget.data) if gadget.logical_z >> i & 1}
        y = {q: "Y" if q in x and q in z else "X" if q in x else "Z" for q in x.keys() | z.keys()}
        lines.append("MPP " + pauli_text({"X": x, "Y": y, "Z": z}[axis]))
    return preparation(gadget, boundary) + "\n".join(lines) + "\n"


def record_for(gadget: SahayGadget, bound: Bound, syndrome: Any) -> list[int]:
    records = list(bound.records)
    indices = {"X": 0, "Z": 0}
    for event in gadget.checks:
        axis = event["axis"]
        i = indices[axis]
        records.append(int(syndrome[i]) if axis == "X" else parity(bound.offset & gadget.z_rows[i]))
        indices[axis] += 1
    return records


def draw_syndrome(
    gadget: SahayGadget, kernel: SahayKernel, bound: Bound, rng: Any, mass: float
) -> tuple:
    syndrome = np.zeros(gadget.rank, dtype=np.uint8)
    for k, original in enumerate(gadget.order):
        zero = kernel.marginal(bound, syndrome, k + 1)
        if zero > mass * (1 + 1e-8):
            raise AssertionError("Prefix probability increased")
        zero = min(zero, mass)
        if rng.random() * mass >= zero:
            syndrome[original] = 1
            mass -= zero
        else:
            mass = zero
    amplitudes = kernel.amplitudes(bound, syndrome)
    norm = float(np.vdot(amplitudes, amplitudes).real)
    # Individual complete records can be exponentially rare at larger distance.
    # Their absolute probability is not a validity threshold.
    if norm <= 0 or abs(norm - mass) > 1e-8 * max(norm, mass):
        raise AssertionError(("Full marginal and output amplitudes disagree", norm, mass))
    return syndrome, amplitudes, norm


def draw_region(
    gadget: SahayGadget, kernel: SahayKernel, boundary: Boundary, faults: dict[int, int], rng: Any
) -> dict:
    syndrome = np.zeros(gadget.rank, dtype=np.uint8)
    bounds = [
        gadget.bind(boundary, faults, outcomes) for outcomes in itertools.product((0, 1), repeat=2)
    ]
    weights = np.array([kernel.marginal(bound, syndrome, 0) for bound in bounds])
    if abs(weights.sum() - 1) > 1e-9:
        raise AssertionError("Body outcomes are not normalized")
    target = rng.random() * weights.sum()
    choice = min(int(np.searchsorted(np.cumsum(weights), target, side="right")), 3)
    bound = bounds[choice]
    syndrome, amplitudes, norm = draw_syndrome(gadget, kernel, bound, rng, weights[choice])
    bloch = gadget.logical_bloch(bound, amplitudes)
    p0 = (1 + (bloch[0] + bloch[1]) / np.sqrt(2)) / 2
    bit = int(rng.random() >= p0)
    return dict(
        records=record_for(gadget, bound, syndrome) + [bit],
        norm=norm,
        amplitudes=amplitudes,
        bloch=bloch,
        offset=bound.offset,
        record_probability=norm * (p0 if bit == 0 else 1 - p0),
    )


def random_faults(gadget: SahayGadget, rng: Any, p: float) -> dict[int, int]:
    faults = {}
    for location in gadget.noise:
        if rng.random() < p:
            channel = gadget.model["channels"][location]
            faults[location] = (
                1
                if channel["kind"] == "flip"
                else int(rng.integers(1, 4 ** len(channel["support"])))
            )
    return faults


def compile_oracle(clifft: Any, text: str) -> Any:
    passes = clifft.default_hir_pass_manager()
    passes.add(clifft.ActiveWidthSchedulePass(search_budget=None))
    hir = clifft.trace(clifft.parse(fresh_reset_wires(text)))
    passes.run(hir)
    program = clifft.lower(hir)
    if program.peak_active_width > 24:
        raise ValueError("Refusing an oversized dense oracle")
    return program


def validate(clifft: Any, gadget: SahayGadget, kernel: SahayKernel) -> dict:
    rng = np.random.default_rng(51026)
    cases = []
    logicals = [
        np.array([1, 0], dtype=complex),
        np.array([0, 1], dtype=complex),
        np.array([1, 1]) / np.sqrt(2),
        np.array([1, 1j]) / np.sqrt(2),
        np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2),
        np.array([np.sqrt(0.3), np.sqrt(0.7) * np.exp(0.37j)]),
    ]
    maximum_error = 0.0
    channels = [gadget.model["channels"][location] for location in gadget.noise]
    pairs = [c for c in channels if len(c["support"]) == 2]
    triples = [c for c in channels if len(c["support"]) == 3]
    reported = [c for c in channels if c.get("role") == "source measurement reported-bit flip"]
    singles = [c for c in channels if c["kind"] == "depolarizing" and len(c["support"]) == 1]
    targeted = [
        {pairs[0]["location"]: 15},
        {triples[0]["location"]: 63},
        {triples[-1]["location"]: 21},
        {reported[3]["location"]: 1},
        {singles[0]["location"]: 3, singles[-1]["location"]: 2},
        {reported[3]["location"]: 1, reported[9]["location"]: 1, triples[0]["location"]: 63},
    ]
    for index in range(18):
        signs = 0 if index < 6 else int(rng.integers(1 << gadget.rank))
        offset = 0
        if index >= 6:
            for dual in gadget.z_duals[:-1]:
                if rng.integers(2):
                    offset ^= dual
        boundary = Boundary(
            offset,
            signs,
            logicals[index % 6],
            0 if index < 6 else int(rng.integers(1 << len(gadget.ancillas))),
        )
        faults = (
            {}
            if index < 6
            else random_faults(gadget, rng, 0.02)
            if index < 12
            else targeted[index - 12]
        )
        selected = [draw_region(gadget, kernel, boundary, faults, rng) for _ in range(4)]
        row: dict[str, Any] = dict(
            case=index,
            faults=faults,
            x_signs=signs,
            offset=offset,
            ancillas=boundary.ancillas,
            logical=[[float(a.real), float(a.imag)] for a in boundary.logical],
            sampled_records=[row["records"] for row in selected],
            axes={},
        )
        for axis in ("X", "Y", "Z", "magic"):
            text = oracle_circuit(gadget, boundary, faults, axis)
            program = compile_oracle(clifft, text)
            records, expected = [], []
            for selected_row in selected:
                base = selected_row["records"][:-1]
                bloch = selected_row["bloch"]
                expectation = (
                    (bloch[0] + bloch[1]) / np.sqrt(2)
                    if axis == "magic"
                    else bloch["XYZ".index(axis)]
                )
                for bit in (0, 1):
                    records.append(base + [bit])
                    expected.append(selected_row["norm"] * (1 + (-1) ** bit * expectation) / 2)
            actual = clifft.record_probabilities(program, np.array(records, dtype=np.uint8))
            error = float(np.max(np.abs(actual - expected)))
            maximum_error = max(maximum_error, error)
            if error > 2e-9:
                raise AssertionError((index, axis, error, actual, expected))
            row["axes"][axis] = dict(
                maximum_absolute_error=error,
                comparisons=len(records),
                peak_active_width=program.peak_active_width,
            )
        cases.append(row)
        print("validated", index, "faults", len(faults), "max error", maximum_error, flush=True)
    return dict(cases=cases, maximum_absolute_error=maximum_error)


def entry_validation(clifft: Any, gadget: SahayGadget) -> dict:
    text, _ = gadget.prefix(0.001)
    program = clifft.compile(text)
    probes = [line for line in program.inspect().splitlines() if "WRITE_EXPECTATION" in line]
    rank = gadget.rank
    classical = probes[: 2 * rank] + probes[2 * rank + 3 :]
    if any(" w1 I sign=" not in line for line in classical):
        raise AssertionError("Entry data checks or ancillas are not symbolic signs")
    if any(" w1 " not in line or " I " in line for line in probes[2 * rank : 2 * rank + 3]):
        raise AssertionError("Logical entry is not confined to one active coordinate")
    result = clifft.sample(program, 1024, seed=510, threads=1, batch_size=1)
    for values in result.exp_vals:
        gadget.boundary(values)
    return dict(
        shots=1024,
        probes=probes,
        prefix_measurements=program.num_measurements,
        prefix_peak_active_width=program.peak_active_width,
        maximum_purity_error=float(
            np.max(np.abs(np.sum(result.exp_vals[:, 2 * rank : 2 * rank + 3] ** 2, axis=1) - 1))
        ),
    )


def complete_record_validation(clifft: Any, gadget: SahayGadget, kernel: SahayKernel) -> dict:
    rng = np.random.default_rng(91027)
    cases = []
    for case in range(5):
        faults = {}
        if case:
            for channel in gadget.model["channels"]:
                if rng.random() < 0.001 * case:
                    faults[channel["location"]] = (
                        1
                        if channel["kind"] == "flip"
                        else int(rng.integers(1, 4 ** len(channel["support"])))
                    )
        prefix, _ = gadget.prefix(0, faults)
        prefix_program = clifft.compile(prefix)
        prefix_result = clifft.sample(prefix_program, 3, seed=810 + case, threads=1, batch_size=1)
        prefix_oracle = compile_oracle(
            clifft,
            "\n".join(line for line in prefix.splitlines() if not line.startswith("EXP_VAL")),
        )
        prefix_probabilities = clifft.record_probabilities(
            prefix_oracle, prefix_result.measurements
        )
        full_text, _ = Exporter(gadget.model, 0, faults).export()
        program = compile_oracle(clifft, full_text)
        records, expected = [], []
        for i, values in enumerate(prefix_result.exp_vals):
            sample = draw_region(gadget, kernel, gadget.boundary(values), faults, rng)
            records.append(prefix_result.measurements[i].tolist() + sample["records"])
            expected.append(prefix_probabilities[i] * sample["record_probability"])
        actual = clifft.record_probabilities(program, np.array(records, dtype=np.uint8))
        relative = float(np.max(np.abs(actual - expected) / np.maximum(expected, 1e-30)))
        if relative > 1e-8:
            raise AssertionError(("Full-record handoff differs", case, actual, expected))
        cases.append(
            dict(
                case=case,
                faults=faults,
                records=records,
                probabilities=actual.tolist(),
                expected=expected,
                maximum_relative_error=relative,
                peak_active_width=program.peak_active_width,
            )
        )
        print("full records", case, "faults", len(faults), "relative error", relative, flush=True)
    return dict(cases=cases)


def benchmark(clifft: Any, gadget: SahayGadget, kernel: SahayKernel) -> dict:
    full_text, meta = Exporter(gadget.model, 0.001).export()
    passes = clifft.default_hir_pass_manager()
    passes.add(clifft.ActiveWidthSchedulePass(search_budget=None))
    start = time.perf_counter()
    full = clifft.compile(
        full_text,
        hir_passes=passes,
        expected_detectors=meta["expected_detectors"],
        expected_observables=meta["expected_observables"],
    )
    full_compile_seconds = time.perf_counter() - start
    prefix_text, prefix_meta = gadget.prefix(0.001)
    start = time.perf_counter()
    prefix = clifft.compile(prefix_text, expected_detectors=prefix_meta["expected_detectors"])
    prefix_compile_seconds = time.perf_counter() - start
    detectors = [
        [meta["source_record_map"][r] for r in event["records"]]
        for event in gadget.model["events"]
        if event["name"] == "DETECTOR"
    ]
    if len(detectors) != full.num_detectors:
        raise AssertionError("Detector mapping differs from the physical model")
    shots = 64
    timings: dict[str, list] = dict(scheduled=[], gadget=[])

    def run(variant: str, count: int, seed: int) -> dict:
        start = time.perf_counter()
        if variant == "scheduled":
            result = clifft.sample(full, count, seed=seed, threads=1, batch_size=1)
            records = result.measurements
            flags = result.detectors
            bad = result.observables[:, 0]
        else:
            rng = np.random.default_rng(seed)
            result = clifft.sample(prefix, count, seed=seed, threads=1, batch_size=1)
            records = np.empty((count, full.num_measurements), dtype=np.uint8)
            for i, values in enumerate(result.exp_vals):
                sample = draw_region(
                    gadget, kernel, gadget.boundary(values), random_faults(gadget, rng, 0.001), rng
                )
                records[i] = result.measurements[i].tolist() + sample["records"]
            flags = np.array(
                [
                    np.bitwise_xor.reduce(records[:, slots], axis=1) ^ parity
                    for slots, parity in zip(detectors, meta["expected_detectors"], strict=True)
                ]
            ).T
            bad = records[:, -1]
        seconds = time.perf_counter() - start
        accepted = ~np.any(flags, axis=1)
        return dict(
            shots=count,
            seconds=seconds,
            accepted=int(accepted.sum()),
            errors=int((accepted & (bad != 0)).sum()),
            record_sha256=hashlib.sha256(records.tobytes()).hexdigest(),
        )

    for name in timings:
        run(name, 2, 519)
    for trial in range(5):
        for name in list(timings) if trial % 2 == 0 else list(reversed(timings)):
            result = run(name, shots, 91000 + trial)
            timings[name].append(result)
            print("benchmark", trial, name, result, flush=True)
    return dict(
        p=0.001,
        mode="full noisy attempts with all records retained",
        shots_per_trial=shots,
        prefix_noise_sites=len(gadget.model["channels"]) - len(gadget.noise),
        region_noise_sites=len(gadget.noise),
        scope="Python research binding plus native arithmetic; production specialization unchanged",
        measurements=full.num_measurements,
        detectors=full.num_detectors,
        compile_seconds=dict(scheduled=full_compile_seconds, prefix=prefix_compile_seconds),
        trials=timings,
        median_ms_per_attempt={
            name: statistics.median(row["seconds"] * 1000 / row["shots"] for row in rows)
            for name, rows in timings.items()
        },
    )


def negative_validation(gadget: SahayGadget, reference: Any) -> dict:
    cases = []
    events = gadget.model["events"]
    first = next(i for i in range(gadget.start, gadget.stop) if events[i]["name"] == "CSX")
    for label in (
        "coherent measurement",
        "data-controlled pair",
        "wrong code translation",
        "missing closing Hadamard",
    ):
        model = copy.deepcopy(gadget.model)
        if label == "coherent measurement":
            model["events"].insert(
                first,
                dict(
                    name="MEASURE",
                    support=[events[first]["support"][0]],
                    axis="Z",
                    inverted=False,
                    record=999,
                ),
            )
        elif label == "data-controlled pair":
            i = next(i for i in range(gadget.start, gadget.stop) if events[i]["name"] == "CCZ")
            model["events"][i]["support"][0] = gadget.data[0]
        elif label == "wrong code translation":
            model["events"][first]["support"][1] = gadget.data[1]
        else:
            i = max(i for i in range(gadget.start, gadget.stop) if events[i]["name"] == "H")
            del model["events"][i]
        try:
            SahayGadget(model, reference)
        except ValueError as error:
            cases.append(dict(mutation=label, rejection=str(error)))
        else:
            raise AssertionError("Unsound structure accepted: " + label)
    return dict(cases=cases)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {0})
    verify_source(args.source)
    sys.path.insert(0, str(args.package.resolve()))
    import clifft

    if Path(clifft.__file__).resolve().parent != args.package.resolve() / "clifft":
        raise ValueError("Selected package was shadowed")
    model_path = args.source / "data/alternatives/inputs/fold-d5/model.json"
    reference = load_reference(args.reference)
    start = time.perf_counter()
    gadget = SahayGadget(json.loads(model_path.read_text()), reference)
    kernel = SahayKernel(gadget, args.reference, args.work)
    setup_seconds = time.perf_counter() - start
    output: dict[str, Any] = dict(
        model_sha256=digest(model_path),
        certificate=gadget.certificate,
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        reference_revision=REFERENCE,
        extension_sha256=digest(next((args.package / "clifft").glob("*.so"))),
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        setup_seconds_including_native_build=setup_seconds,
        native_library_sha256=digest(args.work / "kernel.so"),
        source_sha256={
            name: digest(Path(__file__).with_name(name))
            for name in (
                "sahay_gadget.py",
                "sahay_kernel.py",
                "sahay_contraction_kernel.cc",
                "study_sahay_gadget.py",
            )
        },
    )
    print("certificate", gadget.certificate, flush=True)
    output["validation"] = validate(clifft, gadget, kernel)
    output["negative_validation"] = negative_validation(gadget, reference)
    output["entry_validation"] = entry_validation(clifft, gadget)
    output["complete_record_validation"] = complete_record_validation(clifft, gadget, kernel)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    output["benchmark"] = benchmark(clifft, gadget, kernel)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    kernel.close()


if __name__ == "__main__":
    main()
