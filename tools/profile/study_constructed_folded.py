"""Validate and time explicitly constructed d3/d5/d7/d9 terminal gadgets."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import statistics
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
import quadratic_oracle
from constructed_folded import make_model
from sahay_gadget import Boundary, SahayGadget
from sahay_kernel import SahayKernel
from study_sahay_gadget import compile_oracle, draw_region, oracle_circuit, random_faults
from study_wan_scaling import REFERENCE, family, load_reference, paired_masks
from wan_corpus import digest


def logical_states() -> list:
    return [
        np.array([1, 0], dtype=complex),
        np.array([0, 1], dtype=complex),
        np.array([1, 1]) / np.sqrt(2),
        np.array([1, 1j]) / np.sqrt(2),
        np.array([1, np.exp(1j * np.pi / 4)]) / np.sqrt(2),
        np.array([np.sqrt(0.3), np.sqrt(0.7) * np.exp(0.37j)]),
    ]


def cases(gadget: SahayGadget) -> list:
    rng = np.random.default_rng(270926 + gadget.n)
    result = []
    for i, logical in enumerate(logical_states()):
        offset = 0
        for dual in gadget.z_duals[:-1]:
            if i >= 2 and rng.integers(2):
                offset ^= dual
        signs = int.from_bytes(rng.bytes((gadget.rank + 7) // 8), "little")
        signs &= (1 << gadget.rank) - 1
        boundary = Boundary(
            offset, signs if i >= 2 else 0, logical, int(rng.integers(1 << len(gadget.ancillas)))
        )
        faults = random_faults(gadget, rng, 0.04) if i >= 2 else {}
        if i >= 4:
            triples = [j for j in gadget.noise if len(gadget.model["channels"][j]["support"]) == 3]
            faults[triples[i % len(triples)]] = 63
        result.append((boundary, faults))
    return result


def ideal_check(gadget: SahayGadget, kernel: SahayKernel) -> dict:
    syndrome = np.zeros(gadget.rank, dtype=np.uint8)
    h = np.array([[0, np.exp(-1j * np.pi / 4)], [np.exp(1j * np.pi / 4), 0]])
    error = 0.0
    for state in logical_states():
        for outcomes in itertools.product((0, 1), repeat=2):
            bound = gadget.bind(Boundary(0, 0, state, 0), {}, outcomes)
            actual = kernel.amplitudes(bound, syndrome)
            expected = (np.eye(2) + (-1) ** outcomes[0] * h) @ state / 2
            if outcomes[0] != outcomes[1]:
                expected *= 0
            error = max(
                error,
                float(
                    np.max(
                        np.abs(
                            np.outer(actual, actual.conj()) - np.outer(expected, expected.conj())
                        )
                    )
                ),
                abs(kernel.marginal(bound, syndrome, 0) - float(np.vdot(expected, expected).real)),
            )
    if error > 1e-10:
        raise AssertionError(("Ideal folded logical projector mismatch", error))
    return dict(states=6, outcome_pairs=4, maximum_absolute_error=error)


def gauss_check(gadget: SahayGadget, kernel: SahayKernel) -> dict:
    rng = np.random.default_rng(9027)
    error, comparisons = 0.0, 0
    relative_error, relative_checks, minimum_probability = 0.0, 0, 1.0
    for boundary, faults in cases(gadget):
        sample = draw_region(gadget, kernel, boundary, faults, rng)
        start = len(gadget.record_ids)
        syndrome = np.array(sample["records"][start : start + gadget.rank], dtype=np.uint8)
        for outcomes in itertools.product((0, 1), repeat=2):
            bound = gadget.bind(boundary, faults, outcomes)
            local = bound.values.copy()
            local[:, :, gadget.characters :, 1] *= 1 - 2 * syndrome.astype(int)
            expected_amp = np.array(
                [
                    sum(
                        bound.coefficients[q, i]
                        * quadratic_oracle.mean(gadget.masks, local[q, i], gadget.rank)
                        for i in range(4)
                    )
                    for q in range(2)
                ]
            )
            error = max(
                error, float(np.max(np.abs(kernel.amplitudes(bound, syndrome) - expected_amp)))
            )
            comparisons += 2
            for measured in sorted({0, 1, gadget.rank // 2, gadget.rank}):
                masks = paired_masks(gadget.masks, gadget.rank, measured)
                expected, absolute_sum = 0.0, 0.0
                for q in range(2):
                    for i in range(4):
                        for j in range(i + 1):
                            term = (
                                (1 if i == j else 2)
                                * bound.coefficients[q, i]
                                * bound.coefficients[q, j].conjugate()
                                * quadratic_oracle.mean(
                                    masks,
                                    np.concatenate((local[q, i], local[q, j].conjugate())),
                                    gadget.rank + measured,
                                )
                            ).real
                            expected += term
                            absolute_sum += abs(term)
                difference = abs(kernel.marginal(bound, syndrome, measured) - expected)
                error = max(error, difference)
                # Absolute tolerances alone cannot validate exponentially rare
                # complete records. Also check relative error away from exact
                # destructive cancellation of coherent terms.
                if abs(expected) > 1e-6 * absolute_sum and absolute_sum > 0:
                    relative_error = max(relative_error, difference / abs(expected))
                    relative_checks += 1
                    if expected > 0:
                        minimum_probability = min(minimum_probability, expected)
                comparisons += 1
    if error > 1e-10 or relative_error > 1e-8:
        raise AssertionError(("Independent Gauss sums disagree", error, relative_error))
    return dict(
        cases=6,
        comparisons=comparisons,
        maximum_absolute_error=float(error),
        noncancelling_relative_checks=relative_checks,
        maximum_relative_error=float(relative_error),
        minimum_positive_probability=float(minimum_probability),
    )


def physical_check(clifft: Any, gadget: SahayGadget, kernel: SahayKernel) -> dict:
    rng = np.random.default_rng(491)
    maximum, comparisons, width = 0.0, 0, 0
    for boundary, faults in cases(gadget):
        selected = [draw_region(gadget, kernel, boundary, faults, rng) for _ in range(2)]
        for axis in ("X", "Y", "Z", "magic"):
            program = compile_oracle(clifft, oracle_circuit(gadget, boundary, faults, axis))
            width = max(width, program.peak_active_width)
            records, expected = [], []
            for row in selected:
                bloch = row["bloch"]
                value = (
                    (bloch[0] + bloch[1]) / np.sqrt(2)
                    if axis == "magic"
                    else bloch["XYZ".index(axis)]
                )
                for bit in (0, 1):
                    records.append(row["records"][:-1] + [bit])
                    expected.append(row["norm"] * (1 + (-1) ** bit * value) / 2)
            actual = clifft.record_probabilities(program, np.array(records, dtype=np.uint8))
            maximum = max(maximum, float(np.max(np.abs(actual - expected))))
            comparisons += len(records)
    if maximum > 1e-10:
        raise AssertionError(("Physical Clifft record probabilities disagree", maximum))
    return dict(
        cases=6, comparisons=comparisons, peak_active_width=width, maximum_absolute_error=maximum
    )


def benchmark(gadget: SahayGadget, kernel: SahayKernel, shots: int) -> dict:
    rng = np.random.default_rng(4950927)
    boundary = Boundary(0, 0, logical_states()[4], 0)
    seconds, record_hashes, nonzero_records = [], [], []
    draw_region(gadget, kernel, boundary, {}, rng)
    for _ in range(3):
        digest = hashlib.sha256()
        nonzero = 0
        start = time.perf_counter()
        for _ in range(shots):
            row = draw_region(gadget, kernel, boundary, random_faults(gadget, rng, 0.001), rng)
            digest.update(bytes(row["records"]))
            nonzero += any(row["records"])
        seconds.append((time.perf_counter() - start) / shots)
        record_hashes.append(digest.hexdigest())
        nonzero_records.append(nonzero)
    return dict(
        contract="entire constructed terminal gadget from a supplied encoded boundary",
        probability=0.001,
        shots_per_trial=shots,
        seconds_per_shot=seconds,
        median_seconds_per_shot=statistics.median(seconds),
        record_hashes=record_hashes,
        attempts_with_nonzero_records=nonzero_records,
        early_rejection=False,
    )


def main() -> None:
    import clifft

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--distances", type=int, nargs="+", default=[3, 5, 7, 9])
    parser.add_argument("--shots", type=int, default=8)
    args = parser.parse_args()
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    reference = load_reference(args.reference)
    report: dict = dict(
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        reference_revision=REFERENCE,
        source_sha256={
            name: digest(Path(__file__).with_name(name))
            for name in (
                "study_constructed_folded.py",
                "constructed_folded.py",
                "quadratic_oracle.py",
                "sahay_gadget.py",
                "sahay_kernel.py",
                "sahay_contraction_kernel.cc",
                "study_sahay_gadget.py",
                "study_wan_scaling.py",
                "benchmark_wan_scaling.py",
            )
        },
        cpu_affinity=cpu,
        extension_sha256=digest(Path(clifft.__file__).parent / "_clifft_core.abi3.so"),
        quadratic_oracle_validation=quadratic_oracle.validate(),
        cases=[],
    )
    for distance in args.distances:
        model = make_model(reference.Surface(distance))
        gadget = SahayGadget(model, reference)
        plans = family(gadget.masks, gadget.rank)
        start = time.perf_counter()
        kernel = SahayKernel(gadget, args.reference, args.work / f"d{distance}")
        setup = time.perf_counter() - start
        (args.work / f"d{distance}" / "model.json").write_text(json.dumps(model, indent=2) + "\n")
        row = dict(
            distance=distance,
            provenance=model["provenance"],
            model_sha256=digest(args.work / f"d{distance}" / "model.json"),
            certificate=gadget.certificate,
            records=len(gadget.record_ids) + 2 * gadget.rank + 1,
            syndrome_order=gadget.order,
            syndrome_order_policy="amplitude min-scope elimination order",
            peak_contraction_entries=plans["peak_entries"],
            numeric_payload_bytes=plans["numeric_payload_bytes"],
            setup_seconds=setup,
        )
        print(f"d{distance}: plans built; checking ideal operator", flush=True)
        row["ideal_projector"] = ideal_check(gadget, kernel)
        print(f"d{distance}: checking independent quadratic sums", flush=True)
        row["independent_contraction_check"] = gauss_check(gadget, kernel)
        if distance <= 5:
            print(f"d{distance}: checking physical Clifft probabilities", flush=True)
            row["physical_reference"] = physical_check(clifft, gadget, kernel)
        row["benchmark"] = benchmark(gadget, kernel, args.shots)
        report["cases"].append(row)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        kernel.close()
        milliseconds = row["benchmark"]["median_seconds_per_shot"] * 1000
        print(
            f"d{distance}: {milliseconds:.3f} ms per terminal gadget",
            flush=True,
        )


if __name__ == "__main__":
    main()
