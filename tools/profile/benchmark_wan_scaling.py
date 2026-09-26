"""Time fixed contraction plans after offline syndrome-order selection.

This benchmarks one arithmetic query per prefix, not full physical attempts.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import time
from pathlib import Path
from typing import Any

from study_wan_scaling import (
    REFERENCE,
    bits,
    family,
    load_reference,
    masks_for,
    paired_masks,
    permute_masks,
    sahay_masks,
    scope_plan,
)
from wan_corpus import digest, verify_source


def compact_plan(masks: list[tuple[int, ...]], rank: int) -> str:
    stats = scope_plan(masks, rank)
    if stats["exceeds_existing_limit"]:
        raise ValueError("Refusing a plan beyond the existing reference budget")
    tokens = [2, rank, stats["storage"], len(masks)]
    active = []
    offset = 0
    for terms in masks:
        mask = 0
        for term in terms:
            mask |= term
        variables = bits(mask)
        labels = []
        for local in range(1 << len(variables)):
            expanded = sum(((local >> i) & 1) << v for i, v in enumerate(variables))
            labels.append(
                sum(((expanded & term).bit_count() % 2) << j for j, term in enumerate(terms))
            )
        tokens += [offset, len(labels), *labels]
        active.append((mask, offset))
        offset += len(labels)
    tokens.append(rank)
    for pivot in stats["order"]:
        selected = [(scope, start) for scope, start in active if scope >> pivot & 1]
        union = 0
        for scope, _ in selected:
            union |= scope
        remaining = union ^ (1 << pivot)
        ordered = bits(remaining) + [pivot]
        size = 1 << len(ordered)
        tokens += [offset, size // 2, len(selected)]
        for scope, start in selected:
            positions = [ordered.index(v) for v in bits(scope)]
            low = [
                start + sum(((local >> p) & 1) << j for j, p in enumerate(positions))
                for local in range(min(size, 256))
            ]
            high = (
                []
                if size <= 256
                else [
                    sum(((block >> (p - 8)) & 1) << j for j, p in enumerate(positions) if p >= 8)
                    for block in range(size // 256)
                ]
            )
            tokens += [len(low), len(high), *low, *high]
        active = [(scope, start) for scope, start in active if not scope >> pivot & 1]
        active.append((remaining, offset))
        offset += size // 2
    tokens += [len(active), *(start for _, start in active)]
    if offset != stats["storage"]:
        raise AssertionError("Compact serialization storage differs from scope count")
    return "\n".join(map(str, tokens)) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify_source(args.source)
    args.work.mkdir(parents=True, exist_ok=False)
    ref = load_reference(args.reference)
    header_relative = "tools/profile/gadget_contraction_kernel.h"
    header = args.reference / header_relative
    if header.read_bytes() != subprocess.check_output(
        ["git", "show", REFERENCE + ":" + header_relative]
    ):
        raise ValueError("Unpinned numeric reference")
    harness = Path(__file__).with_name("profile_wan_contraction_scaling.cc")
    binary = args.work / "kernel"
    command = [
        "c++",
        "-std=c++20",
        "-O3",
        "-march=native",
        "-DNDEBUG",
        "-fno-fast-math",
        "-I" + str(header.parent),
        "-I" + str(Path(__file__).resolve().parents[2] / "src"),
        str(harness),
        "-o",
        str(binary),
    ]
    subprocess.run(command, check=True)
    variants = {}
    for d in (3, 7, 9):
        surface = ref.Surface(d)
        n = len(surface.points)
        rows = [sum(1 << q for q in support) for support in surface.checks("X")]
        rank = len(rows)
        masks = masks_for(rows, n, [pair for group in surface.pairs for pair in group])
        if d == 3:
            for k in range(rank + 1):
                paired = paired_masks(masks, rank, k)
                if compact_plan(paired, rank + k) != ref.kernel_plan(
                    ref.FactorPlan(paired, rank + k)
                ):
                    raise AssertionError("Compact compiler differs from expanded reference")
            continue
        if d == 7:
            variants["regular7_current"] = (masks, rank)
            centers = [(x, y) for y in range(0, 2 * d - 1, 2) for x in range(1, 2 * d - 1, 2)]
            order = sorted(range(rank), key=lambda i: (sum(centers[i]), centers[i]))
            variants["regular7_diagonal"] = (permute_masks(masks, order), rank)
        else:
            order = scope_plan(masks, rank)["order"]
            variants["regular9_amplitude"] = (permute_masks(masks, order), rank)
    masks, _ = sahay_masks(args.source)
    variants["sahay5_current"] = (masks, 20)
    variants["sahay5_amplitude"] = (permute_masks(masks, scope_plan(masks, 20)["order"]), 20)
    output: dict[str, Any] = dict(
        reference_revision=REFERENCE,
        harness_sha256=digest(harness),
        header_sha256=digest(header),
        binary_sha256=digest(binary),
        build_command=command,
        compact_validation_d3_plans=7,
        scope="one contraction per syndrome prefix; no physical binding or sampling",
        models={},
    )
    for name, (masks, rank) in variants.items():
        directory = args.work / name
        directory.mkdir()
        start = time.perf_counter()
        (directory / "metadata.txt").write_text(f"{rank}\n{len(masks)}\n")
        hashes = []
        stats = family(masks, rank)
        # Include a large compressed gather: d3 alone cannot exercise the
        # high-block address path because all its joint tables fit 256 entries.
        worst = max(range(rank + 1), key=lambda k: stats["marginals"][k]["peak_entries"])
        paired = paired_masks(masks, rank, worst)
        if compact_plan(paired, rank + worst) != ref.kernel_plan(
            ref.FactorPlan(paired, rank + worst)
        ):
            raise AssertionError("Compact compiler differs on a compressed gather")
        for k in range(rank + 1):
            path = directory / f"marginal_{k}.txt"
            path.write_text(compact_plan(paired_masks(masks, rank, k), rank + k))
            hashes.append(digest(path))
        prepare_seconds = time.perf_counter() - start

        def run(repeats: int) -> dict:
            result: dict = json.loads(
                subprocess.check_output(
                    ["taskset", "-c", "0", str(binary), str(directory), str(repeats)], text=True
                )
            )
            return result

        probe = run(2)
        repeats = max(4, min(512, int(0.3 * 2 / max(probe["seconds"], 1e-6))))
        trials = [run(repeats) for _ in range(5)]
        expected_lookup = sum(p["lookup_bytes"] for p in stats["marginals"])
        expected_workspace = 16 * (
            max(p["storage"] for p in stats["marginals"])
            + max(p["peak_entries"] for p in stats["marginals"])
        )
        for trial in trials:
            if (
                trial["lookup_bytes"] != expected_lookup
                or trial["workspace_bytes"] != expected_workspace
            ):
                raise AssertionError("Native allocation differs from offline prediction")
        result = dict(
            rank=rank,
            compact_validation_prefix=worst,
            compact_validation_peak_entries=stats["marginals"][worst]["peak_entries"],
            plan_sha256=hashes,
            prepare_seconds=prepare_seconds,
            trials=trials,
            median_ms_per_sweep=statistics.median(
                t["seconds"] * 1000 / t["sweeps"] for t in trials
            ),
            lookup_bytes=expected_lookup,
            workspace_bytes=expected_workspace,
        )
        output["models"][name] = result
        args.output.write_text(json.dumps(output, indent=2) + "\n")
        print(name, result["median_ms_per_sweep"], "ms", prepare_seconds, "s prepare", flush=True)


if __name__ == "__main__":
    main()
