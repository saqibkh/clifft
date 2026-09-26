"""Allocation-bounded planning study for larger folded-code distances.

Count contraction scopes and address tables without materializing tensors.
Reference circuits and the experimental planner are pinned external inputs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from wan_corpus import digest, verify_source

REFERENCE = "b87c0da396b4bbde54fbee5535aeac1116170605"
CONTRACTION = "src/python/clifft/_folded/contraction.py"


def bits(value: int) -> list[int]:
    result = []
    while value:
        low = value & -value
        result.append(low.bit_length() - 1)
        value ^= low
    return result


def scope_plan(masks: list[tuple[int, ...]], rank: int) -> dict[str, Any]:
    factors = [0] * len(masks)
    for i, terms in enumerate(masks):
        for term in terms:
            factors[i] |= term
    neighbors = [0] * rank
    leaf_entries = sum(1 << scope.bit_count() for scope in factors)
    for scope in factors:
        for v in bits(scope):
            neighbors[v] |= scope ^ (1 << v)
    remaining = (1 << rank) - 1
    storage = leaf_entries
    work = gathers = compressed_gathers = peak = 0
    order = []
    while remaining:
        pivot = min(bits(remaining), key=lambda v: (neighbors[v].bit_count(), v))
        bit = 1 << pivot
        union = neighbors[pivot] | bit
        selected = [scope for scope in factors if scope & bit]
        if not selected:
            raise ValueError("Elimination variable has no factors")
        size = 1 << union.bit_count()
        peak = max(peak, union.bit_count())
        work += size
        gathers += size * len(selected)
        # Gather addresses are sums of disjoint bit positions, so the reference
        # planner's 256-entry decomposition is exact without enumerating them.
        compressed_gathers += len(selected) * (size if size <= 256 else 256 + size // 256)
        storage += size // 2
        factors = [scope for scope in factors if not scope & bit]
        factors.append(union ^ bit)
        for v in bits(neighbors[pivot]):
            neighbors[v] = (neighbors[v] | union) & ~(bit | (1 << v))
        remaining ^= bit
        order.append(pivot)
    if any(factors):
        raise AssertionError("Incomplete elimination")
    return dict(
        rank=rank,
        peak_scope=peak,
        peak_entries=1 << peak,
        storage=storage,
        work=work,
        gather_entries=gathers,
        compact_gather_entries=compressed_gathers,
        lookup_bytes=leaf_entries + 4 * (compressed_gathers + len(factors)),
        leaf_entries=leaf_entries,
        order=order,
        exceeds_existing_limit=storage > 2_000_000 or gathers > 2_000_000,
    )


def masks_for(x_rows: list[int], n: int, pairs: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    physical = [sum(((row >> q) & 1) << i for i, row in enumerate(x_rows)) for q in range(n)]
    result: list[tuple[int, ...]] = [(m,) for m in physical]
    result.extend((physical[a], physical[b]) for a, b in pairs)
    result.extend((1 << i,) for i in range(len(x_rows)))
    return result


def permute_masks(masks: list[tuple[int, ...]], order: list[int]) -> list[tuple[int, ...]]:
    inverse = {v: i for i, v in enumerate(order)}
    return [tuple(sum(1 << inverse[v] for v in bits(m)) for m in terms) for terms in masks]


def paired_masks(masks: list[tuple[int, ...]], rank: int, measured: int) -> list[tuple[int, ...]]:
    selected = (1 << measured) - 1
    return masks + [
        tuple((m & ~selected) | ((m & selected) << rank) for m in terms) for terms in masks
    ]


def family(masks: list[tuple[int, ...]], rank: int) -> dict[str, Any]:
    start = time.perf_counter()
    amplitude = scope_plan(masks, rank)
    marginals = [scope_plan(paired_masks(masks, rank, k), rank + k) for k in range(rank + 1)]
    plans = [amplitude, *marginals]
    scratch = max(p["storage"] for p in plans)
    product = max(p["peak_entries"] for p in plans)
    return dict(
        amplitude=amplitude,
        marginals=marginals,
        peak_entries=product,
        sum_marginal_work=sum(p["work"] for p in marginals),
        numeric_payload_bytes=sum(p["lookup_bytes"] for p in plans) + 16 * (scratch + product),
        workspace_bytes=16 * (scratch + product),
        lookup_bytes=sum(p["lookup_bytes"] for p in plans),
        plans_exceeding_existing_limit=sum(p["exceeds_existing_limit"] for p in plans),
        planning_seconds=time.perf_counter() - start,
    )


def evaluate_reference(plan: Any, local: Any) -> complex:
    import numpy as np

    scratch = np.empty(plan.storage, dtype=np.complex128)
    for (offset, labels), values in zip(plan.leaves, local, strict=True):
        scratch[offset : offset + len(labels)] = values[labels]
    for offset, size, gathers in plan.steps:
        product = np.ones(2 * size, dtype=np.complex128)
        for gather in gathers:
            product *= scratch[gather]
        scratch[offset : offset + size] = product[:size] + product[size:]
    return complex(np.prod(scratch[plan.outputs]) / (1 << plan.rank))


def validate_reordering(ref: Any, masks: list[tuple[int, ...]], rank: int) -> dict:
    """Compare every prefix marginal with independently enumerated amplitudes."""
    import numpy as np

    rng = np.random.default_rng(49509)
    orders = [list(range(rank)), list(reversed(range(rank))), scope_plan(masks, rank)["order"]]
    max_error = 0.0
    comparisons = 0
    for _ in range(4):
        local = np.exp(1j * np.pi / 4 * rng.integers(0, 8, size=(len(masks), 4)))
        local[-rank:] = 1
        values = np.array(
            [
                np.prod(
                    [
                        table[
                            sum(((g & mask).bit_count() % 2) << k for k, mask in enumerate(terms))
                        ]
                        for table, terms in zip(local, masks, strict=True)
                    ]
                )
                for g in range(1 << rank)
            ]
        )
        amplitudes = np.array(
            [
                sum(values[g] * (-1) ** ((s & g).bit_count() % 2) for g in range(1 << rank))
                / (1 << rank)
                for s in range(1 << rank)
            ]
        )
        probabilities = np.abs(amplitudes) ** 2
        np.testing.assert_allclose(probabilities.sum(), 1, atol=1e-12)
        for order in orders:
            permuted = permute_masks(masks, order)
            for k in range(rank + 1):
                plan = ref.FactorPlan(paired_masks(permuted, rank, k), rank + k)
                for syndrome in range(1 << k):
                    signed = local.copy()
                    for j, original in enumerate(order[:k]):
                        signed[-rank + original, 1] = (-1) ** ((syndrome >> j) & 1)
                    actual = evaluate_reference(plan, np.concatenate([signed, signed.conjugate()]))
                    expected = sum(
                        probabilities[s]
                        for s in range(1 << rank)
                        if all(
                            ((s >> original) & 1) == ((syndrome >> j) & 1)
                            for j, original in enumerate(order[:k])
                        )
                    )
                    max_error = max(max_error, abs(actual - expected))
                    comparisons += 1
    if max_error > 1e-12:
        raise AssertionError(("Reordered Born marginal mismatch", max_error))
    return dict(
        cases=4,
        orders=len(orders),
        marginal_comparisons=comparisons,
        maximum_absolute_error=max_error,
        oracle="exhaustive Fourier amplitudes",
    )


def load_reference(root: Path) -> Any:
    path = root / CONTRACTION
    expected = subprocess.check_output(["git", "show", REFERENCE + ":" + CONTRACTION])
    if path.read_bytes() != expected:
        raise ValueError("Reference planner differs from pinned revision")
    spec = importlib.util.spec_from_file_location("wan_scaling_reference", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def validate_reference(ref: Any, masks: list[tuple[int, ...]], rank: int, result: dict) -> dict:
    selected = [(masks, rank, result["amplitude"])]
    ks = (
        list(range(rank + 1))
        if rank <= 20
        else [max(range(rank + 1), key=lambda k: result["marginals"][k]["peak_entries"])]
    )
    selected += [(paired_masks(masks, rank, k), rank + k, result["marginals"][k]) for k in ks]
    for local, width, expected in selected:
        actual = ref.FactorPlan(local, width)
        for field in ("rank", "storage", "peak_scope", "work", "gather_entries"):
            if getattr(actual, field) != expected[field]:
                raise AssertionError((field, getattr(actual, field), expected[field]))
        # Count the exact compact serialization, including narrow leaf labels.
        lookup = sum(len(parity) for _, parity in actual.leaves) + 4 * len(actual.outputs)
        for _, _, gathers in actual.steps:
            for gather in gathers:
                if len(gather) > 256:
                    blocks = gather.reshape(-1, 256)
                    import numpy as np

                    if not np.array_equal(
                        blocks, gather[:256] + (blocks[:, 0] - gather[0])[:, None]
                    ):
                        raise AssertionError("Gather is not separable")
                    lookup += 4 * (256 + len(gather) // 256)
                else:
                    lookup += 4 * len(gather)
        if lookup != expected["lookup_bytes"]:
            raise AssertionError(("lookup_bytes", lookup, expected["lookup_bytes"]))
    return dict(plans_compared=len(selected), exact_dimensions_and_compact_addresses=True)


def sahay_masks(source: Path) -> tuple[list[tuple[int, ...]], dict]:
    path = source / "data/alternatives/inputs/fold-d5/model.json"
    model = json.loads(path.read_text())
    checks = [e for e in model["events"] if e["name"] == "MEASURE_PRODUCT"]
    data = sorted({q for e in checks for q in e["support"]})
    indices = {q: i for i, q in enumerate(data)}
    rows = {
        axis: [sum(1 << indices[q] for q in e["support"]) for e in checks if e["axis"] == axis]
        for axis in "XZ"
    }
    if any((a & b).bit_count() % 2 for a in rows["X"] for b in rows["Z"]):
        raise AssertionError("Data checks do not commute")
    pairs = sorted(
        {
            tuple(sorted(indices[q] for q in e["support"][1:]))
            for e in model["events"]
            if e.get("native_block") == 2 and e["name"] == "CCZ"
        }
    )
    other = sorted(
        {
            tuple(sorted(indices[q] for q in e["support"][1:]))
            for e in model["events"]
            if e.get("native_block") == 3 and e["name"] == "CCZ"
        }
    )
    if pairs != other:
        raise ValueError("Two checks have different interaction supports")
    return masks_for(rows["X"], len(data), pairs), dict(
        model_sha256=digest(path),
        data_qubits=len(data),
        x_generators=len(rows["X"]),
        interaction_pairs=pairs,
        measurement_order="supplied terminal X checks",
        scope="Factor topology only; complete physical-instrument binding is unimplemented",
    )


def compile_controls(clifft: Any, reference: Path) -> dict:
    rows = {}
    for d in (3, 5, 7):
        path = reference / f"tests/fixtures/folded/f{d}.stim"
        if path.read_bytes() != subprocess.check_output(
            ["git", "show", REFERENCE + ":" + str(path.relative_to(reference))]
        ):
            raise ValueError("Unpinned reconstruction")
        row: dict[str, Any] = dict(
            circuit_sha256=digest(path), scope="paper-guided full reconstruction", plans={}
        )
        for label, budget in (("default", False), ("scheduled", 16), ("unlimited", None)):
            start = time.perf_counter()
            hir = clifft.trace(clifft.parse(path.read_text()))
            passes = clifft.default_hir_pass_manager()
            passes.run(hir)
            if budget is not False:
                scheduler = clifft.ActiveWidthSchedulePass(search_budget=budget)
                scheduling = clifft.HirPassManager()
                scheduling.add(scheduler)
                scheduling.run(hir)
            # lower() avoids the public compile API's default reference shot.
            # No Executor or exponentially sized coefficient state is constructed.
            program = clifft.lower(hir)
            row["plans"][label] = dict(
                peak_active_width=program.peak_active_width,
                coefficient_bytes=16 << program.peak_active_width,
                compile_seconds=time.perf_counter() - start,
                actions=program.num_actions,
            )
        rows[str(d)] = row
        print("compiled", d, row["plans"], flush=True)
    return rows


def save_output(output: dict[str, Any], path: Path) -> None:
    """Keep reproducible per-distance evidence within repository artifact limits."""

    def compact(value: Any) -> Any:
        if isinstance(value, dict):
            return {
                ("order_sha256" if key == "order" else key): (
                    hashlib.sha256(json.dumps(item).encode()).hexdigest()
                    if key == "order"
                    else compact(item)
                )
                for key, item in value.items()
            }
        if isinstance(value, list):
            return [compact(item) for item in value]
        return value

    evidence = compact(output)
    for d, item in evidence["regular"].items():
        detail = path.with_name(f"{path.stem}_d{d}.json")
        detail.write_text(json.dumps(item, indent=2) + "\n")
        evidence["regular"][d] = dict(file=detail.name, sha256=digest(detail))
    path.write_text(json.dumps(evidence, indent=2) + "\n")


def fusion_screen() -> dict[str, Any]:
    from study_wan_bottleneck import rotation_runs

    root = Path(__file__).with_name("research")
    actions_path = root / "wan_bottleneck_full_actions.json"
    profile_path = root / "wan_bottleneck.json"
    actions = json.loads(actions_path.read_text())
    profile = json.loads(profile_path.read_text())["profiles"]["full"]["timed"]
    total = profile["nanoseconds"]
    runs = []
    for run in rotation_runs(actions)[:4]:
        if not run["all_x"]:
            raise ValueError("Rotation cost model requires commuting X masks")
        rank, length = run["x_rank"], run["length"]
        runs.append(
            dict(
                **run,
                fraction_of_full_time=run["nanoseconds"] / total,
                dense_matrix_entries=1 << (2 * rank),
                diagonal_phase_entries=1 << rank,
                separate_real_ops_per_coefficient=6 * length,
                walsh_real_ops_per_coefficient=4 * rank + 6,
                coefficient_bytes=16 << run["active_width"],
            )
        )
    return dict(
        actions_sha256=digest(actions_path),
        profile_sha256=digest(profile_path),
        scope="Arithmetic cost hypothesis; no fused kernel or measured speedup",
        assumptions="6 real operations per complex rotation output; two Walsh transforms "
        "plus complex phase; excludes addressing, normalization setup and phase binding",
        top_runs=runs,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify_source(args.source)
    ref = load_reference(args.reference)
    sys.path.insert(0, str(args.package.resolve()))
    import clifft

    if Path(clifft.__file__).resolve().parent != args.package.resolve() / "clifft":
        raise ValueError("A different Clifft package shadowed the selected build")
    output: dict[str, Any] = dict(
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        reference_revision=REFERENCE,
        reference_planner_sha256=digest(args.reference / CONTRACTION),
        extension_sha256=digest(next((args.package / "clifft").glob("*.so"))),
        compile_only=compile_controls(clifft, args.reference),
        fusion_screen=fusion_screen(),
        regular={},
    )
    for d in (3, 5, 7, 9, 11):
        surface = ref.Surface(d)
        n = len(surface.points)
        rows = [sum(1 << q for q in support) for support in surface.checks("X")]
        rank = len(rows)
        pairs = [pair for group in surface.pairs for pair in group]
        masks = masks_for(rows, n, pairs)
        result = family(masks, rank)
        item: dict[str, Any] = dict(
            data_qubits=n, x_generators=rank, default=result, alternatives={}
        )
        if d <= 7:
            item["validation"] = validate_reference(ref, masks, rank, result)
        if d == 3:
            output["reordering_validation"] = validate_reordering(ref, masks, rank)
        centers = [(x, y) for y in range(0, 2 * d - 1, 2) for x in range(1, 2 * d - 1, 2)]
        orders = dict(
            column=sorted(range(rank), key=lambda i: centers[i]),
            diagonal=sorted(range(rank), key=lambda i: (sum(centers[i]), centers[i])),
            amplitude_elimination=result["amplitude"]["order"],
        )
        if d >= 5:
            for name, order in orders.items():
                item["alternatives"][name] = family(permute_masks(masks, order), rank)
        output["regular"][str(d)] = item
        save_output(output, args.output)
        print(
            "scopes",
            d,
            result["peak_entries"],
            result["numeric_payload_bytes"],
            result["sum_marginal_work"],
            flush=True,
        )
    masks, provenance = sahay_masks(args.source)
    result = family(masks, provenance["x_generators"])
    output["sahay_d5"] = dict(
        provenance=provenance, plans=result, validation=validate_reference(ref, masks, 20, result)
    )
    output["sahay_d5"]["amplitude_order"] = family(
        permute_masks(masks, result["amplitude"]["order"]), 20
    )
    save_output(output, args.output)
    print("Sahay d5", result["peak_entries"], result["numeric_payload_bytes"], flush=True)


if __name__ == "__main__":
    main()
