"""Attribute the scheduled Sahay circuit's execution cost to source regions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Any

from wan_corpus import Exporter, digest, verify_source


def source_regions(source: Path, circuit: Path) -> dict[int, str]:
    model = json.loads((source / "data/alternatives/inputs/fold-d5/model.json").read_text())
    events = model["events"]
    first = next(i for i, e in enumerate(events) if e.get("native_block") == 2)
    second = next(i for i, e in enumerate(events) if e.get("native_block") == 3)
    terminal = next(i for i, e in enumerate(events) if e["name"] == "MEASURE_PRODUCT")
    exporter = Exporter(model, 0.001)
    mapping = {}
    for index, event in enumerate(events):
        begin = len(exporter.lines)
        exporter.dictionary_event(event, index)
        label = (
            "prefix"
            if index < first
            else "first_large_check"
            if index < second
            else "second_large_check"
            if index < terminal
            else "terminal_projection"
        )
        for line in range(begin + 1, len(exporter.lines) + 1):
            mapping[line] = label
    exported, _ = Exporter(model, 0.001).export()
    if exported != circuit.read_text():
        raise ValueError("Source regions do not describe the timed circuit")
    for line in range(len(exporter.lines) + 1, len(exported.splitlines()) + 1):
        mapping[line] = "logical_readout"
    return mapping


def normalized_plan(path: Path) -> str:
    return re.sub(r" plans=\[\d+,\d+\)", "", path.read_text())


def rotation_runs(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Describe fixed-width runs; rank is an exact binary mask calculation."""
    runs = []
    i = 0
    while i < len(rows):
        first = rows[i]
        if not first["instruction"].startswith("ROTATE "):
            i += 1
            continue
        j = i
        basis: dict[int, int] = {}
        total_ns = 0
        has_z = False
        while (
            j < len(rows)
            and rows[j]["instruction"].startswith("ROTATE ")
            and rows[j]["active_width"] == first["active_width"]
        ):
            row = rows[j]
            instruction = row["instruction"]
            x = sum(1 << int(q) for q in re.findall(r"[XY](\d+)", instruction))
            has_z |= bool(re.search(r"[ZY]\d+", instruction))
            total_ns += row["accepted_ns"] + row["rejected_ns"]
            while x:
                pivot = x.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = x
                    break
                x ^= basis[pivot]
            j += 1
        runs.append(
            dict(
                begin=i,
                end=j,
                length=j - i,
                active_width=first["active_width"],
                x_rank=len(basis),
                all_x=not has_z,
                nanoseconds=total_ns,
            )
        )
        i = j
    return sorted(runs, key=lambda r: -r["nanoseconds"])


def collect(prefix: Path, mapping: dict[int, str]) -> dict:
    summary: dict[str, Any] = json.loads(prefix.with_suffix(".json").read_text())
    rows = []
    totals: dict[str, dict[str, dict[str, int]]] = {}
    with prefix.with_suffix(".actions.tsv").open() as stream:
        for raw in csv.DictReader(stream, delimiter="\t"):
            sources = list(map(int, raw["source_lines"].split(","))) if raw["source_lines"] else []
            regions = sorted({mapping[line] for line in sources})
            region = "+".join(regions) if regions else "no_source"
            width = raw["active_width"]
            kind = raw["instruction"].split()[0]
            kernel = re.search(r"kernel=(\w+)", raw["instruction"])
            row: dict[str, Any] = {
                k: int(raw[k])
                for k in (
                    "action",
                    "accepted_ns",
                    "rejected_ns",
                    "accepted_visits",
                    "rejected_visits",
                    "active_width",
                )
            }
            row.update(region=region, source_lines=sources, instruction=raw["instruction"])
            rows.append(row)
            for dimension, value in (
                ("region", region),
                ("kind", kind),
                ("width", width),
                ("kernel", kernel[1] if kernel else "none"),
            ):
                target = totals.setdefault(dimension, {}).setdefault(
                    value, dict(accepted_ns=0, rejected_ns=0, accepted_visits=0, rejected_visits=0)
                )
                for field in target:
                    target[field] += int(row[field])
    summary.update(groups=totals, actions=rows, rotation_runs=rotation_runs(rows))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--perf", type=Path)
    args = parser.parse_args()
    verify_source(args.source)
    args.raw.mkdir(parents=True, exist_ok=False)
    circuit = args.corpus / "fold-d5.stim"
    manifest = json.loads((args.corpus / "manifest.json").read_text())
    meta = manifest["cases"]["fold-d5"]
    if digest(circuit) != meta["circuit_sha256"] or manifest["p"] != 0.001:
        raise ValueError("Unexpected circuit or noise probability")
    references = args.raw / "references.txt"
    references.write_text(
        "".join(map(str, meta["expected_detectors"]))
        + "\n"
        + "".join(map(str, meta["expected_observables"]))
        + "\n"
    )
    mapping = source_regions(args.source, circuit)
    output = dict(
        circuit_sha256=digest(circuit),
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        build=json.loads((args.native / "build.json").read_text()),
        source_regions={
            name: [
                min(k for k, v in mapping.items() if v == name),
                max(k for k, v in mapping.items() if v == name),
            ]
            for name in sorted(set(mapping.values()))
        },
        profiles={},
    )
    import sys

    sys.path.insert(0, str(args.package.resolve()))
    import clifft

    if Path(clifft.__file__).resolve().parent != args.package.resolve() / "clifft":
        raise ValueError("The requested package was shadowed")
    for mode, shots in (("full", 256), ("postselected", 1024)):
        postselect = int(mode == "postselected")
        prefixes = {}
        for variant in ("untimed", "timed"):
            prefix = args.raw / (mode + "-" + variant)
            prefixes[variant] = prefix
            command = [
                "taskset",
                "-c",
                "0",
                str(args.native / variant),
                str(circuit),
                "unlimited",
                str(postselect),
                str(shots),
                "91301",
                str(prefix),
                str(references),
            ]
            if args.perf and variant == "untimed":
                command = [
                    str(args.perf),
                    "record",
                    "-e",
                    "cycles:u",
                    "-F",
                    "499",
                    "-o",
                    str(args.raw / (mode + ".perf")),
                    "--",
                    *command,
                ]
            subprocess.run(command, check=True)
            print(mode, variant, prefix.with_suffix(".json").read_text().strip(), flush=True)
        untimed = json.loads(prefixes["untimed"].with_suffix(".json").read_text())
        timed = collect(prefixes["timed"], mapping)
        for field in ("shots", "accepted", "errors", "record_hash", "peak_active_width"):
            if timed[field] != untimed[field]:
                raise AssertionError(("Instrumentation changed results", field))
        if normalized_plan(prefixes["timed"].with_suffix(".plan.txt")) != normalized_plan(
            prefixes["untimed"].with_suffix(".plan.txt")
        ):
            raise AssertionError("Instrumented and ordinary native plans differ")
        passes = clifft.default_hir_pass_manager()
        passes.add(clifft.ActiveWidthSchedulePass(search_budget=None))
        program = clifft.compile(
            circuit.read_text(),
            hir_passes=passes,
            expected_detectors=meta["expected_detectors"],
            expected_observables=meta["expected_observables"],
            postselection_mask=[1] * len(meta["expected_detectors"]) if postselect else None,
        )
        native_plan = normalized_plan(prefixes["timed"].with_suffix(".plan.txt"))
        if native_plan != program.inspect():
            raise AssertionError("Native profiler and Python public API compile different plans")
        actions_path = args.output.with_name(args.output.stem + "_" + mode + "_actions.json")
        actions_path.write_text(
            "[\n" + ",\n".join(json.dumps(row) for row in timed.pop("actions")) + "\n]\n"
        )
        output["profiles"][mode] = dict(
            untimed=untimed,
            timed=timed,
            actions_file=actions_path.name,
            actions_sha256=digest(actions_path),
            public_api_plan_sha256=hashlib.sha256(native_plan.encode()).hexdigest(),
            instrumentation_wall_ratio=timed["nanoseconds"] / untimed["nanoseconds"],
        )
        if args.perf:
            report = subprocess.check_output(
                [
                    str(args.perf),
                    "report",
                    "--stdio",
                    "--no-children",
                    "--percent-limit",
                    "0.5",
                    "-i",
                    str(args.raw / (mode + ".perf")),
                ],
                text=True,
            )
            (args.raw / (mode + ".perf.txt")).write_text(report)
            output["profiles"][mode]["perf"] = dict(
                event="cycles:u",
                frequency_hz=499,
                scope="Whole process including compilation and startup",
                summary_lines=[
                    line.strip() for line in report.splitlines() if "%  untimed" in line
                ],
            )
        args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
