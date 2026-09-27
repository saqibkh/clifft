"""Replay the pinned reconstructed full f7 protocol against a current-main build."""

from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
from pathlib import Path

from study_wan_scaling import REFERENCE
from wan_corpus import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--core-library", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    reference = args.reference / "tools/profile"
    args.work.mkdir(parents=True, exist_ok=False)
    files = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", REFERENCE, "tools/profile"], text=True
    ).splitlines()
    hashes = {}
    for name in files:
        path = args.reference / name
        if path.read_bytes() != subprocess.check_output(["git", "show", REFERENCE + ":" + name]):
            raise ValueError("Modified reference source: " + name)
        hashes[name] = digest(path)
    commands = []
    for source, binary in (
        ("sample_folded_protocol.cpp", "sampler"),
        ("replay_cultivation.cpp", "replay"),
    ):
        command = [
            "c++",
            "-std=c++20",
            "-O3",
            "-march=native",
            "-fopenmp",
            "-I" + str(root / "src"),
            "-I" + str(reference),
            str(reference / source),
            str(args.core_library),
            "-o",
            str(args.work / binary),
        ]
        commands.append(command)
        subprocess.run(command, check=True)
    # The old physical-history oracle imports Cirq. Its use is confined to
    # validation; the native sampler links only the current Clifft core.
    audit = [
        sys.executable,
        str(reference / "audit_folded_protocol.py"),
        "--sampler",
        str(args.work / "sampler"),
        "--reference",
        str(args.work / "replay"),
        "--distances",
        "7",
        "--probabilities",
        "0",
        "0.001",
        "0.03",
        "--cases",
        "4",
        "--oracle-cases",
        "4",
        "--schedule",
        "unbounded",
        "--native-growth",
        "--output",
        str(args.work / "audit.json"),
    ]
    commands.append(audit)
    subprocess.run(audit, check=True)
    bundle = args.work / "bundle"
    command = [
        sys.executable,
        str(reference / "compile_folded_protocol.py"),
        "--distance",
        "7",
        "--probability",
        "0.001",
        "--native-growth",
        "--output",
        str(bundle),
    ]
    commands.append(command)
    subprocess.run(command, check=True)
    cpu = min(os.sched_getaffinity(0))
    command = [
        "taskset",
        "-c",
        str(cpu),
        str(args.work / "sampler"),
        str(bundle),
        "32",
        "0",
        "18412",
        "0",
        "unbounded",
        "0",
    ]
    commands.append(command)
    benchmark = json.loads(subprocess.check_output(command, text=True))
    benchmark["median_seconds_per_attempt"] = statistics.median(benchmark["seconds"]) / 32
    metadata = json.loads((bundle / "metadata.json").read_text())
    result = dict(
        provenance="pinned local full f7 reconstruction, not an author-supplied d7 export",
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        reference_revision=REFERENCE,
        core_sha256=digest(args.core_library),
        sampler_sha256=digest(args.work / "sampler"),
        replay_sha256=digest(args.work / "replay"),
        driver_sha256=digest(Path(__file__)),
        reference_source_sha256=hashes,
        commands=commands,
        cpu_affinity=cpu,
        early_rejection=False,
        circuit_sha256=metadata["circuit_sha256"],
        visible_records=metadata["visible"],
        hidden_records=metadata["hidden"],
        audit=json.loads((args.work / "audit.json").read_text()),
        benchmark=benchmark,
        ordinary_coefficient_bytes=16 * (1 << benchmark["ordinary_peak_width"]),
        ordinary_sampling="skipped before allocating the dense coefficient array",
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Full reconstructed f7: {benchmark['median_seconds_per_attempt'] * 1000:.3f} ms")


if __name__ == "__main__":
    main()
