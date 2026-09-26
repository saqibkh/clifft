"""Interleaved CPU measurements of current Clifft on the pinned Wan corpus."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

from study_wan_corpus import summarize
from wan_corpus import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    sys.path.insert(0, str(args.package.resolve()))
    import clifft

    if Path(clifft.__file__).resolve().parent != args.package.resolve() / "clifft":
        raise ValueError("The requested Clifft package was shadowed")
    os.sched_setaffinity(0, {0})
    manifest = json.loads((args.corpus / "manifest.json").read_text())
    output = dict(
        main_revision=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], text=True
        ).strip(),
        study_revision=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        extension_sha256=digest(next((args.package / "clifft").glob("*.so"))),
        python=platform.python_version(),
        affinity=sorted(os.sched_getaffinity(0)),
        cpu=next(
            line.strip()
            for line in Path("/proc/cpuinfo").read_text().splitlines()
            if line.startswith("model name")
        ),
        p=manifest["p"],
        threads=1,
        batch_size=1,
        keep_records=True,
        cases={},
    )
    for case in ("fold-d5", "rp2-d5", "chan-four-round-d5"):
        path = args.corpus / (case + ".stim")
        text = path.read_text()
        meta = manifest["cases"][case]
        if digest(path) != meta["circuit_sha256"]:
            raise ValueError("Circuit hash differs from the pinned manifest")
        rows = {}
        programs = {}
        for mode in ("full", "postselected"):
            for label, budget in (("default", False), ("scheduled", 16), ("unlimited", None)):
                key = mode + "/" + label
                options = {}
                scheduler = None
                if budget is not False:
                    passes = clifft.default_hir_pass_manager()
                    scheduler = clifft.ActiveWidthSchedulePass(search_budget=budget)
                    passes.add(scheduler)
                    options["hir_passes"] = passes
                start = time.perf_counter()
                program = clifft.compile(
                    text,
                    expected_detectors=meta["expected_detectors"],
                    expected_observables=meta["expected_observables"],
                    postselection_mask=[1] * len(meta["expected_detectors"])
                    if mode == "postselected"
                    else None,
                    **options,
                )
                compile_seconds = time.perf_counter() - start
                if program.peak_active_width > 24:
                    raise ValueError("Dense worker exceeds the study's memory bound")
                shots = (
                    (32 if label == "default" else 128)
                    if case == "fold-d5"
                    else (1024 if case == "rp2-d5" else 32768)
                )
                row = dict(
                    compile_seconds=compile_seconds,
                    peak_active_width=program.peak_active_width,
                    executable_sha256=hashlib.sha256(program.inspect().encode()).hexdigest(),
                    shots_per_trial=shots,
                    trials=[],
                )
                if scheduler is not None:
                    row["scheduler"] = {
                        name: getattr(scheduler, name)
                        for name in (
                            "applied",
                            "incumbent_peak",
                            "result_peak",
                            "incumbent_dense_work",
                            "result_dense_work",
                            "swept_ops",
                            "classification_probes",
                        )
                    }
                sampler = clifft.sample_survivors if mode == "postselected" else clifft.sample
                kwargs = dict(threads=1, batch_size=1)
                if mode == "postselected":
                    kwargs["keep_records"] = True
                sampler(program, 2, seed=910, **kwargs)
                programs[key] = program, sampler, kwargs
                rows[key] = row
        output["cases"][case] = dict(circuit_sha256=digest(path), rows=rows)
        keys = list(rows)
        for trial in range(args.repeats):
            # Rotate and reverse order so a configuration does not always run first.
            order = keys[trial % len(keys) :] + keys[: trial % len(keys)]
            if trial % 2:
                order.reverse()
            for key in order:
                row = rows[key]
                program, sampler, kwargs = programs[key]
                start = time.perf_counter()
                result = sampler(program, row["shots_per_trial"], seed=91000 + trial, **kwargs)
                seconds = time.perf_counter() - start
                counts = summarize(result, survivors=key.startswith("postselected"))
                row["trials"].append(dict(seconds=seconds, seed=91000 + trial, **counts))
                row["median_us_per_attempt"] = statistics.median(
                    t["seconds"] / t["attempts"] * 1e6 for t in row["trials"]
                )
                args.output.write_text(json.dumps(output, indent=2) + "\n")
                print(case, key, trial, round(seconds, 3), counts["accepted"], flush=True)


if __name__ == "__main__":
    main()
