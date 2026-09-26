"""Build isolated timed/untimed research executors from a current CMake build.

Production sources and the production library are never modified. The single
replacement object adds preallocated, per-action wall-clock counters. All other
objects and compiler flags come from the existing optimized Clifft build.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from wan_corpus import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    directory = Path(__file__).resolve().parent
    compilation = json.loads((args.build / "compile_commands.json").read_text())
    entry = next(e for e in compilation if e["file"].endswith("/sampling/executor.cc"))
    original = Path(entry["file"])
    text = original.read_text()
    before = "            const ExecutablePlan::Action& action = plan_->actions_[action_index];"
    after = "                action);\n            if constexpr (Mode == ShotMode::ReplayRecords)"
    if text.count(before) != 1 or text.count(after) != 1:
        raise ValueError("Executor layout changed; audit instrumentation before rebuilding")
    text = '#include "wan_action_timer.h"\n' + text.replace(
        before, "            const auto wan_start = wan_profile::Clock::now();\n" + before
    ).replace(
        after,
        "                action);\n"
        "            wan_profile::record(action_index, wan_start);\n"
        "            if constexpr (Mode == ShotMode::ReplayRecords)",
    )
    instrumented = args.output / "executor_timed.cc"
    instrumented.write_text(text)
    command = shlex.split(entry["command"])
    # Match the optimized executor's flags, removing only its source/output pair.
    output_index = command.index("-o")
    del command[output_index : output_index + 2]
    command.remove("-c")
    command.remove(str(original))
    command += ["-I" + str(directory)]
    executor_object = args.output / "executor_timed.o"
    subprocess.run(command + ["-c", str(instrumented), "-o", str(executor_object)], check=True)
    harness_object = args.output / "profile_wan_actions.o"
    subprocess.run(
        command + ["-c", str(directory / "profile_wan_actions.cc"), "-o", str(harness_object)],
        check=True,
    )
    library = args.build / "src/clifft/libclifft_core.a"
    for label, objects in (("timed", [executor_object]), ("untimed", [])):
        subprocess.run(
            [
                command[0],
                "-fopenmp",
                str(harness_object),
                *map(str, objects),
                str(library),
                "-o",
                str(args.output / label),
            ],
            check=True,
        )
    (args.output / "build.json").write_text(
        json.dumps(
            dict(
                original_executor_sha256=digest(original),
                instrumented_executor_sha256=digest(instrumented),
                core_library_sha256=digest(library),
                compile_command=entry["command"],
                harness_sha256=digest(directory / "profile_wan_actions.cc"),
                timer_sha256=digest(directory / "wan_action_timer.h"),
                binaries={label: digest(args.output / label) for label in ("timed", "untimed")},
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
