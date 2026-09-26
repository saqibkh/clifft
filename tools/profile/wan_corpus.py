"""Export the pinned Wan-Zapirain frozen models without changing noise sites.

Inputs remain in an external checkout. Generated circuits and their provenance
manifest are research artifacts, not an additional production input format.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

REVISION = "b99452de7e300324e7bbf5cd9347241c3e289bae"
REPOSITORY = "https://github.com/kh428/exact_ler_for_msc"
CASES = (
    "rp2-d3",
    "rp2-d5",
    "fold-d3",
    "fold-d5",
    "chan-v1-d3",
    "chan-v1-d5",
    "chan-four-round-d5",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_source(root: Path) -> None:
    revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True)
    if revision.strip() != REVISION:
        raise ValueError("Use the pinned upstream revision " + REVISION)
    committed_manifest = subprocess.check_output(
        ["git", "-C", str(root), "show", REVISION + ":manifest.json"]
    )
    if (root / "manifest.json").read_bytes() != committed_manifest:
        raise ValueError("The upstream manifest has local edits")
    # Also detect local edits: the commit alone does not identify the input bytes.
    for item in json.loads((root / "manifest.json").read_text())["files"]:
        path = (root / item["path"]).resolve()
        if not path.is_relative_to(root.resolve()) or digest(path) != item["sha256"]:
            raise ValueError("Upstream manifest mismatch: " + item["path"])


def gate_lines(name: str, support: list[int]) -> list[str]:
    if name in {"CSX", "CS_DAG_X"}:
        control, target = support
        first, last = ("T_DAG", "T") if name == "CSX" else ("T", "T_DAG")
        # The conditional phase on the control is essential on a GHZ state.
        return [
            f"{last} {control}",
            f"{first} {target}",
            f"CX {control} {target}",
            f"{last} {target}",
        ]
    if name == "INJECT_HXY":
        return [f"{gate} {support[0]}" for gate in ("S", "T", "H", "T")]
    if name not in {"R", "RX", "H", "CX", "CCZ", "T", "T_DAG", "S", "S_DAG", "X", "Y", "Z"}:
        raise ValueError("Unsupported model gate " + name)
    return [name + " " + " ".join(map(str, support))]


def pauli_text(axes: dict[int, str]) -> str:
    return "*".join(f"{axis}{q}" for q, axis in sorted(axes.items()))


class Exporter:
    def __init__(self, model: dict, p: float, faults: dict[int, int] | None = None):
        self.model, self.p, self.faults = model, p, faults
        self.lines: list[str] = []
        self.records: dict[int, int] = {}
        self.measurements = 0
        self.expected_detectors: list[int] = []
        self.expected_observables = [0]
        self.noise_locations: list[dict] = []

    def rec(self, record: int) -> str:
        return f"rec[{self.records[record] - self.measurements}]"

    def measure(self, text: str, record: int | None = None) -> None:
        if record is not None:
            if record in self.records:
                raise ValueError("Repeated source record")
            self.records[record] = self.measurements
        self.lines.append(text)
        self.measurements += 1

    def detector(self, records: list[int], parity: int = 0) -> None:
        self.lines.append("DETECTOR " + " ".join(map(self.rec, records)))
        self.expected_detectors.append(parity)

    def noise(self, location: int, event_index: int) -> None:
        channel = self.model["channels"][location]
        support = channel["support"]
        if channel.get("source_kind") == "READOUT_FLIP":
            if self.records[channel["record"]] != self.measurements - 1:
                raise ValueError("Readout flip is not adjacent to its measurement")
            gate, target = self.lines[-1].split()
            if gate not in {"MZ", "MX"} or target != str(support[0]):
                raise ValueError("Unexpected readout flip target")
            self.noise_locations.append(
                dict(
                    location=location,
                    event_index=event_index,
                    line=len(self.lines),
                    gate="READOUT_FLIP",
                )
            )
            if self.faults is None:
                self.lines[-1] = f"{gate}({self.p:.17g}) {target}"
            elif location in self.faults:
                self.lines[-1] = f"{gate} !{target}"
            return
        axis = channel.get("axis") or channel.get("source_kind", "X_ERROR").split("_")[0]
        name = axis + "_ERROR" if channel["kind"] == "flip" else f"DEPOLARIZE{len(support)}"
        self.noise_locations.append(
            dict(location=location, event_index=event_index, line=len(self.lines) + 1, gate=name)
        )
        if self.faults is None:
            self.lines.append(f"{name}({self.p:.17g}) " + " ".join(map(str, support)))
        elif location in self.faults:
            code = self.faults[location]
            for j, q in enumerate(support):
                pauli = axis if channel["kind"] == "flip" else "IXZY"[(code >> (2 * j)) & 3]
                if pauli != "I":
                    self.lines.append(f"{pauli} {q}")

    def dictionary_event(self, event: dict, index: int) -> None:
        name = event["name"]
        support = event.get("support", [])
        if name == "NOISE":
            self.noise(event["location"], index)
        elif name == "MEASURE":
            target = ("!" if event["inverted"] else "") + str(support[0])
            self.measure("M" + event["axis"] + " " + target, event["record"])
        elif name == "MEASURE_PRODUCT":
            self.measure(
                "MPP " + pauli_text(dict.fromkeys(support, event["axis"])), event["record"]
            )
        elif name == "POST_Z":
            self.measure(f"M {support[0]}")
            self.lines.append("DETECTOR rec[-1]")
            self.expected_detectors.append(0)
        elif name == "DETECTOR":
            self.detector(event["records"], event["reference_parity"])
        elif name.startswith("FEEDBACK_"):
            self.lines.append(f"C{name[-1]} {self.rec(event['record'])} {support[0]}")
        else:
            self.lines.extend(gate_lines(name, support))

    def list_event(self, event: list, index: int) -> None:
        name, *args = event
        if name == "NOISE":
            self.noise(args[0], index)
        elif name == "RESET":
            q, axis, sign = args
            if sign != 1 or axis not in {"Z", "X"}:
                raise ValueError("Unsupported signed reset")
            self.lines.append(f"R{axis} {q}" if axis == "X" else f"R {q}")
        elif name == "MEASURE":
            q, axis, record = args
            self.measure(f"M{axis} {q}", record)
        elif name == "MEASURE_PAULI":
            x, z, record = args
            axes = {
                q: "IXZY"[((x >> q) & 1) + 2 * ((z >> q) & 1)]
                for q in range((x | z).bit_length())
                if (x | z) >> q & 1
            }
            self.measure("MPP " + pauli_text(axes), record)
        elif name == "FEEDBACK":
            axis, q, record = args
            self.lines.append(f"C{axis} {self.rec(record)} {q}")
        else:
            self.lines.extend(gate_lines(name, args))

    def export(self) -> tuple[str, dict]:
        for index, event in enumerate(self.model["events"]):
            if isinstance(event, dict):
                self.dictionary_event(event, index)
            else:
                self.list_event(event, index)
        if len(self.noise_locations) != len(self.model["channels"]):
            raise ValueError("Noise location coverage differs from the model")
        if self.faults is not None and not self.faults.keys() <= set(
            range(len(self.noise_locations))
        ):
            raise ValueError("Unknown fixed-fault location")
        if "observables" in self.model:
            for detector in self.model["detectors"]:
                self.detector(detector)
            if len(self.model["observables"]) != 1:
                raise ValueError("Expected one logical observable")
            self.lines.append(
                "OBSERVABLE_INCLUDE(0) " + " ".join(map(self.rec, self.model["observables"][0]))
            )
        elif "observable_records" in self.model:
            self.lines.append(
                "OBSERVABLE_INCLUDE(0) " + " ".join(map(self.rec, self.model["observable_records"]))
            )
            self.expected_observables = [self.model["observable_reference_parity"]]
        else:
            # The references evaluate (I - (X_L + Y_L)/sqrt(2))/2 after acceptance.
            # This noiseless logical T-dagger followed by X_L has that same POVM.
            if "logical" in self.model:
                x, z = {self.model["logical"]: "X"}, {self.model["logical"]: "Z"}
            else:
                x = {int(q): axis for q, axis in self.model["logical_x"].items()}
                z = {int(q): axis for q, axis in self.model["logical_z"].items()}
            self.lines.append("R_PAULI(-0.25) " + pauli_text(z))
            self.measure("MPP " + pauli_text(x))
            self.lines.append("OBSERVABLE_INCLUDE(0) rec[-1]")
        metadata = dict(
            measurements=self.measurements,
            expected_detectors=self.expected_detectors,
            expected_observables=self.expected_observables,
            noise_locations=self.noise_locations,
            source_record_map=self.records,
        )
        return "\n".join(self.lines) + "\n", metadata


def export_corpus(root: Path, output: Path, p: float) -> dict:
    verify_source(root)
    output.mkdir(parents=True, exist_ok=False)
    manifest: dict = dict(repository=REPOSITORY, revision=REVISION, p=p, cases={})
    for case in CASES:
        source = root / "data/alternatives/inputs" / case / "model.json"
        text, metadata = Exporter(json.loads(source.read_text()), p).export()
        path = output / (case + ".stim")
        path.write_text(text)
        metadata.update(
            model_path=str(source.relative_to(root)),
            model_sha256=digest(source),
            circuit_sha256=digest(path),
            reference=json.loads(
                (root / "data/alternatives/series" / (case + ".json")).read_text()
            ),
        )
        manifest["cases"][case] = metadata
    source = root / "data/inputs/soft_cultivation_d5_p0005.stim"
    text = re.sub(
        r"((?:DEPOLARIZE[123]|[XYZ]_ERROR|M[XYZ]?|MR[XYZ]?)\()[^)]*(\))",
        lambda m: m[1] + f"{p:.17g}" + m[2],
        source.read_text(),
    )
    path = output / "soft-d5.stim"
    path.write_text(text)
    manifest["cases"]["soft-d5"] = dict(
        source_path=str(source.relative_to(root)),
        source_sha256=digest(source),
        circuit_sha256=digest(path),
        expected_detectors=[],
        expected_observables=[0],
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--p", type=float, default=0.001)
    args = parser.parse_args()
    if not 0 <= args.p <= 1:
        parser.error("p must be a probability")
    manifest = export_corpus(args.source, args.output, args.p)
    print(
        json.dumps(
            {"revision": REVISION, "cases": list(manifest["cases"]), "output": str(args.output)}
        )
    )


if __name__ == "__main__":
    main()
