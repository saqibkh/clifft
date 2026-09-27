"""Offline instrument reference for the supplied Sahay folded region.

Recognize a CSS boundary and a fixed monomial path structure. Physical faults
change local values and record bits, never contraction scopes. This is a
research reference, not a production Clifft executor or recognizer.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from study_wan_scaling import masks_for, permute_masks, scope_plan
from wan_corpus import Exporter, gate_lines, pauli_text

ROOTS = np.exp(1j * np.pi / 4 * np.arange(8))


def parity(value: int) -> int:
    return value.bit_count() & 1


@dataclass
class Boundary:
    offset: int
    signs: int
    logical: Any
    ancillas: int


@dataclass
class Bound:
    coefficients: Any
    values: Any
    offset: int
    records: tuple[int, ...]
    ancillas: int


class SahayGadget:
    def __init__(self, model: dict, reference: Any):
        self.model = model
        events = model["events"]
        checks = [e for e in events if e["name"] == "MEASURE_PRODUCT"]
        self.data = sorted({q for e in checks for q in e["support"]})
        self.n = len(self.data)
        self.index = {q: i for i, q in enumerate(self.data)}
        self.checks = checks
        self.x_rows, self.z_rows = [
            [self.mask(e["support"]) for e in checks if e["axis"] == axis] for axis in "XZ"
        ]
        self.rank = len(self.x_rows)
        if len(self.z_rows) != self.rank or self.n != 2 * self.rank + 1:
            raise ValueError("Expected a one-logical-qubit CSS code")
        for axis in "XZ":
            logical = model["logical_" + axis.lower()]
            if not logical or any(value != axis for value in logical.values()):
                raise ValueError("Expected pure CSS logical Pauli supports")
        self.logical_x, self.logical_z = [
            self.mask(int(q) for q, a in model["logical_" + axis.lower()].items() if a == axis)
            for axis in "XZ"
        ]
        if any(parity(a & b) for a in self.x_rows for b in self.z_rows):
            raise ValueError("Noncommuting code checks")
        if not parity(self.logical_x & self.logical_z) or any(
            parity(row & logical)
            for rows, logical in ((self.x_rows, self.logical_z), (self.z_rows, self.logical_x))
            for row in rows
        ):
            raise ValueError("Incorrect logical Pauli algebra")
        self.x_duals = reference.binary_duals(self.x_rows + [self.logical_x], self.n)
        self.z_duals = reference.binary_duals(self.z_rows + [self.logical_z], self.n)
        native = sorted({e["native_block"] for e in events if "native_block" in e})[-2:]
        if len(native) != 2:
            raise ValueError("Expected two terminal folded checks")
        first = next(i for i, e in enumerate(events) if e.get("native_block") == native[0])
        opening = max(i for i, e in enumerate(events[:first]) if e["name"] == "H")
        self.start = 1 + max(i for i, e in enumerate(events[:opening]) if e["name"] == "R")
        self.stop = next(i for i, e in enumerate(events) if e["name"] == "MEASURE_PRODUCT")
        if any(e["name"] not in {"MEASURE_PRODUCT", "DETECTOR"} for e in events[self.stop :]):
            raise ValueError("Terminal projection is not noiseless")
        used = {q for e in events[self.start : self.stop] for q in e.get("support", [])}
        self.noise = [e["location"] for e in events[self.start : self.stop] if e["name"] == "NOISE"]
        for location in self.noise:
            used.update(model["channels"][location]["support"])
        self.ancillas = sorted(used - set(self.data))
        self.index.update({q: self.n + i for i, q in enumerate(self.ancillas)})
        pairs = sorted(
            {
                tuple(sorted(self.index[q] for q in e["support"][1:]))
                for e in events[self.start : self.stop]
                if e["name"] == "CCZ"
            }
        )
        self.pairs = pairs
        if any(a == b or a >= self.n or b >= self.n for a, b in pairs):
            raise ValueError("Pair factors must act on distinct data qubits")
        self.pair_index = {pair: self.n + i for i, pair in enumerate(pairs)}
        masks = masks_for(self.x_rows, self.n, pairs)
        self.characters = len(masks) - self.rank
        self.order = scope_plan(masks, self.rank)["order"]
        self.masks = permute_masks(masks, self.order)
        self.operations: list[tuple[str, tuple[int, ...], int]] = []
        self.record_ids = []
        for event in events[self.start : self.stop]:
            name = event["name"]
            support = tuple(self.index[q] for q in event.get("support", []))
            if name == "DETECTOR":
                continue
            if name == "NOISE":
                self.operations.append((name, (), event["location"]))
            elif name == "MEASURE":
                if event["axis"] != "Z" or len(support) != 1:
                    raise ValueError("Body measurements must be computational")
                self.operations.append((name, support, int(event["inverted"])))
                self.record_ids.append(event["record"])
            else:
                for line in gate_lines(name, list(support)):
                    gate, *qubits = line.split()
                    if len(qubits) != {"CX": 2, "CCZ": 3}.get(gate, 1):
                        raise ValueError("Unsupported gate arity")
                    self.operations.append((gate, tuple(map(int, qubits)), 0))
        self.shifts = self.certify_paths()

    def mask(self, support: Any) -> int:
        return sum(1 << self.index[q] for q in support)

    def certify_paths(self) -> list[int]:
        """Prove branch-independent records for arbitrary input and fault bits."""
        # Two coherent cat bits, two closing outcomes, then independent entry
        # and physical-fault X bits. Each integer encodes an affine expression.
        bits = [1 << (4 + i) for i in range(len(self.ancillas))]
        symbol = 4 + len(bits)
        flips = [0] * self.n
        roots = []
        hadamards = 0
        for gate, qubits, location in self.operations:
            if gate == "NOISE":
                for target in self.model["channels"][location]["support"]:
                    q = self.index[target]
                    if q < self.n:
                        flips[q] ^= 1 << symbol
                    else:
                        bits[q - self.n] ^= 1 << symbol
                    symbol += 1
                continue
            q = qubits[0]
            a = q - self.n
            if gate in {"R", "MEASURE"}:
                if a < 0 or bits[a] & 3:
                    raise ValueError("Coherent branch leaks into a reset or measurement")
                if gate == "R":
                    bits[a] = 0
            elif gate == "H":
                if a < 0 or hadamards >= 4:
                    raise ValueError("Expected two cat preparation/readout pairs")
                roots.append(q)
                closing = hadamards & 1
                if not closing and bits[a] & 3:
                    raise ValueError("Cat preparation has a coherent input")
                bits[a] = 1 << (hadamards // 2 + 2 * closing)
                hadamards += 1
            elif gate == "CX":
                if a < 0:
                    raise ValueError("Expected classical-path ancilla control")
                target = qubits[1]
                if target < self.n:
                    flips[target] ^= bits[a]
                else:
                    bits[target - self.n] ^= bits[a]
            elif gate == "CCZ":
                if a < 0 or tuple(sorted(qubits[1:])) not in self.pair_index:
                    raise ValueError("Expected a fixed data-pair factor")
            elif gate not in {"T", "T_DAG", "S", "S_DAG", "Z"}:
                raise ValueError("Unsupported monomial path gate " + gate)
        if hadamards != 4 or roots[0] != roots[1] or roots[2] != roots[3]:
            raise ValueError("Incomplete cat preparation/readout")
        if any(value & 3 for value in bits):
            raise ValueError("Unmeasured ancilla retains a coherent branch")
        shifts = []
        for branch in itertools.product((0, 1), repeat=2):
            mask = branch[0] | branch[1] << 1
            delta = sum(parity(value & mask) << q for q, value in enumerate(flips))
            shift = sum(parity(delta & dual) << i for i, dual in enumerate(self.x_duals))
            restored = 0
            for i, row in enumerate(self.x_rows + [self.logical_x]):
                if shift >> i & 1:
                    restored ^= row
            if delta != restored:
                raise ValueError("Coherent paths end in different computational code cosets")
            shifts.append(shift)
        self.certificate = dict(
            symbolic_inputs=symbol,
            coherent_paths=4,
            hadamards=hadamards,
            classical_records=len(self.record_ids),
            branch_independent_records_and_resets=True,
            same_code_coset=True,
            data_qubits=self.n,
            x_generators=self.rank,
            data_pair_factors=len(self.pairs),
            physical_noise_sites=len(self.noise),
            prefix_event_stop=self.start,
            projection_event_start=self.stop,
        )
        return shifts

    def prefix(self, p: float, faults: dict[int, int] | None = None) -> tuple[str, dict]:
        import stim

        exporter = Exporter(self.model, p, faults)
        for i, event in enumerate(self.model["events"][: self.start]):
            exporter.dictionary_event(event, i)
        probes = [
            dict.fromkeys(e["support"], axis)
            for axis in "XZ"
            for e in self.checks
            if e["axis"] == axis
        ]
        width = max(self.index) + 1
        x, z = stim.PauliString(width), stim.PauliString(width)
        for q, axis in self.model["logical_x"].items():
            x[int(q)] = axis
        for q, axis in self.model["logical_z"].items():
            z[int(q)] = axis
        y = 1j * x * z
        if y.sign != 1:
            raise ValueError("Unexpected physical logical-Y phase")
        probes += [{q: "IXYZ"[a] for q, a in enumerate(op) if a} for op in (x, y, z)]
        probes += [{q: "Z"} for q in self.ancillas]
        text = (
            "\n".join(exporter.lines + ["EXP_VAL " + pauli_text(probe) for probe in probes]) + "\n"
        )
        return text, dict(records=exporter.records, expected_detectors=exporter.expected_detectors)

    def boundary(self, values: Any) -> Boundary:
        rank = self.rank
        classical = np.concatenate([values[: 2 * rank], values[2 * rank + 3 :]])
        if np.max(np.abs(np.abs(classical) - 1)) > 1e-10:
            raise ValueError("Prefix is not at the certified classical/code boundary")
        x, y, z = values[2 * rank : 2 * rank + 3]
        if abs(x * x + y * y + z * z - 1) > 1e-10:
            raise ValueError("Logical input is entangled with unrepresented spectators")
        alpha = math.sqrt(max(0.0, (1 + z) / 2))
        logical = (
            np.array([alpha, (x + 1j * y) / (2 * alpha)])
            if alpha > 1e-10
            else np.array([0, 1], dtype=complex)
        )
        offset = 0
        for i, dual in enumerate(self.z_duals[:-1]):
            if values[rank + i] < 0:
                offset ^= dual
        return Boundary(
            offset,
            sum(int(v < 0) << i for i, v in enumerate(values[:rank])),
            logical,
            sum(int(v < 0) << i for i, v in enumerate(values[2 * rank + 3 :])),
        )

    def path(
        self,
        boundary: Boundary,
        faults: dict[int, int],
        branches: tuple[int, ...],
        outcomes: tuple[int, ...],
    ) -> tuple:
        local = np.ones((self.characters, 4), dtype=complex)
        bits, flips, scalar = boundary.ancillas, 0, 1 + 0j
        records = []
        hadamards = 0

        def phase(q: int, exponent: int) -> None:
            nonlocal scalar
            if q < self.n:
                flip = (flips >> q) & 1
                local[q, 0] *= ROOTS[(exponent * flip) % 8]
                local[q, 1] *= ROOTS[(exponent * (flip ^ 1)) % 8]
            else:
                scalar *= ROOTS[(exponent * ((bits >> (q - self.n)) & 1)) % 8]

        def fault(q: int, axis: str) -> None:
            nonlocal scalar, bits, flips
            if axis in {"Y", "Z"}:
                phase(q, 4)
            if axis == "Y":
                scalar *= 1j
            if axis in {"X", "Y"}:
                if q < self.n:
                    flips ^= 1 << q
                else:
                    bits ^= 1 << (q - self.n)

        for gate, qubits, location in self.operations:
            if gate == "NOISE":
                code = faults.get(location, 0)
                if code:
                    channel = self.model["channels"][location]
                    for j, target in enumerate(channel["support"]):
                        axis = (
                            channel["axis"]
                            if channel["kind"] == "flip"
                            else "IXZY"[(code >> (2 * j)) & 3]
                        )
                        fault(self.index[target], axis)
                continue
            q = qubits[0]
            a = q - self.n
            if gate in {"T", "T_DAG", "S", "S_DAG", "Z"}:
                phase(q, {"T": 1, "T_DAG": -1, "S": 2, "S_DAG": -2, "Z": 4}[gate])
            elif gate == "CX":
                target = qubits[1]
                if bits >> a & 1:
                    if target < self.n:
                        flips ^= 1 << target
                    else:
                        bits ^= 1 << (target - self.n)
            elif gate == "CCZ":
                x, y = sorted(qubits[1:])
                if bits >> a & 1:
                    table = local[self.pair_index[x, y]]
                    table[3 ^ (((flips >> x) & 1) | (((flips >> y) & 1) << 1))] *= -1
            elif gate == "R":
                bits &= ~(1 << a)
            elif gate == "H":
                bit = (outcomes if hadamards & 1 else branches)[hadamards // 2]
                scalar *= (-1 if bits >> a & bit & 1 else 1) / math.sqrt(2)
                bits = (bits & ~(1 << a)) | (bit << a)
                hadamards += 1
            elif gate == "MEASURE":
                records.append(((bits >> a) & 1) ^ location)
        return scalar, flips, local, tuple(records), bits

    def bind(self, boundary: Boundary, faults: dict[int, int], outcomes: tuple[int, ...]) -> Bound:
        paths = [
            self.path(boundary, faults, b, outcomes) for b in itertools.product((0, 1), repeat=2)
        ]
        if any(path[3:] != paths[0][3:] for path in paths):
            raise AssertionError("Certified branch-independent outputs changed")
        offset = boundary.offset ^ paths[0][1]
        coefficients = np.empty((2, 4), dtype=complex)
        values = np.ones((2, 4, len(self.masks), 4), dtype=complex)
        for k, (scalar, flips, local, _, _) in enumerate(paths):
            shift = self.shifts[k]
            for logical in (0, 1):
                base = offset ^ flips ^ (self.logical_x if logical else 0)
                target = values[logical, k]
                for q in range(self.n):
                    bit = (base >> q) & 1
                    target[q, :2] = local[q, [bit, bit ^ 1]]
                for i, (a, b) in enumerate(self.pairs):
                    target[self.n + i] = local[
                        self.n + i, np.arange(4) ^ (((base >> a) & 1) | (((base >> b) & 1) << 1))
                    ]
                for i in range(self.rank):
                    target[self.characters + i, 1] = (-1) ** ((boundary.signs >> i) & 1)
                coefficients[logical, k] = (
                    scalar
                    * boundary.logical[logical ^ (shift >> self.rank)]
                    * (-1) ** parity(boundary.signs & shift)
                )
        return Bound(coefficients, values, offset, paths[0][3], paths[0][4])

    def logical_bloch(self, bound: Bound, amplitudes: Any) -> Any:
        a, b = amplitudes
        norm = float(np.vdot(amplitudes, amplitudes).real)
        cross = a.conjugate() * b
        sign = (-1) ** parity(bound.offset & self.logical_z)
        return (
            np.array([2 * cross.real, sign * 2 * cross.imag, sign * (abs(a) ** 2 - abs(b) ** 2)])
            / norm
        )
