"""Independent checks for the research model adapter's nontrivial boundaries."""

import numpy as np
import pytest
import stim
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from study_wan_corpus import fresh_reset_wires
from wan_corpus import Exporter, gate_lines

import clifft


@pytest.mark.parametrize("name", ["CSX", "CS_DAG_X", "INJECT_HXY"])
def test_native_gate_on_entangled_input(name):
    if name == "INJECT_HXY":
        w = np.exp(1j * np.pi / 4)
        matrix = np.array([[1, w**3], [w, 1]]) / np.sqrt(2)
        support = [0]
    else:
        matrix = np.eye(4, dtype=complex)
        matrix[1, 1] = matrix[3, 3] = 0
        matrix[1, 3] = 1
        matrix[3, 1] = 1j if name == "CSX" else -1j
        support = [0, 1]
    n = len(support)
    qc = QuantumCircuit(2 * n)
    prefix = []
    for q in support:
        qc.h(q)
        qc.cx(q, q + n)
        prefix.extend([f"H {q}", f"CX {q} {q + n}"])
    qc.unitary(matrix, support)
    qc.save_statevector()
    reference = np.asarray(AerSimulator(method="statevector").run(qc).result().get_statevector())
    program = clifft.compile("\n".join(prefix + gate_lines(name, support)))
    actual = clifft.get_statevector(program)
    # Choi states check the entire controlled map, including relative phases.
    assert abs(np.vdot(reference, actual)) == pytest.approx(1, abs=1e-12)


def test_reported_fault_preserves_postmeasurement_state():
    model = dict(
        channels=[dict(kind="flip", support=[0], source_kind="READOUT_FLIP", record=0)],
        events=[["X", 0], ["MEASURE", 0, "Z", 0], ["NOISE", 0], ["MEASURE", 0, "Z", 1]],
        detectors=[[0, 1]],
        observables=[[0]],
    )
    text, _ = Exporter(model, 0, {0: 1}).export()
    actual = clifft.sample(clifft.compile(text), 16, seed=12)
    reference = stim.Circuit("X 0\nM(1) 0\nM 0").compile_sampler(seed=12).sample(16)
    np.testing.assert_array_equal(actual.measurements, reference)
    np.testing.assert_array_equal(actual.measurements, np.tile([0, 1], (16, 1)))


def test_detector_reference_and_record_feedback():
    model = dict(
        channels=[],
        events=[
            dict(name="X", support=[0]),
            dict(name="MEASURE", support=[0], axis="Z", inverted=False, record=0),
            dict(name="FEEDBACK_X", support=[1], record=0),
            dict(name="MEASURE", support=[1], axis="Z", inverted=False, record=1),
            dict(name="DETECTOR", records=[1], reference_parity=1),
        ],
        observable_records=[0, 1],
        observable_reference_parity=0,
    )
    text, meta = Exporter(model, 0).export()
    program = clifft.compile(text, expected_detectors=meta["expected_detectors"])
    actual = clifft.sample(program, 16, seed=12)
    reference = stim.Circuit("X 0\nM 0\nCX rec[-1] 1\nM 1").compile_sampler().sample(16)
    np.testing.assert_array_equal(actual.measurements, reference)
    assert not actual.detectors.any()
    assert not actual.observables.any()


def test_logical_readout_measures_magic_axis():
    model: dict = dict(
        channels=[], logical=0, events=[dict(name="H", support=[0]), dict(name="T", support=[0])]
    )
    text, _ = Exporter(model, 0).export()
    assert clifft.record_probabilities(
        clifft.compile(fresh_reset_wires(text)), ["0", "1"]
    ) == pytest.approx([1, 0])
    model["events"].append(dict(name="Z", support=[0]))
    text, _ = Exporter(model, 0).export()
    assert clifft.record_probabilities(
        clifft.compile(fresh_reset_wires(text)), ["0", "1"]
    ) == pytest.approx([0, 1])


def test_fresh_reset_wires_preserve_entangled_marginal():
    text = fresh_reset_wires("H 0\nCX 0 1\nR 0\nM 0\nM 1\n")
    actual = clifft.record_probabilities(clifft.compile(text), ["00", "01", "10", "11"])
    assert actual == pytest.approx([0.5, 0.5, 0, 0], abs=1e-12)
