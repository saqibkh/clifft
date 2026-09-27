# Larger-distance folded controls

2026-09-27. This study continues without waiting for larger circuits from the
authors. There are two separate experiments: a **complete reconstructed f7
protocol**, and **constructed terminal gadgets at code distances 3, 5, 7 and 9**.
Neither is an author-supplied larger-distance export.

Current main was checked remotely at
`63fc77ab7753220f396a87c274e9d2252a4b1b71`. The full f7 sampler links that main's
core library and enables the opt-in active-width scheduler with unlimited
search. The constructed controls use the same main for physical probability
references at d3/d5, and the existing experimental native contraction kernel
for gadget sampling. These are research programs, not a new production pass.

## Complete reconstructed f7

The measured median is **172 ms per complete attempt**, over three trials of
32 attempts at `p=0.001` on CPU 1. Each attempt produces 351 visible records and
436 hidden reset records. [Full f7 evidence](reconstructed_f7.json) contains the
trials, circuit hash, physical fault histories and independent replay results.

Twelve complete histories were checked at `p=0`, `0.001` and `0.03`, including
histories with 50--56 physical faults. The independent coherent-Clifford
physical oracle (Stim/Cirq tableaux) agrees with the sampler's conditional log
probabilities to `3.56e-15`; physical current-main Clifft supplies the prefix
probability. Detector and logical bits also agree with the full records.

The input and older research compiler/sampler are pinned to
`b87c0da396b4bbde54fbee5535aeac1116170605`. The new driver checks the reference
source bytes against that commit before using them, builds the native sampler
and physical-history replay executable against current main, and regenerates
the circuit and compiled plans.

This circuit includes injection, the smaller folded checks, 5-to-7 growth,
both d7 checks, noisy CSS extraction, final syndrome extraction and logical
readout. The run retains every visible and hidden measurement record on both
accepted and rejected attempts. It uses a current-main width-eight prefix,
fixed contractions for the folded blocks, and the existing certified growth
transfer. The terminal sampler uses the older row ordering; it has not received
the syndrome-order optimization used for the new constructed controls.

The f7 cat preparation contains an additional verification GHZ state. The
narrow Sahay model recognizer correctly does not accept that three-Hadamard
preparation/readout pattern. The older f7 compiler already eliminates the extra
coherent bit by conditioning on its verification measurement and proves the
remaining four-path structure. The experiment uses that existing certificate;
it does not weaken the Sahay recognizer's guard.

Full dense execution remains outside the allocation budget: even with the
opt-in scheduler, current main reaches active width 44, implying 256 TiB for
complex-double coefficients alone. Only its plan is constructed. Consequently
there is no measured dense-versus-specialized f7 speedup ratio.

## Constructed terminal gadgets

The measured medians below are from three trials of 16 samples per distance,
at `p=0.001`, with early rejection disabled. All four cases use CPU 0.
[Raw evidence](constructed_folded.json) includes each trial and source hashes.

| Constructed gadget | Data qubits | Output records | Median per gadget | Numeric plans and workspace |
| --- | ---: | ---: | ---: | ---: |
| d3 | 13 | 21 | 1.92 ms | 0.019 MiB |
| d5 | 41 | 53 | 7.81 ms | 0.45 MiB |
| d7 | 85 | 101 | 81.5 ms | 6.60 MiB |
| d9 | 145 | 165 | 560 ms | 28.66 MiB |

The memory column counts contraction payload, not process RSS, compiler memory
or Python object overhead. It is derived from the fixed serialized plans.
The largest d9 contraction has 8,192 entries. These timings cannot be used as
speedup ratios against either the different full f7 circuit or external d5.

[The generator](../constructed_folded.py) uses the pinned regular surface-code
geometry, alternating native controlled-SX factors along the fold, and CCZs
between reflected data pairs. Each of two checks prepares a chain GHZ cat,
measures a two-endpoint parity verifier, applies the folded factors, uncomputes
the cat, and measures all its qubits. Noiseless CSS projection and a logical
magic-axis measurement finish the instrument.

Noise is explicit: independent probability-p reset/measurement X flips and
uniform nonidentity one-, two- or three-qubit Pauli faults after each gate.
There is no idle-noise schedule. A native controlled gate has one noise site;
its exact T decomposition does not create extra sites. The cat schedule is a
simple constructed control, **not a validated fault-tolerant d7/d9 schedule**.
These models differ from both the external Sahay d5 model and the reconstructed
full f7 circuit.

All distances enter from a supplied signed encoded logical state. They include
both checks, all physical faults, every body and CSS record, and the final
logical readout. They exclude injection, growth, noisy CSS extraction, escape
and decoding. Timings include Python physical-fault binding and sampling plus
the native fixed-plan arithmetic. These are full **terminal-gadget** timings,
not full cultivation attempts and not the earlier one-contraction-per-prefix
kernel timings.

The model-level structural certificate uses actual gates and Pauli supports
after selecting candidates from native-block annotations. It proves
branch-independent records/resets and compatible data-code translations for
arbitrary physical fault bits. All controls have four coherent paths, including
the 145-data-qubit, 72-X-generator d9 case. The amplitude elimination order is
used for syndrome sampling; physical record order remains unchanged.

## Correctness checks

- At every distance, six logical states and all four cat-outcome pairs are
  checked against the ideal logical magic-axis projectors. This checks that
  the constructed larger gadgets implement the intended encoded operation.
- At d3/d5, physical current-main Clifft computes complete record probabilities
  with logical X, Y, Z and magic-axis endings: 96 comparisons at each distance,
  including random signed boundaries and multi-qubit faults.
- At every distance, a separate offline quadratic Gauss-sum implementation
  checks amplitudes and selected Born marginals for six boundary/fault cases:
  144 comparisons per distance. This oracle uses algebraic elimination rather
  than the tested tensor-contraction algorithm. It shares the physical local
  factor binding, so it is an independent arithmetic reference, not a second
  independent d9 physical-circuit implementation. Its polynomial elimination
  is separately checked against 216 exhaustive small binary sums.
- Relative errors are checked for noncancelling marginals, so tiny complete
  record probabilities cannot pass merely because an absolute tolerance is
  larger than the answer.

The largest d9 relative discrepancy is `1.89e-15`; positive probabilities down
to `3.47e-18` were checked. The ideal projector comparisons are within
`1.34e-15` at all distances, and the d3/d5 physical reference comparisons are
within `5.00e-16`.

Large signed/faulted d9 cases exposed an overly strong assertion in the earlier
research sampling helper: it rejected a valid probability of `2^-54` even
though the amplitude norm and marginal agreed. The helper now requires a
positive probability and relative agreement. It no longer assumes individual
records have probability at least `1e-14`. No production behavior changes.

The external Sahay d5 validation was rerun with the stricter relative guard:
all 576 region-probability comparisons and 15 complete-record comparisons still
pass (maximum absolute region error `1.45e-15`, maximum relative full-record
error `1.34e-15`). See [the regression evidence](sahay_gadget_guard_regression.json).
The seven existing adapter tests against Aer/Stim also pass.

## Reproduction

Use an isolated Python environment with NumPy, Stim and current-main Clifft.
The older full-f7 physical-history oracle additionally requires Cirq. A separate
interpreter avoids a shared virtual environment's editable Clifft import hook.

Extract the pinned reference files into a separate directory, then run:

```sh
python tools/profile/study_constructed_folded.py \
  --reference /path/to/pinned-reference \
  --work /tmp/constructed-folded \
  --output /tmp/constructed-folded.json \
  --distances 3 5 7 9 --shots 16

python tools/profile/study_reconstructed_f7.py \
  --reference /path/to/pinned-reference \
  --core-library /path/to/current-main-build/src/clifft/libclifft_core.a \
  --work /tmp/reconstructed-f7 \
  --output /tmp/reconstructed-f7.json
```

The work directories must be new. Source, circuit, library and executable
hashes are recorded in the evidence. Setup/build/plan generation is excluded
from sampling timings and recorded separately. The measurements use scalar
sampling on one EPYC 9554P core and three timing trials. Small shot counts are
performance checks, not logical-error-rate estimates.

## What this supports next

Additional author circuits are not required to investigate scaling. We now
have a full f7 control and a d9 terminal instrument with explicit provenance,
records and probability checks. The remaining practical work is to unify the
guarded recognition and continuation path for ordinary Clifft inputs and move
research binding into compiler-precomputed, preallocated native execution.
A full d9 cultivation claim would additionally require constructing and
validating the 7-to-9 growth/entry handoff and its noise schedule. The older
full-protocol harness's 128-bit physical masks and 64-bit syndrome masks also
need deliberate extension; the constructed Python-bound d9 experiment does
not exercise those limits.
