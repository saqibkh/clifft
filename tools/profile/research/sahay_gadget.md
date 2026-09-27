# Complete Sahay d5 gadget experiment

2026-09-27. **The supplied Sahay Fold d5 model now has a working experimental
replacement for its two large checks and terminal code projection.** Ordinary
Clifft runs the original prefix; the replacement samples the remaining physical
records and logical readout. This is a full noisy-attempt experiment, unlike
the earlier contraction-kernel benchmark.

Main was rechecked at `63fc77ab7753220f396a87c274e9d2252a4b1b71`. The input is
`data/alternatives/inputs/fold-d5/model.json` from Wan's repository at
`b99452de7e300324e7bbf5cd9347241c3e289bae`. Its hash and the current-main extension
hash are in [the evidence](sahay_gadget.json). This remains the externally
supplied **native-T reconstruction**, not a matched author-supplied noisy T
export. It includes injection, growth and cultivation but excludes escape and
decoding. The [corpus report](wan_corpus.md) describes those qualifications.

## What runs

The full circuit has 123 visible measurement records, 118 detectors and 2,116
physical noise locations. The prefix retains its first 70 records and 1,418
noise locations. The gadget includes:

- Both cat preparations and their verification records.
- Both large native checks, including control phases, controlled-gate faults
  and idle noise.
- All 12 body records, including the supplied reported-bit noise dilation.
- All 40 terminal data-code measurements and the final logical readout.
- All 698 remaining physical noise locations.

Neither simulator exits early when a detector rejects an attempt in this
comparison. Both produce every visible record, including rejected attempts.
The gadget also computes the conditional output logical amplitudes, their
probability, and the signed code frame. That representation is checked against
physical logical X/Y/Z measurements. It is not yet restored into a production
Clifft executor for arbitrary further circuit execution.

## The structure that was recognized

Candidate region boundaries are selected using the input model's native-block
annotations. The subsequent algebraic checks use the actual operations and
Pauli supports, not the block names or a circuit hash. This is an offline
model-level recognizer; it is not a general pattern matcher for arbitrary
Clifft syntax.

The entry boundary is moved before the first cat preparation, after its initial
resets: event 2,275. This avoids assuming an already prepared GHZ state is a
set of classical ancillas. Current Clifft's prefix plan reduces all 40 data
stabilizers and seven required ancilla Z operators to tracked signs. The three
physical logical Pauli operators act on its one active coordinate. All 1,024
noisy prefix samples satisfy those conditions; the largest logical-purity
error is `8.88e-16`. The prefix's peak width is eight; its exit width is one.

Within the region, fixing the two cat branches makes the data evolution a
permutation of computational basis states with local phases. The compiler
checks:

1. The terminal checks define an independent, commuting one-logical-qubit CSS
   code, with valid logical X/Z operators.
2. The body has two complete cat preparation/readout pairs. Data phases and
   pair interactions fit fixed one- or two-parity factors.
3. Every visible measurement and discarded reset bit is independent of the
   unresolved coherent branches. This is checked symbolically for arbitrary
   entry-ancilla and physical-fault X bits, not only on sampled faults.
4. The differences between coherent paths' final data flips belong to the
   span of the X stabilizers and logical X. Their code translations and
   associated sign rules are computed before sampling.
5. The terminal projection is noiseless; the supplied measurement-record
   ordering is retained.

There are **four coherent paths**, 41 data qubits, 20 X generators and 16
pair factors. The symbolic check uses 811 input/branch/outcome/fault variables.
After a physical fault pattern is drawn, its Pauli effects change numerical
factor values and record bits. They do not change contraction topology.

The evaluator first samples the joint outcomes of the two cat readouts. It
then samples the terminal X syndrome with the improved offline ordering from
the [scaling study](wan_scaling.md), computes the two output logical amplitudes,
and samples the supplied logical readout. Z syndromes follow from the signed
computational code frame. Coherent paths are summed with their relative phases;
they are not sampled as classical alternatives.

## Correctness checks

The main comparison uses ordinary Clifft's exact measurement-record probability
API on the physical circuit. Fixed Pauli faults replace stochastic channels for
these oracle queries. Resets are purified with fresh zero wires and old wires
are traced, using the already tested oracle adapter; the timed simulations
retain the original physical resets and wire count.

- **18 region cases, 576 probability comparisons.** Six logical input states
  include computational, X/Y, magic and generic complex states. Cases include
  random signed code spaces and entry ancilla values, multiple physical faults,
  and targeted correlated two-/three-qubit faults and reported-bit flips.
  Each selected body/projection record is checked with both outcomes of
  physical logical X, Y, Z and magic-axis readout. The maximum absolute error
  is `1.45e-15`. These are selected record queries, with some repetition, not
  exhaustive enumeration of the d5 record space.
- **15 complete-record comparisons across five full-circuit fault patterns.**
  The original prefix produces a boundary and visible prefix record. The
  product of its exact record probability and the gadget's conditional suffix
  probability matches ordinary Clifft on the complete 123-bit record. The
  maximum relative error is `1.34e-15`. This checks the connection to the actual
  injection/growth prefix, rather than only an ideal encoded input.
- **Four malformed structures rejected.** Measuring an unresolved coherent
  branch, using a data qubit as a pair-factor control, changing a required data
  translation, or deleting a closing Hadamard each fails the algebraic checks.
- **Independent primitive checks retained.** All seven existing adapter tests
  pass, including Qiskit Aer checks of the controlled native gates on entangled
  inputs, Stim checks of reported measurement flips, and reset/readout checks.

These checks support the complete physical instrument on this pinned model.
They do not constitute a formal machine-checked proof for all possible models
or an exhaustive enumeration of all noise histories. The original circuit-level
experimental recognizer still falls back on this input; no guard was removed to force it to
select the gadget.

## Full-attempt timing

| Method on the supplied Sahay d5 model | Median ms per complete attempt | Trial range |
| --- | ---: | ---: |
| Current main with unlimited active-width scheduling | 85.97 | 84.79-88.88 |
| Original Clifft prefix plus experimental gadget | 8.56 | 8.40-9.06 |

This is approximately **10x faster** for the measured scalar CPU configuration.
It includes prefix execution, all physical noise draws, coherent interference,
syndrome sampling, terminal readout, record assembly and detector evaluation.
The gadget's Python binding overhead is included. The earlier 0.131 ms
contraction-only number is not used as an attempt time.

Both methods use p=0.001, one worker, scalar batches, all records retained, and
CPU 0 on the EPYC 9554P VM. Five interleaved trials of 64 attempts give 320
attempts per method. Compilation and two warmup attempts per method are outside
timing. The ordinary extension is the existing GCC 13.3 Release/native build;
the new arithmetic wrapper uses GCC 13.3, `-O3 -march=native -fno-fast-math`.

The ordinary full plan takes 0.526 s to compile; the prefix takes 0.017 s.
Experimental model/plan preparation plus invoking the C++ compiler to build
the arithmetic wrapper takes 1.62 s. That setup number is not a production
compile-time estimate. The numeric contraction family is the same approximately
0.46 MiB plan/workspace payload measured previously, plus bound values and
Python/prefix storage; this is not a measured whole-process memory figure.

The methods use different random streams. Their records are therefore not
expected to match shot for shot. They accept 53/320 and 43/320 attempts
respectively, with no observed accepted logical errors; these small samples
neither establish rare-error rates nor replace the exact probability checks.
The comparison does not establish the best possible batched, threaded or GPU
performance, and postselected throughput was not measured in this experiment.

## Consequence for larger distances

We now have an external model demonstrating that the gadget idea extends beyond
our reconstructed fixtures. This also improves d5 time in the tested setup,
although that was not required for the structural experiment to succeed.

The algebraic conditions above are not tied to distance five or particular
qubit identifiers. An equivalent two-check circuit can still have four coherent
paths as its data code grows. However, the sums within each path still require
contraction. The earlier d7/d9 size results remain encouraging; the d11 planner
limits remain relevant. Four paths alone do not establish polynomial scaling.

**No new complete d7/d9 physical circuit was benchmarked here.** Transferring the
result requires a candidate with a certified entry boundary, compatible gate
and measurement structure, correct physical fault binding and affordable
contraction plans. A larger geometric code with no validated injection/growth
protocol is not evidence of a complete fault-tolerant protocol at that distance.

The next implementation step is to port this verified model-level binding and
record/state interface into the experimental native gadget path, with all
storage preallocated and planning outside execution. Then test a larger
physical control against these same conditions. The Python research driver is
an allocating reference harness, not a production hot executor. No production
source, public API, Stim source or architecture was changed in this study.

## Reproduction

Use the pinned external model and isolated current-main package from the prior
reports. Extract the existing experimental contraction reference without
switching the working branch:

```bash
mkdir -p /tmp/sahay-reference
git archive b87c0da396b4bbde54fbee5535aeac1116170605 tools/profile src/python/clifft/_folded | tar -x -C /tmp/sahay-reference
python tools/profile/study_sahay_gadget.py --source /tmp/wan-exact --reference /tmp/sahay-reference --package /tmp/clifft-main --work /tmp/sahay-gadget --output /tmp/sahay-gadget.json
python -m pytest tools/profile/test_wan_corpus.py -q
```

The work directory must not already exist. The driver pins itself to CPU 0 and
verifies the external manifest, reference planner/kernel and imported package.
It records source, model, extension and native-library hashes, the boundary and
rejection checks, exact probability comparisons, all trial times and record
hashes. Do not share CPU 0 with another benchmark.
