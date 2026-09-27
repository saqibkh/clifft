# Larger-distance cultivation scaling screen

2026-09-26. Current main is `63fc77ab7753220f396a87c274e9d2252a4b1b71`,
freshly fetched for this study. Production sources and the isolated installed
extension are the same as in the [bottleneck study](wan_bottleneck.md).

Follow-up: the [complete Sahay d5 gadget experiment](sahay_gadget.md) implements
and validates the missing physical binding, connects it to the original prefix,
and measures full noisy attempts on the supplied external model.

**Prioritize a complete check/projection contraction with offline selection of
syndrome-sampling order.** Larger distances strengthen the case for replacing
the dense representation. The current experimental contraction also needs
improvement: its syndrome-prefix order unnecessarily enlarges intermediate
tables. Changing that order gives a measured 3.8x arithmetic-kernel gain at d7
and makes the d9 geometry fit the existing numeric planner limits. This is not
yet an end-to-end sampler for the newly supplied circuits.

## What the papers establish

The papers describe several different kinds of speedup:

- Takada, Bartlett and Williamson give an exact, fault-conditioned reduction to
  Clifford stabilizer simulation for admissible cultivation protocols. They
  simulate fault distance f=7, with polynomial per-shot scaling; rare logical
  errors still require substantial sampling. This is evidence that exponential
  dense-state growth can be avoided for a suitable protocol, not a claim about
  arbitrary logical inputs or every cultivation circuit.
  [Paper](https://arxiv.org/html/2609.16929v1).
- Hartweg and Pineiro Orioli report residual MPS bond dimension five for Reg5,
  compared with active width 22 in Clifft. Their 90.5 ms versus 1.2 s comparison
  uses Clifft 0.7.0 and p=0 with noise instructions present. Separately, a
  prescreen rejects 99.7% of weight-five patterns, reducing calls to exact
  simulation by about 1000x. These are different benefits; neither number is
  a matched benchmark against our newly scheduled noisy baseline. Their
  reported circuit benchmarks cover d3/d5.
  [Paper](https://arxiv.org/html/2609.19116v1).
- Wan and Zapirain evaluate exact noise-averaged quantities for the supplied
  d3/d5 models, including series coefficients and selected all-order
  evaluations. Their computational task differs from generating complete
  noisy measurement records. The supplied appendix does not provide an
  independently validated d9 physical circuit.
  [Updated paper](https://arxiv.org/html/2609.18922v2).

Thus there is direct motivation to study d7 and beyond, but this study does not
reproduce a published d9 speedup. Distance labels also need care: a geometric
surface-code distance alone does not establish the fault distance of a whole
protocol with injection, growth, checks, escape and decoding.

## Inputs and comparison contract

The author repository remains pinned at
`b99452de7e300324e7bbf5cd9347241c3e289bae`. The Sahay d5 factor topology is
extracted from its two native folded blocks and terminal data checks. Both
blocks have the same interaction supports; the data checks commute. This is
the existing native-T reconstruction, with the qualifications in the
[corpus report](wan_corpus.md), not a matched noisy T export.

The experimental reference planner, numeric kernel and reconstructed f3/f5/f7
circuits come from `b87c0da396b4bbde54fbee5535aeac1116170605`. Their bytes are
checked against that commit. The regular geometry has 13/41/85/145/221 data
qubits and 6/20/42/72/110 X generators at d3/5/7/9/11. The d9/d11 cases extend
this code geometry and folded interaction pattern only. They are not complete
validated physical protocols. In particular, the older physical sampler uses
128-bit masks, which cannot hold the 145-data-qubit d9 model. The numeric
contraction has a separate representation and can be exercised independently.

The earlier small postselected overlap is insufficient to cost full sampling.
For r checks, the current full-syndrome scheme uses r+1 bra/ket contraction plans
for successive syndrome-prefix marginals, in addition to an amplitude plan.
This study counts those plans. Full physical execution additionally binds
faults, sums coherent-term/logical-sector contributions, samples outcomes,
maintains records and hands off the state. Those costs remain unmeasured here.

## Latest-main dense baseline and fusion alternative

Compile-only inspection uses current main's default HIR passes, then either
no scheduling, `ActiveWidthSchedulePass(search_budget=16)`, or unlimited
search. It calls `lower()` directly to avoid allocating an executor or running
the public compile API's default reference shot.

| Reconstructed circuit | Peak active width in all three configurations | Coefficient array alone |
| --- | ---: | ---: |
| f3 | 8 | 4 KiB |
| f5 | 22 | 64 MiB |
| f7 | 44 | 256 TiB |

The scheduler substantially improved measured d5 time in the prior study, but
it does not remove this d7 capacity barrier. We did not allocate the d7 state.

Broader fusion is still a plausible d5 optimization. The four expensive pure-X
runs in the actual Sahay profile account for 54.7% of full execution time:

| Run length m | X-mask rank r | Dense fused matrix entries | Diagonal phase entries | Separate arithmetic | Walsh/phase arithmetic |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 3 | 64 | 8 | 24 | 18 |
| 11 | 5 | 1,024 | 32 | 66 | 26 |
| 3 | 3 | 64 | 8 | 18 | 18 |
| 8 | 3 | 64 | 8 | 48 | 18 |

Arithmetic is a rough real-operation count per complex coefficient: 6m for
ordinary rotations, versus 4r+6 for two Walsh transforms and a complex phase.
The latter can diagonalize commuting X translations on each orbit. Counts
exclude addressing, phase binding, normalization setup and memory traffic;
they are hypotheses, not measured speedups. A dense fused matrix would cost
more than this structured transform. Noise-dependent signs must still be
bound correctly. Neither variant reduces the exponential coefficient array.
Even removing all rotation cost would cap the measured d5 overall gain at
about 3.8x if everything else stayed fixed. Fusion should therefore be a
secondary experiment for the user's larger-distance goal.

## An effective offline choice: syndrome-sampling order

The reference plans eliminate the variable with smallest neighboring scope,
breaking ties by variable index. Changing the internal order in which the
commuting terminal syndrome bits are sampled changes both the prefix networks
and those tie breaks. We compared row, column, diagonal and amplitude-plan
elimination orders. These are four inexpensive candidates, not an optimal
contraction-order search.

This is an internal probability factorization. Integration must map the
sampled bits back to their original physical record positions. It cannot
reorder arbitrary circuit measurements, discard rejected outcomes, or assume
commutation across intervening operations.

| Regular geometry | Order | Largest intermediate entries | Total marginal work | Numeric payload MiB |
| --- | --- | ---: | ---: | ---: |
| d5 | Current row | 256 | 23,792 | 0.63 |
| d5 | Amplitude elimination | 64 | 16,224 | 0.45 |
| d7 | Current row | 65,536 | 1,596,316 | 10.43 |
| d7 | Diagonal | 2,048 | 355,164 | 6.52 |
| d9 | Current row | 4,194,304 | 102,789,272 | 185.89 |
| d9 | Amplitude elimination | 8,192 | 4,706,712 | 28.66 |
| d11 | Current row | 1,073,741,824 | 17,659,477,356 | 38,133.99 |
| d11 | Diagonal | 16,777,216 | 224,693,156 | 630.96 |
| d11 | Amplitude elimination | 1,048,576 | 399,321,740 | 165.59 |

Work is the sum of joint-assignment counts over elimination steps and all
prefix marginals, not FLOPs or a runtime estimate. Payload includes lookup
arrays for all plans and reusable complex128 scratch/product arrays; it
excludes container overhead, physical bindings, local tables and surrounding
Clifft state. It includes the amplitude plan. The scope-only counter computes
these sizes without allocating the exponential tables.

At d7 the diagonal order reduces peak entries 32x and counted work 4.5x. At d9
the amplitude order reduces them 512x and 21.8x. All reordered d9 plans fit the
reference planner's two-million-entry storage and uncompressed-gather limits;
37 current-order plans fail those limits. We did not bypass the limits to time
the current d9 order.

At d11 both improved families still fail those limits. The diagonal order has
less counted work, while amplitude order has much less storage. This is not
evidence of polynomial scaling. Merely lifting limits would not address the
growth. The fault-conditioned Clifford approach deserves a separate
large-distance comparison even if a smaller-distance proxy was slower.

For the actual supplied Sahay d5 topology, amplitude order reduces the largest
table from 512 to 64 entries, marginal work from 23,120 to 16,224, and the total
numeric payload from 0.62 to 0.46 MiB. A full physical-instrument binding for
these blocks remains to be implemented and validated.

## Native arithmetic measurements

The pinned existing C++ contraction kernel evaluates directly emitted compact
plans. The new exporter creates separable low/high gather maps without first
allocating their full expanded counterparts. Its bytes are compared with the
reference exporter, including large compressed gathers.

| Geometry and order | Median ms per marginal sweep | Numeric payload MiB |
| --- | ---: | ---: |
| d7 current | 8.50 | 10.35 |
| d7 diagonal | 2.23 | 6.45 |
| d9 amplitude elimination | 27.15 | 28.46 |
| Sahay d5 current | 0.169 | 0.61 |
| Sahay d5 amplitude elimination | 0.131 | 0.44 |

**One sweep means one contraction query per syndrome prefix, not one physical
circuit attempt.** It omits the amplitude plan, physical-fault binding,
coherent-term pairs, logical sectors, outcome sampling and state/record handoff.
Its storage is correspondingly a little lower than the preceding table. In
particular, dividing the current-main 93.7 ms full-attempt baseline by 0.131 ms
would not be a valid speedup claim. These results establish a 3.8x d7 and 1.3x
Sahay d5 arithmetic-sweep improvement for the selected ordering experiment.

Five batches per model use one worker pinned to CPU 0 on the same EPYC 9554P
VM, GCC 13.3, `-O3 -march=native -fno-fast-math`, complex128, preallocated
workspace and changing local phases. Models run sequentially, not in an
interleaved comparison; use the raw trials to assess timing variation. Plan
construction and parsing are outside timing. Preparation took about 1.0-1.6 s
for d7 and 4.7 s for d9, including the reference-exporter validation. No
special-purpose new arithmetic kernel was required.

## Validation and next decision

- The scope-only counter matches 54 reference plans across regular
  d3/d5, selected d7 plans and the Sahay d5 topology, including exact compact
  address sizes.
- Reordering matches independent exhaustive Fourier/Born probabilities for
  four d3 phase instances and three orders: all 1,524 prefix probabilities,
  maximum absolute error `4.49e-16`.
- Direct compact serialization matches all seven d3 prefix plans and the
  largest-intermediate plan in each of the five native benchmark variants.
  The latter exercise joint tables of 64 through 65,536 entries, including
  the large low/high gather path.
- The native kernel's allocated lookup/workspace payload matches the offline
  prediction in every trial. This validates numeric-plan sizing and execution,
  not the unimplemented d9 physical protocol or new Sahay instrument binding.

The next experiment should add offline order selection to the experimental
full-syndrome contraction, then validate complete physical records and state
handoff on the reconstructed d7 control and supplied Sahay d5 blocks. The
latter's signed-code entry boundary is already established by the earlier
study; the whole two-check/projection instrument remains the missing proof.
Preserve every physical fault location, relative phase, accepted/rejected
branch and terminal observable. Compare full and postselected attempts
separately against current main with its opt-in scheduler.

Only after that bridge should d9 be presented as a full-circuit benchmark:
obtain or validate the complete protocol, extend physical masks, and check the
binding at that distance. Meanwhile the geometry screen can reject poor
orders cheaply through d11. Any eventual production implementation must keep
planning offline and hot execution preallocated, consistent with Clifft's
architecture. No production recognizer or executor was changed in this study.

## Reproduction and evidence

Use the isolated current-main package and research Python environment from
the preceding report; it needs NumPy and the package's dependencies. Extract
the pinned experimental reference without switching the working branch:

```bash
mkdir -p /tmp/wan-scaling-reference
git archive b87c0da396b4bbde54fbee5535aeac1116170605 tools/profile tests/fixtures/folded src/python/clifft/_folded | tar -x -C /tmp/wan-scaling-reference
python tools/profile/study_wan_scaling.py --reference /tmp/wan-scaling-reference --source /tmp/wan-exact --package /tmp/clifft-main --output /tmp/wan-scaling.json
python tools/profile/benchmark_wan_scaling.py --reference /tmp/wan-scaling-reference --source /tmp/wan-exact --work /tmp/wan-scaling-kernels --output /tmp/wan-scaling-kernels.json
```

The benchmark work directory must not already exist. Do not share CPU 0 with
another benchmark. The scope driver records extension and input hashes,
compile-only results, probability checks, fusion assumptions and per-distance
evidence. The numeric benchmark records plan/kernel hashes and all five trials.

Evidence: [provenance and Sahay plans](wan_scaling.json), regular
[d3](wan_scaling_d3.json), [d5](wan_scaling_d5.json), [d7](wan_scaling_d7.json),
[d9](wan_scaling_d9.json), [d11](wan_scaling_d11.json), and
[native measurements](wan_scaling_kernels.json). Elimination orders are stored
as hashes to keep artifacts small; the deterministic planner reconstructs them.
