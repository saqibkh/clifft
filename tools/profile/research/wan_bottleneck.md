# Sahay d5 bottleneck study on current main

2026-09-26. Main revision `63fc77ab7753220f396a87c274e9d2252a4b1b71`, freshly
fetched before the study. Production sources match that revision. This study
adds research tools and evidence only.

**The dominant remaining cost is the two large folded checks and their terminal
code projection.** The latest opt-in scheduler substantially improves the
baseline, but does not remove this bottleneck. Continue with an offline cost
comparison of a complete check/projection contraction and broader commuting
rotation fusion before implementing either production extension.

This follows the [corpus validation study](wan_corpus.md). Its pinned models,
physical noise, reference parities, terminal readout, and limits still apply.
In particular, the supplied Sahay d5 case is a native-T reconstruction with
2,116 physical noise sites, not a matched author-supplied noisy T export. It
excludes escape and decoding.

## Latest-main baseline and opt-in pass

The recently landed pass is `ActiveWidthSchedulePass` from main's `e401286a`
(PR 507). It runs after the default HIR passes. Compare ordinary compilation,
the shipped default search budget of 16, and `search_budget=None`:

| Sahay d5 configuration | Compile ms | Full ms/attempt | Postselected ms/attempt |
| --- | ---: | ---: | ---: |
| Default pipeline | 35 | 686.3 | 240.4 |
| Add scheduler with default budget | 113 | 116.9 | 37.6 |
| Add scheduler with unlimited search | 519 | 93.7 | 29.8 |

The scheduler gives approximately 5.9x and 7.3x full-sampling gains respectively.
Unlimited search improves on the default scheduler by about 1.25x. It costs
roughly another 0.4 seconds to compile, recovered after about 18 full attempts
or 55 postselected attempts at the measured median rates. These crossover
estimates depend on workload and timing variation.

All three plans peak at width 22. The scheduler's estimated dense work falls
from 719,447,418 to 152,183,576 with its default budget, and to 128,363,144 with
unlimited search. Width alone therefore misses the main scheduling gain.
The coefficient array alone is 64 MiB; this is not total worker memory.

Five trials per configuration use rotating/reversed order, one CPU worker,
scalar batches, retained records, and p=0.001. Full and postselected programs
use identical input bytes and expected output parities. Each scheduled mode
has 640 attempts; each unscheduled mode has 160. Compilation and two warmup
attempts are outside timing. Full-sampling trial ranges are 653-697, 116-123,
and 91-97 ms respectively. Postselected ranges are 117-321, 30-42, and 24-33 ms;
the unscheduled postselection estimate remains particularly noisy.

The hardware is the same EPYC 9554P VM, pinned CPU 0, GCC 13.3 Release/native
AVX-512, Python 3.13.11. Results compare these CPU configurations, not every
available batch/thread setting or accelerator. The build was installed into
an isolated package directory, and the driver checks its imported path and
records the extension hash. The only upstream change since the prior study's
`a8350dfe` baseline concerns Sinter dependencies; it changes no sampling kernel.

Controls retain the same distinction between useful scheduling and a need for
a new representation:

| Circuit | Default full us | Default-budget scheduler full us | Unlimited scheduler full us |
| --- | ---: | ---: | ---: |
| RP2 d5 | 3,618 | 1,605 | 1,614 |
| Chan four-round d5 | 18.2 | 16.5 | 16.6 |

These use 5,120 and 163,840 attempts per configuration. Unlimited search has no
measured throughput advantage on these controls; Chan's two scheduled executable
plans are identical. The Sahay result does not justify using unlimited search
for every circuit.

Raw compile times, plan hashes, scheduler statistics, all trial times and counts
are in [wan_main_scheduling.json](wan_main_scheduling.json).

## Where the optimized time goes

The source-attributed profiler runs the unlimited-search plan. An isolated
copy of the executor translation unit adds steady-clock counters around each
existing action. Counter storage is preallocated, and all other production
objects and compiler flags come from the current-main build. The ordinary
library, executor source, public API, and hot-path architecture are unchanged.

Native plans retain source provenance solely for attribution. After removing
those debug annotations, both native executables' inspection output exactly
matches the Python public API's plan. Instrumented and uninstrumented executors
produce identical measurement-record hashes, acceptance counts and error counts
for 256 full attempts and 1,024 postselected attempts with the same seed. A
separate d3 control also matches the record hash over 256 attempts.

| Source region | Full execution time | Postselected execution time |
| --- | ---: | ---: |
| Prefix, including injection/growth and smaller checks | 0.02% | 0.04% |
| First large check and intervening ancillas | 42.8% | 47.5% |
| Second large check and readout | 46.2% | 42.5% |
| Terminal code projection | 11.0% | 9.9% |

These percentages use instrumented executor wall time, including the sampled
attempts that reject. Source regions are attributed through the actual
optimized plan, not by separately compiling and subtracting circuit prefixes.
The JSON records the exported source-line ranges. The first-check interval
also includes preparation for the next check; Clifford-only gates absorbed by
the compiler have no separate runtime charge.

Full sampling spends 73.5% in rotations, 18.3% in active measurements, and 8.1%
in promotion of new active coordinates. Postselection gives 72.8%, 18.5%, and
8.5%. About 86% of execution time is at widths 21 and 22. Less than 0.1% of
executor wall time is outside the instrumented actions. This residual includes
lazy noise handling, resets and dispatch/timer overhead; it is not an isolated
measurement of noise-generation time.

Independent `perf` cycle sampling on the uninstrumented executable corroborates
the result: approximately 71% rotations, 18% measurement probability/collapse,
and 8% promotions for full sampling. It collected 12,536 full and 16,530
postselected samples without lost samples. This profile includes compilation
and process startup, unlike the action timers.

The native postselected run keeps 167/1,024 attempts. In the instrumented run,
survivors average 92.9 ms and rejected attempts average 19.3 ms, for 31.3 ms per
attempt overall. The uninstrumented run averages 31.8 ms. Full instrumented and
uninstrumented runs average 89.1 and 95.9 ms. They ran sequentially, with `perf`
attached only to the uninstrumented run; these differences do not establish
timer overhead or a speedup from instrumentation. Use the separate interleaved
public-API benchmark for performance comparisons. No rare logical-error-rate
claim follows from these samples.

## Consequences for the next experiment

1. **Target the complete pair of large checks and terminal projection.**
   Optimizing the prefix or noise plumbing has negligible measured headroom.
   Replacing just the first check has an idealized full-sampling speedup ceiling
   of about 1.75x. Replacing both checks while leaving terminal projection at its
   current cost caps the gain near 9x. These are zero-cost-replacement bounds,
   not predicted speedups; an integrated replacement can also change the cost
   of the projection.
2. **Keep a cheaper compiler alternative in the comparison.** There are zero
   static or dynamic fused rotation groups in the optimized plan. Several
   expensive fixed-width runs contain only X Paulis, so their rotations commute
   even when noise changes their signs. The four largest runs account for about
   55% of full execution time. Their lengths are 4, 11, 3 and 8, with exact
   binary X-mask ranks 3, 5, 3 and 3. Current fusion caps this rank at two;
   dynamic fusion additionally limits independent signs and requires eight
   actions. Merely switching on existing fusion does not cover these runs.
   Offline orbit/table-size and arithmetic estimates should determine whether
   a broader commuting kernel deserves a small benchmark. This observation
   does not establish that a larger dense fused matrix would be faster.
3. **Cost the logical-block contraction over the same measured interval.**
   The earlier signed data-code boundary makes it a plausible way to avoid
   these dense excursions. Construct its elimination plan and bound the
   largest intermediate, total work, and per-shot noise/record binding. Retain
   all physical fault locations, accepted/rejected outcomes and relative phases.
   Compare against 93.7 ms full or 29.8 ms postselected attempts, including the
   different entry/exit and rejection behavior. A small logical boundary alone
   does not prove affordable intermediate contractions.

This completes the bottleneck decision step. It does not extend the recognizer,
implement a new contraction/fusion kernel, or establish performance of either
candidate. It narrows the next comparison to two concrete ways of reducing the
same dominant dense work.

## Validation and reproduction

The seven existing adapter tests and all 18 current active-width scheduler
tests pass against the latest-main build, including the independent Aer
unitary/noisy-scheduling references. Record-hash and plan-identity checks above
validate the profiling instrumentation against the unchanged executor.

Build current main with the options from the prior study and
`-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`. Keep the package installation and the
research Python environment isolated from editable-import hooks. Then run:

```bash
python tools/profile/measure_wan_scheduling.py --package /tmp/clifft-main --corpus /tmp/wan-corpus --output /tmp/wan-main-scheduling.json
python tools/profile/build_wan_action_profile.py --build /tmp/clifft-build --output /tmp/wan-native
python tools/profile/study_wan_bottleneck.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/clifft-main --native /tmp/wan-native --raw /tmp/wan-profile --output /tmp/wan-bottleneck.json --perf /path/to/perf
```

The native builder creates both timed and untimed executables and records
source, library and executable hashes. It refuses to patch an unexpected
executor layout. `--perf` is optional. Run the performance measurements
sequentially, with no other benchmark sharing CPU 0.

Evidence: [profile summaries and provenance](wan_bottleneck.json),
[full action counters](wan_bottleneck_full_actions.json), and
[postselected action counters](wan_bottleneck_postselected_actions.json).
