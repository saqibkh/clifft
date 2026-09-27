# Wan-Zapirain cultivation corpus study

Later audit: [structural applicability across all supplied families](cultivation_applicability.md)
separates the pattern recognizer's fallback from local algebraic coverage.

2026-09-26. The released code unblocks independent circuit validation and
provides a concrete expensive folded workload. It does not yet demonstrate
generalization of the existing folded specialization: all eight imported
cases fall back. Continue a bounded Sahay d=5 certificate investigation;
defer generalizing or merging the current experimental specialization.

Later result: the [complete Sahay d5 gadget experiment](sahay_gadget.md) validates
a model-level replacement for both large checks and terminal projection. The
original circuit-level experimental recognizer described below remains unchanged.

Follow-up: the [latest-main bottleneck study](wan_bottleneck.md) compares the
shipped scheduler budget with unlimited search and attributes the remaining
cost to the large checking region. Its decision replaces the immediate next
experiment below with an offline comparison against broader commuting fusion.

## Pinned inputs and transformations

- Upstream: [exact_ler_for_msc](https://github.com/kh428/exact_ler_for_msc/tree/b99452de7e300324e7bbf5cd9347241c3e289bae),
  revision `b99452de7e300324e7bbf5cd9347241c3e289bae`.
- Ordinary Clifft: `a8350dfe027bd08a58525c85da7692d3c5b52d59`.
- Folded experiment: `b87c0da396b4bbde54fbee5535aeac1116170605`.
  All four installed `_folded/*.py` files match that revision byte for byte;
  its extension hash matches the previously recorded integration benchmark.
- [Exporter](../wan_corpus.py) checks both the committed upstream manifest and
  all 402 file hashes. The generated manifest records model/circuit hashes,
  source-to-output record mapping, every physical noise location, expected
  output parities, and the published reference coefficients.

The corpus contains RP2 d3/d5, Sahay Fold d3/d5, Chan flagged d3, Chan flagged
three-round/four-round d5, and the original SOFT d5 circuit. Input files remain
in an external checkout; this change does not redistribute the source circuits.
The package is Apache-2.0, but its upstream circuit notices still apply. Source
attributions and limits are in the upstream
[model descriptions](https://github.com/kh428/exact_ler_for_msc/blob/b99452de7e300324e7bbf5cd9347241c3e289bae/docs/OTHER_CULTIVATION.md).

The adapter reads the already converted actual-T gate/noise models. It does
not run the supplied S-proxy `.stim` files as T circuits. In particular:

- Native `CSX`, `CS_DAG_X`, and injection gates are decomposed exactly. Native
  gate noise remains after the complete gate; no decomposition noise is added.
  The conditional phase on a cat control is retained.
- Every noise site uses the requested common p. Idle noise and single-, two-,
  and three-qubit channel supports remain intact.
- Chan's `READOUT_FLIP` becomes noisy measurement reporting. It does not apply
  a Pauli to the measured quantum state. The Fold d5 model's explicit reported
  bit register is retained as supplied.
- Recorded feedback, inverted measurements, and detector reference parities
  are preserved. `POST_Z` becomes measurement plus a zero-result detector.
- The RP2/Fold references' final `(I-(X_L+Y_L)/sqrt(2))/2` error effect is
  implemented as noiseless logical T-dagger followed by logical X measurement.
  Chan's supplied logical-observable record parity is retained.

These models exclude escape and decoding. Sahay Fold d5 is a declared native-T
reconstruction with idle noise and a different schedule from our Takada-derived
f5 fixture. It is not an author-supplied noisy T export. There is still no new
independent d7 artifact in this release.

## Correctness evidence

Seven focused adapter tests pass, including entangled-input comparisons against
Qiskit Aer for the three custom gates, a Stim measurement-reporting comparison,
feedback/parity checks, and exact logical-readout/reset-marginal checks.
Every imported circuit has trivial detectors and no logical errors in the
bounded p=0 runs, under both default and scheduled compilation.

The original SOFT input has SHA256
`859d925c74712efaae8bab1d944f9e5be4e042e226eeafcc6d7bef9d6fbfbae5`,
identical to `tests/fixtures/cultivation_d5.stim`. At p=0.001, scheduled
postselected sampling gives acceptance 0.145239 over 51,591 attempts, versus
the exact 0.143955584096: a difference of 0.83 binomial standard errors. The
default-pipeline estimate is 1.70 standard errors from the same reference.
These sample sizes cannot resolve its approximately 3.33e-9 conditional LER.

The strongest logical-error check sums the complete 32-record space satisfying
Chan's detector constraints, including both logical outcomes. It uses full
injection, growth, and checking source models with the specified four faults
and all other faults absent. Current Clifft's exact record-probability API gives:

| Circuit and condition | Acceptance A | Accepted error B |
| --- | ---: | ---: |
| Three-round, fault-free | 1 | 0 |
| Three-round, published four faults | 0.2500000000000001 | 0.2500000000000001 |
| Four-round, fault-free | 1 | 0 |
| Four-round, same translated faults | 0 | 0 |

For this oracle only, resets allocate fresh initially-zero wires and discarded
wires are traced by the probability query. This avoids the API's hidden-reset
restriction without conditioning on arbitrary reset outcomes. It preserves
the visible quantum instrument. Timed sampling uses the original physical
resets and wire count. A test covers resetting half of an entangled pair.

Single-fault stratified sampling checks the first acceptance coefficient for
all seven additional models. When a reference sums identity-response noise
sites out, the physical conditional acceptance is
`(a_1 + N_physical - N_reference) / N_physical`. All seven estimates are within
2.22 standard errors and produce no accepted logical errors. Sample sizes are
65,536 for the small cases and Chan d5 variants, 2,048 for RP2 d5, and 128 for
Fold d5. The last is a coarse check, not a precision estimate.

Both upstream certificate verifiers pass. A fresh Fold d3 polynomial contraction
through degree five, for one prime and one sqrt(2) embedding, also reproduces
all eight slices of each A/X/Y observable and the assembled B coefficients
(126.3 seconds including kernel preparation). This reruns the supplied physics
implementation; it is not an independent rational reconstruction or a full d5
campaign. The additional variants have no supplied finite-series remainder
bounds; finite-series LER values are not treated as all-orders ground truth.

## Bounded CPU measurements

GCC 13.3 Release, native AVX-512, AMD EPYC 9554P, Python 3.13.11, CPU 0,
one worker, scalar batches. Each entry is the median of three timed batches,
excluding compilation and warmup. Postselection keeps survivor records.
Scheduled means the default passes plus `ActiveWidthSchedulePass(search_budget=None)`.
All timings include noise draws and rejected attempts. Both pipelines use
identical circuit bytes, output parities, p=0.001, and thread budget.

| Circuit | Peak active width | Default full us/attempt | Scheduled full us/attempt | Scheduled postselected us/attempt |
| --- | ---: | ---: | ---: | ---: |
| RP2 d3 | 8 | 2.79 | 2.43 | 1.75 |
| RP2 d5 | 18 | 3,559 | 1,620 | 368 |
| Sahay Fold d3 | 8 | 5.80 | 2.83 | 2.42 |
| Sahay Fold d5 | 22 | 633,530 | 84,443 | 29,551 |
| Chan flagged d3 | 4 | 0.99 | 0.97 | 0.82 |
| Chan flagged three-round d5 | 10 | 17.51 | 15.71 | 6.72 |
| Chan flagged four-round d5 | 10 | 18.46 | 18.21 | 7.64 |
| SOFT d5 | 10 | 13.02 | 9.23 | 3.46 |

Widths do not change with scheduling. Its reduction of dense work accounts for
the observed approximately 7.5x Fold d5 full-sampling gain; this is an existing
compiler optimization, not a representation/specialization improvement.
Scheduled Fold d5 compilation is 520 ms versus 34 ms by default. The coefficient
array at width 22 alone occupies 64 MiB; whole-process RSS is recorded separately
and is a cumulative high-water mark across cases, not per-worker storage.

The adaptive batch sizes keep this assessment bounded. Fold d5 has only 24-27
attempts per mode/pipeline; its postselection timing is particularly sensitive
to which trajectories survive. Small differences are not statistically
established. Raw JSON retains every timing, batch size, accepted/error count,
compile cost, scheduler work estimate, and package/circuit hash. No rare-LER
throughput or best-in-class claim follows from these runs.

## Eligibility and next decision

All eight cases are rejected by the experimental recognizer, and actual
`compile(..., specialize_folded=True)` confirms ordinary fallback. The original
f3 fixture still selects specialization with peak active width one as a
positive control.

The immediate rejection is the recognizer's global exclusion of record and
Pauli-product targets: these inputs contain terminal MPP/readout operations,
and several contain feedback. Removing that guard alone would be unsound and
would not establish schedule compatibility. The Sahay model also has native
controlled-gate noise/idle layers, cat phases, measurement-reporting behavior,
and terminal code projection that differ from the exact template accepted by
the current terminal-region certificate. The colour-code and RP2 cases require
different block families altogether.

There is nevertheless a concrete favorable diagnostic for Sahay Fold d5.
Immediately before native block 2, after its cat preparation, ordinary Clifft
ends with one active coordinate. All 40 terminal data-code stabilizer probes
compile to identity on that coordinate with tracked signs. In 1,024 noisy
prefix shots their absolute expectation values are exactly one. The next check
contains 5 CSX, 4 CS_DAG_X, and 16 CCZ gates. The source record and probe details
are preserved in `wan_boundary.json`.

This establishes the signed **data-code** boundary, not classical ancillas or a
complete two-check instrument. A useful next experiment is an offline,
fault-conditioned contraction for this precise independently selected model,
including its actual cat preparation, both checks, all record outcomes, idle
noise, and final projection. Validate accepted and rejected record
probabilities and relative phases against ordinary Clifft before attempting a
production recognizer extension. Compare end-to-end cost against the scheduled
baseline above. If supporting this requires broad protocol-specific machinery,
retain it as a reference case and defer that extension.

No executor, kernel, public API, or architectural invariant changes are made
by this study.

## Reproduction and artifacts

Use an isolated Python 3.12+ environment with NumPy, Stim, pytest, Qiskit and
Qiskit Aer. Build the specified Clifft revisions into separate non-editable
package directories. The driver checks the imported package path and records
the extension hash. It must not load a different editable installation.

```bash
git clone https://github.com/kh428/exact_ler_for_msc.git /tmp/wan-exact
git -C /tmp/wan-exact checkout b99452de7e300324e7bbf5cd9347241c3e289bae
python tools/profile/wan_corpus.py --source /tmp/wan-exact --output /tmp/wan-corpus --p 0.001
python tools/profile/study_wan_corpus.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/current-clifft --output /tmp/baseline.json
python tools/profile/study_wan_corpus.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/current-clifft --schedule --output /tmp/scheduled.json
python tools/profile/study_wan_corpus.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/current-clifft --strata-only --output /tmp/strata.json
python tools/profile/study_wan_corpus.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/current-clifft --boundary-only --output /tmp/boundary.json
git show b87c0da3:tests/fixtures/folded/f3.stim > /tmp/folded-control.stim
python tools/profile/study_wan_corpus.py --source /tmp/wan-exact --corpus /tmp/wan-corpus --package /tmp/folded-clifft --eligibility-only --control-circuit /tmp/folded-control.stim --output /tmp/eligibility.json
python -m pytest tools/profile/test_wan_corpus.py -q
```

From the upstream checkout, run `python -B verify_results.py`,
`python -B verify_other.py`, and
`python -B recompute_other.py polynomial fold-d3 --degree 5 --output /tmp/fresh-fold3`.
The last command also requires a C++17 compiler and NumPy.

Recorded evidence: [default baseline](wan_baseline.json),
[scheduled baseline](wan_scheduled.json), [one-fault checks](wan_strata.json),
[eligibility](wan_eligibility.json), [boundary probes](wan_boundary.json),
[original certificates](wan_certificates_original.json),
[additional certificates](wan_certificates_alternatives.json), and
[fresh modular contraction](wan_fold_d3_recomputed.json).
