# Cultivation gadget applicability audit

2026-09-27. **The underlying small-branch mechanism is present in more families
than the current fast folded sampler handles.** The audit finds qualifying
local inverse-T parity instruments in SOFT, Chan and RP2, as well as the
already studied folded checks in Sahay. This supports a broader research
direction. It does not establish a broadly useful production optimization or
speedups on those additional families.

The earlier result that every external circuit fell back described one exact
folded-pattern recognizer. It understated structural applicability when used
as evidence about the mechanism itself.

## Scope and definitions

The audit covers **nine supplied circuits across four families**: the previous
eight-case Wan corpus plus the supplied SOFT d3 file. Larger reconstructions
and constructed controls are listed separately and do not count as independent
external families. Wan's frozen models are pinned to
`b99452de7e300324e7bbf5cd9347241c3e289bae`; the manifest and all distributed file
hashes are verified before use. Sahay d5 remains a supplied native-T
reconstruction, not a matched author-supplied noisy T export.

Main was checked remotely at `63fc77ab7753220f396a87c274e9d2252a4b1b71`. The
ordinary compiler and boundary probes use that build. The older experimental
folded pattern/certificate is pinned to
`b87c0da396b4bbde54fbee5535aeac1116170605`. No production source or optimizer
behavior changes in this audit.

There are three different questions:

1. **Recognition:** does an existing frontend select a replacement from the
   supplied representation?
2. **Algebra:** does a selected local region have a certified small coherent
   expansion, and can its input be represented compactly?
3. **Useful execution:** can we sample the complete noisy output contract,
   including continuation, more efficiently than current Clifft?

A success at the second level is not reported as a success at the third.
[The evidence](cultivation_applicability.json) records source locations,
certificates, actual rejection reasons, prefix probes and size bounds.

## Results by supplied family

| Supplied circuits | Local structure established | Current fast folded path | Remaining work |
| --- | --- | --- | --- |
| Sahay Fold d3 | Two cats, four paths; signed code entry verified | Passes the existing model certificate after an equivalent decoder rewrite | Bind/validate the replacement on this representation; no new timing |
| Sahay Fold d5 | Two cats, four paths; signed code entry verified | Existing validated research replacement | Ordinary circuit recognition and native integration |
| SOFT d3/d5 | One/two inverse-T parity regions, two terms per region; largest-code entry verified | Not handled by the folded path | Separate parity lowering, physical noise/records and terminal continuation |
| Chan d3, three-round d5, four-round d5 | One/two inverse-T parity regions; 2/13 flags preserve two terms per region; largest-code entry verified | Not handled by the folded path | Flag-aware parity lowering, reported-bit noise, full continuation and validation |
| RP2 d3/d5 | One/two inverse-T parity regions; 9/15 flags preserve two terms per region | Not handled by the folded path | Intermediate morphed-code boundary, parity lowering and full continuation |

The inverse-T scan finds **11 qualifying local regions in seven circuits**.
Repeated smaller-distance regions in different files are not independent
methods or extra families. The two-term bound is local and conditioned on
the measurement record and physical Pauli fault pattern. Composing many such
regions can still increase coherent rank unless a valid boundary compresses it.

The data dimensions make the RP2 limitation concrete. Its final CSS codes
have 9/25 data qubits with independent X/Z ranks 4/12. The selected T layers
instead act on 15/35 qubits. This audit does not infer a compact intermediate
code by substituting the final code. Chan and SOFT use phase supports matching
their 7/19-data-qubit final codes, with X/Z ranks 3/9.

## The existing folded contract

The Sahay model certificate requires:

- An independent, balanced one-logical-qubit CSS basis and explicit logical
  X/Z supports at the terminal boundary.
- A signed encoded logical state and classical ancillas at entry.
- Two cat preparation/readout pairs. With the cat branches fixed, the body
  must use supported monomial gates, computational ancilla records/resets,
  and fixed one- or two-parity factors.
- Records and discarded reset bits independent of unresolved coherent paths.
- Branch-dependent data flips differing only by X stabilizers and logical X.
- A supported terminal measurement/continuation contract and affordable
  contraction plans.

The narrow implementation additionally uses native-block annotations to find
candidates, assumes four Hadamards in its model region, and expects noiseless
terminal CSS products. These are implementation constraints, not universal
conditions for a useful coherent decomposition.

All nine supplied circuits still fail the older exact-pattern recognizer,
whose global target guard rejects Pauli-product/record targets before examining
the region. The model-level Sahay recognizer accepts Fold d5 as supplied.
Literal failures elsewhere include:

- Fold d3 expresses the final code through a noiseless decoder and single-qubit
  postselection, rather than explicit CSS product measurements.
- RP2 has redundant CSS check records and different intermediate check support.
  The raw check counts violate the recognizer's independent-basis assumption;
  rank reduction shows that its final code still encodes one logical qubit.
- Chan d3 includes its logical measurement among product measurements; blindly
  treating every product as a code generator is incorrect. Chan d5 uses another
  event schema. Their actual measured-and-reset parity circuits also differ
  from the cat-controlled monomial body expected by the folded certificate.
- SOFT has no native-block model annotations and uses the inverse-T parity
  construction and its own continuation.

Qubit permutation and exact controlled-gate decomposition both preserve Fold
d5 model eligibility. Removing only its native-block annotations causes
rejection while leaving its quantum circuit unchanged. Thus the existing
candidate finder is not a representation-independent structural optimizer.

The pinned reconstructed f3/f5/f7 fixtures pass their separate template and
region certificate. Constructed d3/d5/d7/d9 models pass the model certificate.
The earlier full-f7 growth validation and constructed-d9 sampling remain
useful scaling controls, not additional evidence of external-family coverage.

## Why the other parity circuits have two coherent terms

[The local audit](../parity_envelope_audit.py) checks actual operations:

1. The two T layers are distinct, signed inverses on the same qubits.
2. There is exactly one central X measurement followed by an X reset.
3. Extra flags are prepared on their first use, measured on their last use,
   use matching preparation/readout bases, and are outside the T support.
4. After those flag endpoints are separated, the two CNOT networks are exact
   inverses as Clifford tableaux, including flags and spectator qubits.
5. The pulled-back central X is all X on the T support and commutes with every
   prepared flag's basis. Therefore the coherent branch does not distinguish
   the flag readout. Flags can still report faults or depend on classical
   measurement outcomes.

Let the first phase layer be D and the first CNOT network A. For central
outcome m, define P = A^-1 X_root A and Q = A^-1 Z_root A. Before contracting
the prepared/measured flags, the ideal local Kraus operator is:

```text
K_m = (D^-1 Q^m D + (-1)^m D^-1 Q^m P D) / 2.
```

Each term is a tensor product of single-qubit Clifford operators: conjugating
a Pauli by a single T or T-dagger gives a Clifford operator. The two terms
must be added coherently. They are not two classical alternatives to sample.
The reset correction Q is necessary for the conditional output state.

For a fixed physical Pauli fault pattern, middle Clifford faults commute to
Pauli frames and signs; faults in the T layers can also introduce single-qubit
Clifford corrections. The local branch count remains bounded. This algebraic
observation does not implement the fault-to-frame maps or prove a complete
noisy sampler. No non-Pauli noise channel is covered by this argument.

The existing earlier unflagged parity reference already covered SOFT's form.
The useful new audit result is that explicitly checking flag endpoints and
the full Clifford inverse relation also covers the supplied Chan and RP2
checking regions. Four malformed controls are rejected: a missing inverse
CNOT, a wrong flag readout basis, an unclosed flag, and a noninverse T layer.

## Boundary and representation checks

At the largest-code parity region of each of the three Chan and two SOFT
inputs, current main compiles all terminal CSS stabilizers and a single-qubit
axis for every involved spectator, including entry flags, into tracked signs.
These probes act at width one. In each case, 128 noisy prefix samples have exactly unit
absolute expectation for those signed probes. The prefix peaks at width one
for d3 and four for d5. Logical X/Y/Z purity is checked as well, to detect
entanglement with unrepresented spectators; its maximum error is `8.89e-16`.
This is encouraging entry evidence; it does not validate the region's full
output instrument or its subsequent continuation.

Single-bit-factor contraction planning on those code bases has a peak of
eight entries and numeric plan/workspace payloads of 4,040 bytes at d3 and
28,868 bytes at d5. These are a **scope-only screen** for a prospective parity
implementation, excluding physical fault binding, actual continuation and
runtime overhead. They are not measured sampler memory or throughput.

For Fold d3, the audit pulls its 24-gate noiseless decoder back into code
measurements and physical logical axes. It retains all original measurement
records and noise sites. Complete physical current-main Clifft probabilities
agree on 32 selected records across four fault patterns, including rejected
records. The rewritten model passes the existing folded certificate and all
1,024 noisy entry checks. Fold d5's 1,024 entry checks were also rerun.

## Release decision and next experiment

The defensible release claim today remains the narrow, already measured
folded-cultivation case. A broad speedup claim is premature. Only external
Sahay d5 has a matched full-attempt performance comparison for this research
replacement: about 10x in the recorded scalar CPU configuration. The new
local certificates do not extend that speedup measurement to Chan, SOFT or RP2.
Current-main Chan d5 already runs in roughly 16 microseconds per complete
attempt; Python research binding cannot be treated as competitive with it.

However, the research pitch can now be more specific:

> Recognize cultivation checks that reduce to a small coherent sum of Clifford
> actions, and compile their physical faults, records and code transitions
> into fixed native sampling plans.

**Prioritize one complete flagged Chan four-round d5 instrument next**, using
this two-term certificate, before expanding the folded-specific production
path. This is an externally supplied distinct family with a verified compact
entry and a concrete flag/record challenge. Preserve its exact noise and
logical-readout conventions, compare both accepted and rejected complete
records, and include the known four-fault witness. Then measure native
end-to-end cost against current main's opt-in scheduler. A slower result is
useful evidence against a general performance release, even if the algebra
works. Keep RP2 conditional until its intermediate boundary is certified.

More author-provided large circuits are not the immediate blocker. The audit
identifies a concrete extension to test on existing inputs, and a clear gap
between structural applicability and useful released performance.

## Reproduction

Use the isolated current-main Clifft package with NumPy and Stim. Extract the
pinned experimental `src/python/clifft/_folded` and fixture files into the
reference directory. Run with a new work directory:

```sh
python tools/profile/audit_cultivation_applicability.py \
  --source /path/to/exact-ler-pinned-checkout \
  --reference /path/to/pinned-folded-reference \
  --work /tmp/cultivation-applicability \
  --output /tmp/cultivation-applicability.json
```

The report records source hashes, reference revisions, the main extension
hash, line-local candidate outcomes, exact Clifford identities, symbolic
entry probes, selected probability comparisons and negative controls. No new
performance measurements are part of this audit.
