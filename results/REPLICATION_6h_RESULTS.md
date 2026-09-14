# §6h — replicating the differentiation adjustment in independent series

Plan fixed in PREREGISTRATION §6h **before any replication number existed**,
including the cohorts (by accession), the decision rule, and the acceptance in
advance of an "uninformative" outcome.

**The claim under test** is not H2 and not the composite. It is the §6g module
split: after adjusting for hepatocyte identity, does the REDUCTION rise in
tumor survive substantially better than the DRAIN fall?

---

## Summary

| | GSE76427 (§6g, original) | GSE164760 (cohort 1) | GSE14520 (cohort 2) |
|---|---|---|---|
| platform | Illumina HT-12 V4 | Affymetrix HG-U219 | Affymetrix HG-U133A |
| etiology | mixed | NASH | HBV |
| design | 115 T / 52 adjacent, **52 pairs** | 53 T / 29 adjacent, unpaired | 225 T / 220 adjacent, **213 pairs** |
| REDUCTION, unadjusted | +0.728 | +0.294 | **+0.836** |
| DRAIN, unadjusted | −0.712 | +0.048 ✗ | **−0.836** |
| REDUCTION, D1-adjusted | +0.517 (**71%**) | +0.248 (**84%**) | **+0.717 (84%)** |
| DRAIN, D1-adjusted | −0.295 (41%) | +0.209 (n/a) | **−0.197 (25%)** |
| REDUCTION, D2-adjusted | +0.505 (69%) | +0.242 (82%) | +0.692 (81%) |
| DRAIN, D2-adjusted | −0.145 (20%, *n.s.*) | +0.251 (n/a) | −0.066 (8%, *n.s.*) |
| **§6h verdict** | — | **uninformative** | **REPLICATED** |

D2 is the deliberately over-adjusted sensitivity model registered in §6g: it
conditions on the hepatic P450s, the very enzymes POR donates electrons to, so
it partly conditions on the drain itself. Under it the two arms separate
completely — REDUCTION holds ~80% of its shift in both paired cohorts while
DRAIN loses significance in both. That is the predicted behavior, not a
surprise, and it is reported as a bound rather than as the primary.

---

## Cohort 1 — GSE164760: uninformative by the registered rule, but not silent

The rule required the unadjusted shifts to be present *in the registered
directions* before adjustment means anything. They are not: **DRAIN is +0.048**
— essentially nil, and the wrong sign. So §6h returns **UNINFORMATIVE**, which
counts as neither support nor refutation. That is the call, and it stands.

Two things are nonetheless worth stating, because they are informative even
though the verdict is not:

1. **REDUCTION does replicate here.** Unadjusted +0.294 → D1-adjusted **+0.248
   (95% CI +0.008 to +0.488), P = 0.043, 84% retained.** Same direction as
   §6g, and a *higher* retention.
2. **The absent DRAIN shift has an obvious candidate explanation.** GSE164760's
   comparator is *non-tumoral NASH liver adjacent to HCC* — already diseased,
   already partly dedifferentiated. If the DRAIN fall in GSE76427 is largely
   loss of hepatocyte identity (which §6g found: only 41% survived), then
   against an already-dedifferentiated comparator it should shrink toward zero.
   It does. This is consistent with §6g rather than in tension with it — but it
   is a *post-hoc* reading of an uninformative cohort and is labeled as such.

The confound itself is present but much weaker here: D1 lower in tumor,
P = 0.033, against P = 2.8 × 10⁻¹⁴ in GSE76427 — exactly what an
already-diseased comparator predicts.

⚠ This cohort is **not outcome-blind**: its H2 was tested and was null
(P = 0.44), and its components were computed in PHASE C. §6h said so in advance.

---

## Cohort 2 — GSE14520: the module split replicates

445 samples (225 tumor, 220 adjacent), **213 patient-matched pairs**, HBV
etiology, a third platform. Pairing was recovered from GEO's own `Individual`
field and independently cross-checked against the LCS specimen code (213 vs
212 pairs — agreement).

**Unadjusted, this is the strongest replication in the study:**

| | estimate | P |
|---|---|---|
| H2 composite, Hodges–Lehmann | **+1.345** | 9.5 × 10⁻⁵⁶ |
| H2 composite, Cliff's δ | **+0.862** | — |
| paired ΔRSI, 213 pairs | **+1.396** | 8.5 × 10⁻⁵² |
| paired ΔREDUCTION | **+0.857** | 6.0 × 10⁻³⁹ |
| paired ΔDRAIN | **−0.800** | 3.9 × 10⁻²⁵ |

Compare GSE76427: HL +1.227, δ +0.759, paired ΔRSI +1.058. **The effect
reproduces at nearly identical magnitude in a different etiology, on a different
platform, with four times the pairs.** H2 is now supported in two independent
cohorts and the direction of both modules is reproduced in three.

### The registered adjustment, now run

The D1 covariate is complete — **22 of 22 genes**, D2 **26 of 26** — and §6i
records how (including the acquisition error that had made HNF1A look absent:
GPL3921 annotates it under the retired symbol `TCF1`, probesets `210515_at` and
`216930_at`, Entrez 6927). Adding HNF1A was **asserted numerically** not to move
a single index value: the composite H2 and both module deltas reproduce the
already-reported unadjusted numbers to four decimal places.

**Step 1, positive control.** D1 is far lower in tumor: **−0.364 vs +0.372,
P = 9.7 × 10⁻⁵⁹**. The dedifferentiation confound is not merely present here, it
is the strongest of the three cohorts — which makes this a hard test, not an
easy one.

**Step 2, primary — paired-difference regression over 213 pairs.**

| | intercept (RSI rise at zero change in hepatocyte identity) | 95% CI | P | retained |
|---|---|---|---|---|
| adjusted for D1 | **+0.792** | +0.659 to +0.924 | 1.2 × 10⁻²⁴ | **63%** |
| adjusted for D2 (over-adjusted) | +0.656 | +0.511 to +0.802 | 4.0 × 10⁻¹⁶ | 52% |

Unadjusted mean ΔRSI was +1.396 (median +1.261). **Step 3, supporting:** the
unpaired tumor coefficient over all 445 samples is +0.790 (+0.667 to +0.913),
P = 2.6 × 10⁻³¹.

**Step 4, module level — the actual §6h claim.**

| | unadjusted | D1-adjusted | 95% CI | P | retained |
|---|---|---|---|---|---|
| REDUCTION | +0.857 | **+0.717** | +0.580 to +0.853 | 2.2 × 10⁻²⁰ | **84%** |
| DRAIN | −0.800 | **−0.197** | −0.327 to −0.068 | 0.0032 | **25%** |

**Step 5, D2 sensitivity:** REDUCTION +0.692 (81%, P = 1.3 × 10⁻¹⁶); DRAIN
−0.066 (8%, **P = 0.36 — no longer distinguishable from zero**).

**§6h verdict: REPLICATED.** The unadjusted shifts are present in the registered
directions, REDUCTION's adjusted CI excludes zero, and REDUCTION's retention
(84%) exceeds DRAIN's (25%). **§6g verdict: NOT explained by dedifferentiation**
— the intercept is significantly positive and retains 63%, above the 0.5
threshold fixed in advance.

Two things are worth saying plainly about how this compares to §6g. First, the
split is **wider** here, not narrower: 84% vs 25%, against 71% vs 41% in
GSE76427. Second, GSE76427's composite retention was **50.3%** — a knife-edge on
the registered 0.5 threshold, and flagged as such at the time. GSE14520's 63% is
not knife-edge, and it was produced by a rule that had already been written down
when 50.3% was the only number in hand.

### Two platform limitations, stated not buried

- **AIFM2 is absent from HG-U133A**, so REDUCTION is a 3-gene version in
  GSE14520. This matters in a specific direction — and the direction first
  recorded here was wrong (corrected in PREREGISTRATION §6k, 2026-08-27). AIFM2
  is the *strongest* member of REDUCTION where it is measured (GSE76427
  Δz = +1.292, P = 2.31 × 10⁻¹⁶), so the platform drops the module's most
  informative gene from this cohort. The +0.836 is not strictly comparable to
  the 4-gene versions elsewhere, and the non-comparability runs against the
  hypothesis, not for it.
- COQ5 absent (SUPPLY 14/15) — immaterial.

---

## Where the module-split claim now stands

- **Direction replicates in three cohorts.** REDUCTION up and DRAIN down in
  tumor, in GSE76427 and GSE14520 decisively; REDUCTION up in GSE164760 too.
- **The adjustment now replicates in two independent cohorts** — different
  patients, different platforms, different etiology (mixed vs HBV) — under a
  rule fixed before the second cohort was acquired. REDUCTION retains 71% and
  84%; DRAIN retains 41% and 25%. Under the over-adjusted D2 model DRAIN loses
  significance in both (P = 0.22, P = 0.36) while REDUCTION does not.
- **GSE164760 remains uninformative by the registered rule** and is counted as
  neither support nor refutation. Its REDUCTION retention (84%, P = 0.043) is
  reported but not credited.
- **PHASE E §3a can now be upgraded** from a single-cohort observation to a
  two-cohort replicated finding — with the wording changed from "about half
  survives" to what the two cohorts actually show: **the composite retains 63%
  and 50% of its rise after adjustment for hepatocyte identity, and the
  retention is carried by the reduction arm.**

### What this does not license

- It does not make the decomposition causal. These are cross-sectional
  transcriptomes; "adjusting for hepatocyte identity" is a covariate model, not
  an experiment.
- **AIFM2's absence from HG-U133A remains a limitation in the cohort that
  gives the sharpest result, but it does not favor the hypothesis** (corrected
  in PREREGISTRATION §6k). The 84% retention is a 3-gene REDUCTION achieved
  without the module's strongest measured member. This belongs in the
  limitations, not in a footnote.
- GSE164760 is not outcome-blind and GSE14520's pairing was recovered by us from
  GEO metadata, not supplied as an analysis-ready design.

## Next

1. Rewrite PHASE E §3a as a two-cohort finding (above), and carry the D2
   dissociation into the manuscript — it is the cleanest single piece of
   evidence that the drain fall is largely identity loss while the reduction
   rise is not.
2. Remaining open items are unchanged and unrelated: H5 needs transcript-level
   quantification, H3 needs a larger adjacent-tissue series, and H6 needs a
   genome-wide re-prep of GSE164760.
