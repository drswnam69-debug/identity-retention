# PHASE B results — discovery and replication

**Run:** 2026-08-25 · **Index lock hash:** `7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`
**Cohort:** 216 liver biopsies — control 10, NAFL 51, MASH F0–1 34, F2 53, F3 54, F4 14
(group counts read from GEO and verified against the series totals)
**Outcome variables remain unopened.** Survival and tumor status are PHASE E.

---

## 1. H1 is supported

RSI rises monotonically across the MASLD spectrum.

| Test | Statistic | P |
|---|---|---|
| **Jonckheere–Terpstra, increasing** (pre-specified primary) | z = 5.98 | **1.1 × 10⁻⁹** |
| Spearman ρ vs ordered stage | 0.397 (95% CI 0.270–0.509) | — |
| Stage term, adjusted linear model | β = +0.218 (SE 0.032) | 6.7 × 10⁻¹¹ |
| RSI vs NAS score | ρ = +0.418 | 1.6 × 10⁻¹⁰ |
| RSI vs fibrosis stage | ρ = +0.394 | 2.0 × 10⁻⁹ |

**The trend is not merely disease versus normal.** Controls are only 10 samples
and sit well below the rest, so the obvious objection is that the whole signal
is control-vs-disease. It is not:

| Subset | n | JT z | P | ρ |
|---|---|---|---|---|
| All six groups (pre-specified) | 216 | 5.98 | 1.1 × 10⁻⁹ | 0.397 |
| **Controls excluded (NAFL → F4)** | 206 | 4.82 | **7.2 × 10⁻⁷** | 0.331 |
| **MASH only, by fibrosis stage (F0–1 → F4)** | 155 | 3.05 | **0.0012** | 0.242 |
| Control vs any disease (Mann–Whitney) | 216 | −4.04 | 5.4 × 10⁻⁵ | — |

Both effects are present: a step at the onset of disease, and a graded rise
with fibrosis within it. The second is the one that matters for the hypothesis,
and it survives on its own.

## 2. Sensitivity analyses — the trend is robust, the *interpretation* is not

| Variant | All 216 | Controls excluded |
|---|---|---|
| Locked (pre-specified) | z = 5.98, ρ = +0.40 | z = 4.82, ρ = +0.33 |
| **S2 + cholesterol-arm genes** | z = 5.95, ρ = +0.39 | z = 4.76, ρ = +0.33 |
| S3 CYB5R3 alone as REDUCTION | z = 4.56, ρ = +0.30 | z = 3.15, ρ = +0.22 |
| S4 POR removed from DRAIN | z = 5.55, ρ = +0.37 | z = 4.00, ρ = +0.28 |
| post hoc: NQO1 removed | z = 5.09, ρ = +0.34 | z = 3.73, ρ = +0.26 |

**S2 must be reported as specified.** The pre-specification stated that if
adding cholesterol-arm genes (SC5D, DHCR7, DHCR24, MSMO1) left the result
unchanged, the claim that the index is CoQ-specific rather than
cholesterol-driven is **weakened**. It did (z 5.98 → 5.95). That is what will
be reported.

## 3. What actually moves — and it is not what the hypothesis emphasized

Component trends, **controls excluded** (the honest comparison):

| Module | z | P | ρ |
|---|---|---|---|
| SUPPLY | 1.35 | 0.089 | +0.097 |
| REDUCTION | 1.80 | 0.036 | +0.127 |
| **DRAIN** | **−4.25** | — | **−0.292** |

Within the disease spectrum the index is carried mainly by the **drain falling**,
not by mevalonate supply rising. The framing that motivated the study — supply
expands — is at best weakly supported in human tissue.

### Per-gene, controls excluded (BH-FDR across 33 genes)

**Rising, q < 0.05:** NQO1 +0.412 (q = 2.4 × 10⁻⁸) · IDI1 +0.201 · PDSS1 +0.199 ·
ACSL4 +0.188
**Falling, q < 0.05:** LPCAT3 −0.287 (q = 4.7 × 10⁻⁴) · MTARC2 −0.243 · PMVK −0.197
**Directionally consistent but not significant:** HMGCR +0.160 · AIFM2 +0.156 ·
SMG9 +0.149 · MTARC1 −0.156 · SMG1 −0.132

### Three findings that need stating plainly

1. **CYB5R3 itself does not move** (ρ = −0.074, q = 0.38). The project's central
   gene is flat across the entire spectrum. This is the single most important
   negative result of PHASE B.
   *It is also what the reserve model in the concept paper predicts* — an enzyme
   in catalytic excess need not change its expression; what changes is demand.
   But that prediction was written for the Perspective, not pre-registered here,
   so it is a consistency check, not a confirmation. It must not be presented as
   though the analysis confirmed it.

2. **The REDUCTION module's signal is almost entirely NQO1**, whose mean
   expression is low (mean log₂CPM 1.39, near the filter threshold) and which is
   a canonical NRF2 target. Its rise is at least as consistent with an NRF2
   oxidative-stress response as with expanded CoQ-reducing capacity. Removing it
   weakens but does not abolish the trend (ρ 0.33 → 0.26).

3. **DRAIN is not a coherent module** — Cronbach α = −0.101. MTARC1, MTARC2 and
   POR do not co-vary; they are a functional grouping, not an expression
   program. The composite score is therefore dominated by whichever member has
   the most variance (here MTARC2). Report α alongside the module result.

## 4. Internal validity (PHASE A, for the record)

| Module | genes present | Cronbach α | PC1 variance |
|---|---|---|---|
| SUPPLY | 15/15 | 0.833 | 34.1% |
| REDUCTION | 4/4 | 0.715 | 57.2% |
| DRAIN | 3/3 | **−0.101** | 49.7% |

No specified gene was missing from the matrix. PC1 of the whole expression
matrix explains 28.9%; no sample exceeded |z| = 4 on PC1.

## 5. Where this leaves the paper

The pre-specified primary hypothesis passed, on real human tissue, with a
graded trend that survives removing controls and restricting to MASH. That is a
publishable result and it was obtained without touching an outcome variable.

But the mechanism the index was built to capture is only partly borne out.
Honest framing for the manuscript:

- **Claim:** the balance between reductive supply and pathogenic drain shifts
  with MASLD severity. **Supported.**
- **Claim:** the shift is driven by expanded mevalonate/CoQ supply.
  **Not supported within the disease spectrum.**
- **Claim:** the index is CoQ-specific rather than cholesterol-driven.
  **Weakened by S2, as pre-specified.**
- **Claim:** CYB5R3 is the operative hub. **Its expression is flat.** The
  reserve model explains this, but the analysis does not test it.

## 6. Replication — three independent cohorts, added 2026-08-25

Axes were fixed in PREREGISTRATION.md §6a **before** any replication cohort's
RSI was computed, because each cohort annotates severity differently.

| Cohort | Axis | n | JT z | P (increasing) | Spearman ρ (95% CI) |
|---|---|---|---|---|---|
| GSE135251 (discovery) | group | 216 | 5.98 | 1.1 × 10⁻⁹ | 0.397 (0.270–0.509) |
| GSE130970 | fibrosis F0–F4 | 78 | 3.56 | 1.9 × 10⁻⁴ | 0.389 (0.169–0.572) |
| GSE167523 | NAFL < MASH | 98 | 2.95 | 1.6 × 10⁻³ | 0.299 (0.095–0.479) |
| **Pooled (random effects)** | — | **392** | — | — | **0.372 (0.277–0.459)** |

**Heterogeneity Q = 0.75 on 2 df, I² = 0%.** The three cohorts — different
platforms, different units, different severity axes — give effectively the same
effect size.

Pre-specified amendment item 3 (one comparison on an identical axis):
GSE135251 re-run on the fibrosis axis gives ρ = 0.394 (0.267–0.507),
against GSE130970's ρ = 0.389 on the same axis. The two agree to the third
decimal.

### The component pattern depends on the axis, and that is informative

| Cohort | Axis | supply z | reduction z | drain z |
|---|---|---|---|---|
| GSE135251 | group (fibrosis-based) | +2.56 | +3.59 | **−3.84** |
| GSE130970 | fibrosis | +0.27 | +1.56 | **−2.75** |
| GSE167523 | NAFL → MASH | **+2.45** | **+3.00** | +0.19 |

Along a **fibrosis** axis the index moves because the drain falls. Along the
**NAFL → MASH transition** it moves because supply and reduction rise, with the
drain flat. These are not contradictory: they say the balance shifts by
different means at different points of the disease. The supply-side mechanism
that motivated the study appears at the steatosis-to-steatohepatitis step; the
drain-side change accompanies fibrosis progression.

In GSE167523, CYB5R3 itself does rise (ρ = +0.294, q = 0.0045) — unlike
GSE135251, where it was flat. SMG9 is the top-ranked gene in that cohort
(ρ = +0.430, q = 1.4 × 10⁻⁴), with HMGCR, SREBF2 and AIFM2 also passing FDR.

### One near-miss that must be recorded

GSE167523 uses **retired gene symbols**: MARC1 and MARC2 for MTARC1 and
MTARC2. On the first run those two genes were simply absent, the DRAIN module
collapsed to POR alone, and the result was **ρ = −0.265 — an apparent failed
replication in the opposite direction.** Resolving the aliases from the
Synonyms column of NCBI's annotation table restored the module and the result
became ρ = +0.299.

The "no specified gene silently missing" check is what caught this. Without it
the manuscript would have reported a contradiction that does not exist. Symbol
harmonization is now applied to every symbol-keyed cohort.

NQO1 is genuinely absent from GSE167523 (FPKM ≈ 0.19), not a naming problem, so
that cohort's REDUCTION module has three of four genes. This is reported rather
than patched.

### Cohort-specific cautions

- **GSE167523** has PC1 = 55% of total variance, and module α values of 0.96
  (SUPPLY), 0.92 (REDUCTION) and 0.89 (DRAIN) — far above the other cohorts.
  A single dominant axis inflates all within-module correlations, so the high
  α there is not evidence of genuine co-regulation.
- **DRAIN coherence is cohort-dependent**: α = −0.10 in GSE135251, +0.58 in
  GSE130970, +0.89 in GSE167523. Report α with every module result.
- PC1 outliers at |z| > 4: GSM3758009 (GSE130970), GSM5106632 (GSE167523).
  No exclusion rule was pre-specified, so neither was removed.
- GSE167523 is FPKM. Rescaling each sample to a fixed total converts FPKM to
  TPM, which is what the pipeline does before the log. Length normalization
  therefore differs from the other two cohorts; because the index is z-scored
  within cohort, this does not affect the comparison.
- GSE130970 F4 has n = 2 and GSE167523 has only two ordered groups, so its
  Jonckheere–Terpstra test reduces to a Mann–Whitney test. Both were declared
  in advance.

## 7. Next

1. ~~Replication in GSE130970 and GSE167523~~ — **done, §6. Both replicate.**
2. PHASE C (GSE164760) for H2 and H3, the field-effect question.
3. PHASE D: NMD-target burden vs RSI (H5). SMG9 rising and SMG1 falling are both
   directionally consistent here but neither is significant, so H5 is still open.
4. Do **not** redefine the index. The lock holds.
