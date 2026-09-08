# PHASE E — prognosis (H4), ferroptosis (H6), and an H2 replication

Index locked 2026-08-25, hash `7a2bf934fe6e5118…`. Analysis plan fixed in
PREREGISTRATION §6d (H4/H6) and §6f (H2 replication), each recorded before the
corresponding numbers were computed. Acquisition route for GSE76427 recorded in
§6e.

---

## 1. Cohort

**GSE76427** — 115 HCC tumors and 52 patient-matched adjacent non-tumor livers,
115 patients, 23 overall-survival events. Mixed etiology, Illumina HumanHT-12
V4. Outcome-blind to this project: no prior analysis in this program has touched
it (unlike TCGA-LIHC, see §6b).

All 22 panel genes present; 42 symbols from 59 probes. Internal validity:
SUPPLY α=0.672, REDUCTION α=0.615, DRAIN α=0.755. No PC1 outliers.

---

## 2. H4 — RSI and overall survival: **not supported**

Registered direction: HR > 1.

| model | HR per SD | 95% CI | P |
|---|---|---|---|
| Cox, continuous (primary) | **0.986** | 0.625 – 1.555 | 0.95 |
| Cox, adjusted (age, sex, BCLC) | 0.951 | 0.607 – 1.491 | 0.83 |
| Log-rank across RSI tertiles | — | — | 0.83 |

Tertiles carried 8 / 7 / 8 events — no gradient. Proportional hazards holds
(Schoenfeld ρ = −0.31, P = 0.16), so the null is not an artifact of a
time-varying effect.

**The model is not blind to real prognostic signal.** In the same adjusted fit,
BCLC C/D carries **HR 4.56 (1.44 – 14.45), P = 0.0099**. A 23-event cohort that
recovers stage as a strong predictor and returns 0.99 for RSI is giving an
informative answer, not a powerless one.

**What this null does and does not exclude.** With 23 events, 80% power reaches
HR ≥ 1.79 (or ≤ 0.56) per SD. The observed CI (0.625 – 1.555) sits inside that
window: a *strong* prognostic effect is excluded, a *modest* one (HR ≈ 1.3) is
not. This is reported as a bounded null, not as proof of no association.

⚠ The adjusted model runs at 5.75 events per variable, below the conventional
10. It is reported as pre-specified and read as supporting, not primary.

**Global-axis check (standing rule, §6).** Within tumors, RSI vs mean expression
ρ = +0.145, P = 0.12 — below the materiality threshold, so the additionally
adjusted survival model was correctly not triggered.

---

## 3. H2 replication — **supported, and strongly** (§6f)

Registered direction (§1, fixed 2026-08-25): **tumor > adjacent**. In
GSE164760 this was not supported (shift +0.15, δ = +0.10, P = 0.44) —
directionally consistent but far short of significance. GSE76427 retests it in
an independent cohort with patient-matched pairs.

| test | estimate | P |
|---|---|---|
| Mann–Whitney, 115 vs 52 | Hodges–Lehmann **+1.227** (95% CI +1.002 to +1.439) | **4.6 × 10⁻¹⁵** |
| Cliff's delta | **+0.759** | — |
| Wilcoxon signed-rank, **52 patient-matched pairs** | median paired difference **+1.129** | **8.1 × 10⁻⁸** |

Mean RSI: tumor **+0.378**, adjacent **−0.835**.

**It is not the global expression axis.** Residualizing RSI on each sample's
mean expression leaves the contrast essentially intact: shift +1.163,
P = 2.7 × 10⁻¹⁵.

**Which arm carries it — the informative part.** Decomposing into the three
modules (within-cohort z, tumor vs adjacent):

| module | tumor | adjacent | P |
|---|---|---|---|
| SUPPLY (mevalonate/CoQ) | −0.009 | +0.021 | 0.49 |
| REDUCTION (CYB5R3, CYB5R1, AIFM2, NQO1) | **+0.212** | −0.469 | 1.2 × 10⁻¹¹ |
| DRAIN (MTARC1, MTARC2, POR) | **−0.276** | +0.611 | 6.1 × 10⁻¹⁴ |

The separation is produced entirely by **reduction up and drain down** — the two
terms that encode retained reductive capacity — with the biosynthetic supply arm
flat. This is the first tissue-level contrast in the project where the
reduction arm moves in the hypothesized direction, and it is the CYB5R3
"reductive reserve" axis specifically, not the mevalonate axis.

### 3a. The dedifferentiation confound — tested, and it is real but partial (§6g)

The threat named above was pre-registered as §6g and tested. A hepatocyte-identity
score **D1** was specified in advance from 22 genes, deliberately containing no
redox enzyme and no P450 (`HNF4A HNF1A FOXA1 FOXA2 NR1H4 · ALB TTR TF SERPINA1
AHSG APOH FGA FGB FGG F2 · CPS1 OTC ARG1 TAT G6PC1 PCK1 · ASGR1`), plus a
deliberately over-adjusted **D2** = D1 + `CYP2E1 CYP3A4 CYP1A2 CYP2C9`. All 26
genes were present on GPL10558 (G6PC1 via its retired name G6PC).

**The confound exists.** D1 is markedly lower in tumor (−0.150 vs +0.331,
P = 2.8 × 10⁻¹⁴). This is not a hypothetical objection.

**Primary test — paired-difference regression over the 52 matched pairs.**
Regressing ΔRSI on ΔD1, the intercept is the RSI rise at zero change in
hepatocyte identity:

| | intercept (ΔRSI at ΔD1 = 0) | 95% CI | P | fraction of the unadjusted shift retained |
|---|---|---|---|---|
| unadjusted | +1.129 (median) | — | 8.1 × 10⁻⁸ | — |
| **adjusted for D1** | **+0.568** | +0.292 to +0.843 | **1.3 × 10⁻⁴** | **50%** |
| adjusted for D2 (over-adjusted) | +0.413 | +0.120 to +0.707 | 6.8 × 10⁻³ | 37% |

Unpaired model over all 167 samples agrees: tumor coefficient +0.689
(+0.462 to +0.915), P = 1.2 × 10⁻⁸.

**The honest reading: about half the tumor rise is loss of hepatocyte identity,
and about half is not.** The surviving half is statistically robust and survives
deliberate over-adjustment with the P450s.

⚠ **The registered verdict is knife-edge and should not be leaned on.** The rule
fixed in §6g called ≥50% retention "not explained by dedifferentiation". The
observed retention is **50.3%** on the registered denominator (the paired
median) and 53.7% on the mean — i.e. it clears the threshold by three tenths of
a percentage point. The verdict label is therefore reported, but the substantive
claim is the one stated above: **halved, but robust**. No one should read
"NOT explained" as though the confound were absent.

**Module level is where this gets informative.** Repeating the same adjustment
separately:

| module | unadjusted Δ | D1-adjusted intercept | 95% CI | P | retained |
|---|---|---|---|---|---|
| REDUCTION | +0.728 | **+0.517** | +0.252 to +0.783 | 2.8 × 10⁻⁴ | **71%** |
| DRAIN | −0.712 | −0.295 | −0.521 to −0.069 | 0.012 | 41% |

**The drain fall is majority dedifferentiation; the reduction rise is not.**
This is the sharpest result in the study. The reviewer's objection was correct
about DRAIN — mARC1/mARC2/POR do largely track hepatocyte identity — and wrong
about REDUCTION, which retains 71% of its effect after adjustment and remains
highly significant. The CYB5R3-side signal is the durable one.

The paper should therefore lead on **reduction**, not on the composite, and
state the DRAIN caveat rather than bury it.

### 3b. Replicated in an independent cohort (§6h/§6i)

The §6g result above rested on one cohort, and its 50.3% composite retention sat
on the registered threshold by three tenths of a point. It has now been retested
in **GSE14520** — 445 samples, **213 patient-matched pairs**, HBV rather than
mixed etiology, Affymetrix HG-U133A rather than Illumina — under a rule fixed in
§6h before that cohort was acquired.

| | GSE76427 (§6g) | GSE14520 (§6h) |
|---|---|---|
| pairs | 52 | **213** |
| D1 lower in tumor | P = 2.8 × 10⁻¹⁴ | **P = 9.7 × 10⁻⁵⁹** |
| composite, D1-adjusted intercept | +0.568 (**50%**) | **+0.792 (63%)** |
| composite, D2-adjusted | +0.413 (37%) | +0.656 (52%) |
| REDUCTION retained | **71%** | **84%** |
| DRAIN retained | 41% | **25%** |

**The module split replicates and widens.** Under the deliberately over-adjusted
D2 model the two arms separate completely in both cohorts: REDUCTION keeps
69% / 81% of its shift and stays highly significant, while DRAIN falls to
20% / 8% and **loses significance in both** (P = 0.22, P = 0.36). That is what
the mechanism predicts — POR exists to donate electrons to those very P450s, so
D2 conditions on part of the drain itself.

So the substantive claim above is upgraded from single-cohort to **replicated in
two independent cohorts**: the tumor RSI rise is roughly half to two-thirds
independent of hepatocyte identity, and the independent part is carried by the
**reduction** arm. GSE164760 remains **uninformative** by the registered rule
(its unadjusted DRAIN shift is absent) and is counted as neither support nor
refutation.

### 3c. Narrowed by the composition test (§6j)

§3b said the module split replicates. **§6j narrows that.** After adjusting for
cell composition *and* hepatocyte identity together, the split holds decisively
in GSE14520 (REDUCTION 89% vs DRAIN 31%) but the two arms are **no longer
distinguishable** in GSE76427, where the retentions are 51.20% and 51.12% — a
margin of 0.0008. The composite rise survives in both (65% and 40% retained,
P = 1.8 × 10⁻²² and 0.014), but the *dissociation* — the study's sharpest claim —
now rests on the 213-pair cohort alone.

Composition by itself explains none of the rise (101% and 116% retained), so the
"it's infiltrating cells" objection is wrong as stated; it is the **combination**
with identity loss that bites, and only in the smaller cohort. With 52 pairs
that may be power rather than substance — GSE76427's intervals are wide enough
to contain GSE14520's split — but that is an explanation offered, not a result.

**§3b's wording must therefore change** from "the module split replicates and
widens" to: *the split is decisive in the larger paired cohort and not
demonstrable in the smaller once both confounds are removed.* Full detail in
`results/COMPOSITION_6j_RESULTS.md`.

⚠ **One limitation must travel with this number.** AIFM2 is absent from
HG-U133A, so GSE14520's REDUCTION is a 3-gene score — and AIFM2 is the
*strongest* member of REDUCTION where it is measured (GSE76427 Δz = +1.292,
P = 2.31 × 10⁻¹⁶). The 84% retention is therefore not strictly comparable to the
71% from a 4-gene REDUCTION, and the non-comparability runs against the
hypothesis rather than for it. (An earlier version of this paragraph asserted the
opposite; see PREREGISTRATION §6k, 2026-08-27.) It belongs in the limitations
section, not a footnote.

---

## 4. H6 — ferroptosis signatures: unchanged from the earlier run

RSI-panel genes removed from every set before scoring (4–8 genes per set).

| cohort | NEG. reg. | POS. reg. | GOBP_FERROPTOSIS | WP_FERROPTOSIS |
|---|---|---|---|---|
| GSE135251 (n=216) | +0.439 (1.4e-11) | +0.502 (3.3e-15) | +0.551 (1.5e-18) | +0.350 (1.3e-07) |
| GSE130970 (n=78) | +0.209 (0.066) | +0.190 (0.095) | +0.256 (0.023) | +0.195 (0.087) |
| GSE167523 (n=98) | +0.171 (0.092) | +0.217 (0.032) | +0.209 (0.039) | +0.075 (0.46) |

Registered direction was **positive with suppressors, negative with drivers**.
RSI is positive with *both*. The association is real and survives global-axis
adjustment, but it is **non-directional** — consistent with an NRF2-loaded
regulatory program rather than with ferroptosis resistance as specified.
**H6 not supported directionally.**

⚠ **Coverage limitation, stated rather than omitted.** GSE164760 and GSE76427
were prepared with panel-restricted probe maps (42 genes), so ferroptosis sets
score only 0–2 genes there and both cohorts are **skipped**, not silently
absent. H6 rests on the three RNA-seq cohorts only.

---

## 5. Where the six hypotheses now stand

| | hypothesis | verdict |
|---|---|---|
| H1 | RSI rises across the MASLD severity gradient | **supported, replicated** — pooled ρ = 0.372 (0.277–0.459), I² = 0%, n = 392 |
| H2 | RSI higher in tumor than adjacent liver | **supported on replication** — GSE76427 δ = +0.76, paired P = 8.1e-08 (null in GSE164760) |
| H3 | field effect in adjacent liver | not supported (P = 0.060) |
| H4 | RSI predicts overall survival | **not supported** — HR 0.99 per SD; strong effects excluded, modest ones not |
| H5 | NMD burden vs RSI | **not testable** — the NMD proxy tracks mean expression at ρ 0.82–0.98 |
| H6 | RSI vs ferroptosis resistance | not supported directionally |

Per §7, the index is **not** redefined.

---

## 6. What the study now is

A staging-associated index, not a prognostic one. RSI tracks MASLD severity
consistently across three independent cohorts, and separates tumor from adjacent
liver decisively in a fourth — but carries no independent survival information
once stage is in the model. The honest framing is a **tissue-state** index whose
tumor signal is carried by the reduction-minus-drain arm, with the
dedifferentiation confound named as the primary threat to that interpretation.

## 7. Next

1. ~~Differentiation-adjusted H2~~ — **done** (§3a) and ~~replicated~~ —
   **done** (§3b, GSE14520, 213 pairs). The remaining cohort-level gap is a
   genome-wide re-prep of GSE164760, which would unlock H6 there and let §3a be
   repeated on a full matrix.
3. H5 still needs transcript- or junction-level quantification (check recount3
   for `SRP217231` first).
4. A larger adjacent-tissue series to give H3 adequate power.
