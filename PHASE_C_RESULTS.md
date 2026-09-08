# PHASE C results — GSE164760, the HCC transition and the field effect

**Run:** 2026-08-25 · **Index lock hash:** `7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`
**Cohort:** 170 Affymetrix HG-U219 arrays — healthy 6, NASH without HCC 74,
cirrhotic 8, non-tumor NASH adjacent to HCC 29, NASH-HCC tumor 53
(group sizes read from GEO and verified before analysis)
**Analyzed as an independent arm**, never merged with the RNA-seq cohorts.
Survival remains unopened.

---

## 1. Both pre-specified hypotheses failed

| | Hypothesis | Cliff's δ | Hodges–Lehmann shift | P |
|---|---|---|---|---|
| **H2** | tumor > adjacent non-tumor | +0.103 | +0.150 | **0.44** |
| **H3** | adjacent > NASH without HCC (field effect) | +0.239 | +0.285 | **0.060** |

H2 is flatly not supported. H3 points in the predicted direction with a
small-to-moderate effect but does not reach significance at α = 0.05.

The paired sensitivity — 21 tumor/adjacent pairs inferred from consecutive
source identifiers, since GEO declares no patient ID — also fails
(median difference +0.441, P = 0.13).

Both results are robust to the normalization choice. The array data are on a
linear scale and are not equalized across samples (total-intensity CV ≈ 8%),
so log₂ followed by quantile normalization was used. Repeating everything with
global scaling instead gives H2 δ = +0.110, P = 0.41 and H3 δ = +0.228,
P = 0.072 — the same conclusion.

## 2. The pattern that is there, stated carefully

The pre-specified contrasts are not the only ones, and the others are
informative:

| Contrast | Cliff's δ | P |
|---|---|---|
| NASH vs healthy | −0.004 | 0.99 |
| **cirrhotic vs NASH** | **+0.696** | **0.0013** |
| **tumor vs NASH without HCC** | **+0.322** | **0.0020** |
| tumor vs adjacent (H2) | +0.103 | 0.44 |

Group medians: healthy −0.10, NASH −0.18, cirrhotic +0.50, adjacent +0.21,
tumor +0.23.

Read together: RSI is clearly higher in the **cancer-bearing or cirrhotic
liver** than in NASH without HCC, but **adjacent non-tumor tissue and tumor are
indistinguishable from each other**. One reading is that H2 fails *because* the
field is already shifted — the adjacent liver has moved most of the way to the
tumor's value. That is the field-effect story, and it is consistent with the
data.

**But H3, the test written to detect exactly that, returned P = 0.060.** The
supporting contrast (tumor vs NASH without HCC) was listed in the analysis plan
as context, not as a test of H3. So the honest statement is: *the pattern is
consistent with a field effect, and the pre-specified test of it did not reach
significance.* It should be reported that way, not as a confirmed field effect.

## 3. Component and gene level

Only the REDUCTION module moves, and it moves in both contrasts:

| Module | H2 δ (P) | H3 δ (P) |
|---|---|---|
| SUPPLY | +0.019 (0.89) | +0.017 (0.90) |
| **REDUCTION** | **+0.320 (0.017)** | **+0.294 (0.021)** |
| DRAIN | −0.005 (0.97) | −0.135 (0.29) |

⚠ Six component tests were run. Under Bonferroni (α = 0.0083) neither survives.
These are exploratory.

### Gene level, field-effect contrast (BH-FDR within the contrast)

- **CYB5R3 δ = +0.470, q = 0.0029** — higher in adjacent non-tumor liver than
  in NASH without HCC.
- **AIFM2 (FSP1) δ = −0.410, q = 0.0081** — lower in adjacent liver.

Nothing else passes. The two significant genes move in **opposite** directions,
which is why the REDUCTION module — which contains both — shows only a diluted
composite signal, and why the full index shows less still.

**The CYB5R3 result deserves emphasis and caution in equal measure.** CYB5R3 is
flat across the entire MASLD spectrum in two independent cohorts and two
independent pipelines (§PHASE_B, and the earlier work disclosed in
PREREGISTRATION §6b). Here, in the one comparison that asks about tissue
adjacent to cancer, it moves — and it is the largest gene-level effect in the
panel. That is the first positive signal for CYB5R3 anywhere in this project.

It is a single gene in a single cohort, in a contrast whose composite primary
test did not reach significance, and it was not a pre-specified gene-level
hypothesis. It is a finding to pursue, not to claim.

Note on scope: the disclosed earlier CYB5R3 work never used GSE164760, so this
result belongs to the present manuscript rather than to that one.

## 4. Internal validity, and why it is weaker here

| Module | genes present | Cronbach α | PC1 variance |
|---|---|---|---|
| SUPPLY | 15/15 | 0.564 | 27.7% |
| REDUCTION | 4/4 | 0.133 | 34.6% |
| DRAIN | 3/3 | 0.095 | 53.2% |

Module coherence is markedly lower than in the RNA-seq cohorts (SUPPLY α was
0.83 in both GSE135251 and GSE130970). Array measurement is noisier, and the
index was computed from 121 probesets collapsed to 42 genes by taking the
highest-expressing probeset per gene. Lower coherence means lower power, and it
is part of why H3 may have missed.

## 5. Identifier work that was necessary and is easy to get wrong

GPL13667's annotation dates from 2010 and uses retired symbols. Matching the
locked gene list against it naively loses genes silently:

| Array symbol | Current symbol |
|---|---|
| MOSC1 / MOSC2 | MTARC1 / MTARC2 |
| C17orf71 / C19orf61 | SMG8 / SMG9 |
| SC5DL / SC4MOL | SC5D / MSMO1 |
| SCD | SCD1 |

Without this, the whole DRAIN module would have been absent — the same failure
mode that briefly inverted GSE167523's result in PHASE B.

⚠ **PRG3 was deliberately excluded.** Some databases list PRG3 as an AIFM2
alias, but this array annotates probeset 11737402_at as *proteoglycan 3*, a
different gene. Accepting the alias would have injected an unrelated gene into
the REDUCTION module. AIFM2 has its own probesets, so nothing was lost.

## 6. Where this leaves the study

| Hypothesis | Status |
|---|---|
| H1 spectrum trend | **supported**, three cohorts, pooled ρ = 0.372, I² = 0% |
| H2 tumor vs adjacent | **not supported** |
| H3 field effect | **not supported** (P = 0.060, direction as predicted) |
| H5 NMD burden vs RSI | not yet tested |
| H4, H6 | PHASE E, not opened |

H3 was one of the two hypotheses carrying the study's originality. It did not
reach significance, and the pre-specification says explicitly that the index
will not be redefined to rescue it. The lock holds.

What the paper can honestly claim after PHASE C: the supply-versus-drain
balance tracks MASLD severity, robustly and reproducibly; it is further
elevated in cirrhotic and cancer-bearing liver; and adjacent non-tumor tissue
resembles tumor more than it resembles NASH without HCC — a pattern consistent
with a field effect that this study was not powered to establish.

## 7. Next

1. PHASE D — H5, NMD-target burden versus RSI. Still open, and now the last
   untested hypothesis before the prognostic phase.
2. PHASE E — H4 and H6. **ICGC LIRI-JP should be the primary prognostic
   cohort**, because TCGA-LIHC's CYB5R3 survival relationship has already been
   seen (PREREGISTRATION §6b) and LIRI-JP has not.
3. Consider whether a larger adjacent-tissue series exists that could test H3
   with adequate power. GSE164760's 29 adjacent samples are the constraint.
