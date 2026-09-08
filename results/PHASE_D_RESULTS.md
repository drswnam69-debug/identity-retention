# PHASE D results — H5, NMD-target burden versus RSI

**Run:** 2026-08-25 · **Index lock hash:** `7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`
**Status: EXPLORATORY, and the test does not settle H5.** See PREREGISTRATION §6c.
Survival remains unopened.

---

## 1. Two things were fixed before computing anything

**The registered direction of H5 was wrong.** H5 read: *"NMD-target burden (a
proxy for SMG1 activity) is inversely correlated with RSI (ρ < 0)."* Burden is
the accumulation of NMD substrates, which happens when NMD activity is **low**
— so burden is a proxy for the *inverse* of activity, not for activity. The
framework (low SMG1 activity → higher biosynthetic set-point → higher RSI)
predicts **ρ > 0**. Because the registered sign is an error rather than a
prediction, H5 was tested **two-sided**, with both wordings on the record.

**The gene set was fixed and checked.** NMD Consensus set, 130 genes
(Palou-Márquez G & Supek F, *Genome Biol* 2025;26:316), from Zenodo
doi:10.5281/zenodo.15974216. Verified before use: **no overlap with the RSI
panel**, and **no NMD machinery genes** — so the score is not circular with
either the index or the upstream module.

## 2. The nominal result

Score = mean within-cohort z of the set. Higher = more substrate accumulation =
lower inferred NMD activity.

| Cohort | genes matched | ρ (RSI) | 95% CI | P |
|---|---|---|---|---|
| GSE135251 | 128/130 | **+0.259** | +0.122 to +0.386 | 1.2 × 10⁻⁴ |
| GSE130970 | 123/130 | −0.003 | −0.238 to +0.232 | 0.98 |
| GSE167523 | 120/130 | +0.150 | −0.062 to +0.349 | 0.14 |
| **Pooled** | — | **+0.160** | **+0.010 to +0.303** | I² = 45% |

GSE164760 was skipped: its matrix retains only the 42 panel genes, so none of
the set is present.

The pooled interval barely excludes zero, the sign matches the framework rather
than the registered wording, and only the discovery cohort is significant.
Sensitivities N1 (drop 10 mitochondrial/OXPHOS members) and N2 (drop 4
ISR/stress members) change almost nothing.

**None of this should be believed, for the reason in §3.**

## 3. The proxy fails its own validity checks

| Cohort | score vs **mean expression** | vs supply | vs reduction | vs drain |
|---|---|---|---|---|
| GSE135251 | **+0.815** | +0.511 | +0.659 | +0.308 |
| GSE130970 | **+0.921** | +0.479 | +0.423 | +0.469 |
| GSE167523 | **+0.979** | +0.828 | +0.894 | +0.793 |

The score correlates with each sample's **average expression level** at ρ =
0.82 to 0.98. It is very largely a re-measurement of a global expression axis,
not a readout of NMD activity. That single fact explains the rest: it
correlates positively with *all three* modules including the drain, and because
RSI subtracts the drain, most of the shared component cancels — leaving the
small residual ρ ≈ 0.16 that the primary test picked up.

**N3, the pre-specified orthogonal check, fails.** If the score measured NMD
activity it should correlate *negatively* with NMD machinery expression (more
machinery → more decay → less burden). Instead:

| Cohort | SMG1 | SMG8 | SMG9 | UPF1 | UPF2 |
|---|---|---|---|---|---|
| GSE135251 | −0.22 | −0.28 | +0.76 | +0.50 | −0.01 |
| GSE130970 | +0.35 | +0.54 | +0.52 | +0.41 | +0.53 |
| GSE167523 | +0.69 | +0.83 | +0.78 | +0.81 | +0.75 |

Mostly positive, and incoherent between cohorts.

### Conclusion on H5

**H5 is not tested by this analysis, and remains open.** The gene-level proxy
does not behave like an NMD readout. This is what PREREGISTRATION §6c predicted
on methodological grounds; the diagnostics now show it empirically. The
correlation with RSI is reported for completeness and must not be cited as
support for the SMG1 arm of the framework.

A real test needs transcript- or junction-level quantification, so that each
gene's NMD-target isoform can be contrasted against its matched control isoform
— which is what the source method actually does, and what collapsing to gene
level destroys.

## 4. The failure produced something worth more than the test

The same diagnostic, applied to the index itself, answers a question that
should have been asked earlier: **is the RSI also just re-reading global
expression?**

| Cohort | supply | reduction | drain | → **RSI** |
|---|---|---|---|---|
| GSE135251 | +0.45 | +0.58 | +0.17 | **+0.32** |
| GSE130970 | +0.40 | +0.48 | +0.42 | **+0.02** |
| GSE167523 | +0.80 | +0.89 | +0.84 | **+0.09** |
| GSE164760 | +0.75 | +0.25 | +0.36 | **−0.02** |

Every module is heavily loaded on the global axis. **The index is not** — the
drain subtraction cancels it, which is the balance construction doing exactly
what it was designed to do. This is the strongest argument yet for a contrast
index over any single module or single gene.

GSE135251 retains ρ = +0.32, so H1 was re-tested with the sample's mean
expression partialled out:

| Cohort | raw ρ | adjusted ρ | P (adjusted) |
|---|---|---|---|
| GSE135251 | +0.397 | **+0.345** | 2.2 × 10⁻⁷ |
| GSE130970 | +0.389 | **+0.388** | 4.9 × 10⁻⁴ |
| GSE167523 | +0.299 | **+0.289** | 4.1 × 10⁻³ |

**H1 survives in all three cohorts.** The main finding is not a global
expression artifact. This check is now part of `02_compute_rsi.py` and prints
for every cohort.

One caveat carried forward to PHASE C: mean expression does differ between the
field-effect groups in GSE164760 (adjacent lower than NASH, P = 0.023), though
RSI is essentially uncorrelated with it there (ρ = −0.02).

## 5. Status

| Hypothesis | Status |
|---|---|
| H1 spectrum trend | **supported**; 3 cohorts, pooled ρ = 0.372, I² = 0%; survives global-axis adjustment |
| H2 tumor vs adjacent | not supported |
| H3 field effect | not supported (P = 0.060) |
| **H5 NMD burden vs RSI** | **not testable with gene-level data; open** |
| H4, H6 | PHASE E, not opened |

## 6. Next

1. **PHASE E** — H4 and H6. ICGC LIRI-JP is the primary prognostic cohort;
   TCGA-LIHC is not outcome-blind (PREREGISTRATION §6b).
2. H5 stays open. If it is worth pursuing, it needs transcript-level
   requantification, and that is a separate piece of work.
3. Note for the manuscript: a prior script on this machine
   (`SMG1_NMD_activity_signature.R`, for a different SMG1/LIHC manuscript) uses
   an ad-hoc NMD-substrate signature for the same purpose and carries the same
   limitation. It also includes SMG5 and SMG7, which are NMD-autoregulated and
   therefore partly circular, and SLC7A11 and SESN2, which are in this study's
   effector module. Whatever is published should use one defensible set, and
   the two manuscripts should not report conflicting NMD proxies.
