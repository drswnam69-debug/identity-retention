# PHASE E, part 1 — H6: RSI and ferroptosis signatures

**Run:** 2026-08-25 · **Index lock hash:** `7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`
**Plan:** PREREGISTRATION §6d, fixed before any outcome was read.
H6 uses expression only; no survival variable was opened for this part.

---

## 1. Gene sets, and the overlap that had to be removed

FerrDb V2 was the pre-specified first choice; its host is unreachable, so the
pre-specified fallback — MSigDB — was used. Two sets carry the direction H6
needs:

| Set | genes | dropped for overlapping the RSI panel | kept |
|---|---|---|---|
| GOBP_NEGATIVE_REGULATION_OF_FERROPTOSIS (**suppressors**) | 19 | AIFM2, GPX4, NQO1, SLC7A11 | **15** |
| GOBP_POSITIVE_REGULATION_OF_FERROPTOSIS (**drivers**) | 10 | ACSL4 | **9** |
| GOBP_FERROPTOSIS | 31 | ACSL4, AIFM2, GPX4, NQO1, SLC7A11 | 26 |
| WP_FERROPTOSIS | 64 | ACSL4, AIFM2, COQ2, GPX4, HMGCR, LPCAT3, POR, SLC7A11 | 56 |

The overlap is large — up to 8 of 64 genes — and every shared gene was removed
before scoring, as §6d requires. Without that step the test would have been
partly circular.

## 2. H6 is not supported in the direction it predicted

H6 predicted RSI **positive** with a ferroptosis-resistance signature, i.e.
positive with suppressors and **negative** with drivers.

| Cohort | suppressors ρ | drivers ρ |
|---|---|---|
| GSE135251 (n=216) | **+0.439** (P = 1.4 × 10⁻¹¹) | **+0.502** (P = 3.3 × 10⁻¹⁵) |
| GSE130970 (n=78) | +0.209 (P = 0.066) | +0.190 (P = 0.095) |
| GSE167523 (n=98) | +0.171 (P = 0.092) | +0.217 (P = 0.032) |

RSI is positively correlated with **both** arms in every cohort, and in the
discovery cohort slightly *more* with drivers than with suppressors. The
predicted dissociation does not appear.

## 3. It is not the global-expression artifact

The failure mode found in PHASE D was checked for here. Both signature scores
load heavily on each sample's mean expression (ρ = 0.63 to 0.95), so the raw
correlations were repeated with mean expression partialled out:

| Cohort | set | raw ρ | adjusted ρ | P (adjusted) | set vs mean expression |
|---|---|---|---|---|---|
| GSE135251 | suppressors | +0.439 | **+0.318** | 2.0 × 10⁻⁶ | +0.696 |
| GSE135251 | drivers | +0.502 | **+0.410** | 3.9 × 10⁻¹⁰ | +0.687 |
| GSE130970 | suppressors | +0.209 | **+0.252** | 0.027 | +0.627 |
| GSE130970 | drivers | +0.190 | **+0.277** | 0.015 | +0.777 |
| GSE167523 | suppressors | +0.171 | **+0.283** | 0.0049 | +0.952 |
| GSE167523 | drivers | +0.217 | **+0.318** | 0.0015 | +0.902 |

The associations survive adjustment, and in the two replication cohorts they
strengthen. So the relationship is real — it is simply not the relationship H6
predicted.

## 4. What this most likely means

RSI tracks the **ferroptosis regulatory apparatus as a whole**, not
ferroptosis resistance specifically. Both GOBP sets are heavily loaded with
NRF2/oxidative-stress genes (the suppressor set retains HMOX1, NFE2L2, GSTP1,
EGLN3; the stripped genes were themselves NRF2 targets), and in progressing
liver disease that whole program moves together. A signature built this way
cannot separate "more defense" from "more threat".

Two honest consequences:

1. **H6 as written is not supported.** It is reported as a failed directional
   prediction, not reinterpreted into a success.
2. The finding that RSI rises alongside ferroptosis-regulatory gene expression
   is consistent with the framework's *setting* — a liver under increasing
   lipid-peroxidation pressure — but it does not establish the direction of the
   threshold, which is what the framework actually claims. Establishing that
   needs a functional readout (lipid peroxidation, CoQH₂/CoQ), not a
   transcriptional signature.

## 5. Status

| Hypothesis | Status |
|---|---|
| H1 spectrum trend | supported; 3 cohorts, pooled ρ = 0.372, I² = 0%; survives global-axis adjustment |
| H2 tumor vs adjacent | not supported |
| H3 field effect | not supported (P = 0.060) |
| H5 NMD burden | not testable with gene-level data; open |
| **H6 ferroptosis direction** | **not supported**; a real but non-directional association |
| H4 survival | pending — awaiting the GSE76427 expression matrix |

## 6. H4 is set up and blocked only on one file

Everything else for H4 is in place and verified:

- `survival.py` — Cox (Breslow ties), Kaplan–Meier, k-sample log-rank and a
  Schoenfeld proportional-hazards check, in numpy alone. `test_survival.py`
  matches statsmodels to machine precision on coefficients, standard errors,
  hazard ratios and p-values, with and without heavy ties, and the PH test has
  both a positive and a negative control.
- **GSE76427** phenotype built from GEO metadata and asserted against GEO's own
  totals: 115 tumors, 52 adjacent, 115 patients, **23 overall-survival events**.
- GPL10558 probe map: 59 probes covering all 42 panel genes, with the 2010-era
  retired symbols resolved and PRG3 excluded.

⚠ **23 events is a hard constraint.** The pre-specified primary model — Cox
with RSI as a continuous variable — is appropriate, but the adjusted model
(RSI + age + sex + BCLC) will run at under 6 events per variable and is
underpowered. It will be reported with that warning attached rather than
treated as confirmatory.

The only missing input is `GSE76427_non-normalized.txt.gz` (47.5 MB). Chrome
began blocking repeated automatic downloads from NCBI after the earlier
cohorts, so this one needs a manual click.
