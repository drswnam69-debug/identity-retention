# §6m — cell-type attribution and protein corroboration: results

**Run:** 2026-08-31 · **Plan:** PREREGISTRATION §6m, recorded before any value below was read.
**Source:** Human Protein Atlas — liver single-cell RNA (11 liver cell types, 6,404 cells)
and normal-liver immunohistochemistry. **Index lock hash unchanged.**

This is **supporting evidence about the genes**. Per §6j-S and §6m it **cannot change
any verdict in §6g–§6l**, and it has not been fed back into the index.

---

## 1. §6m-A — hepatocyte share, by the rule fixed in advance

Registered statistic: hepatocyte nCPM ÷ summed nCPM across the 11 liver cell types.
Registered threshold: ≥ 0.50 for the majority of a module's genes.

| Gene | Module | Hepatocyte nCPM | **Hepatocyte fraction** | Weighted¹ | Top cell type |
|---|---|---|---|---|---|
| MTARC1 | DRAIN | 1.0 | 0.385 | 0.768 | Monocytes (1.2) |
| **MTARC2** | DRAIN | 155.3 | **0.549** | 0.919 | **Hepatocytes** |
| **POR** | DRAIN | 709.7 | **0.602** | 0.909 | **Hepatocytes** |
| AIFM2 | REDUCTION | 31.0 | 0.274 | 0.778 | Plasma cells (34.8) |
| CYB5R1 | REDUCTION | 49.7 | 0.131 | 0.533 | Kupffer cells (60.5) |
| CYB5R3 | REDUCTION | 3.7 | 0.041 | 0.174 | Vascular endothelium (30.1) |
| NQO1 | REDUCTION | 1.4 | 0.006 | 0.049 | Cholangiocytes (119.1) |

¹ Secondary and descriptive only, not the registered statistic: nCPM weighted by the
atlas's own cell counts, which approximates each type's contribution to a bulk sample
**given this atlas's cell proportions**.

### Verdict under the registered rule

- **DRAIN — premise SUPPORTED.** 2 of 3 genes clear 0.50; both *MTARC2* and *POR* have
  hepatocytes as their single highest cell type.
- **REDUCTION — premise NOT SUPPORTED.** 0 of 4 genes clear 0.50.

§6m fixed in advance that a split verdict is reported as a split verdict. It is.

---

## 2. §6m-B — protein in normal liver, and what it says about §6m-A

| Gene | Module | Hepatocyte IHC | Cholangiocyte IHC | Hepatocyte RNA fraction |
|---|---|---|---|---|
| MTARC1 | DRAIN | Medium | Not detected | 0.385 |
| MTARC2 | DRAIN | **High** | Medium | 0.549 |
| POR | DRAIN | **High** | Low | 0.602 |
| AIFM2 | REDUCTION | Not detected | Not detected | 0.274 |
| CYB5R1 | REDUCTION | Medium | Medium | 0.131 |
| CYB5R3 | REDUCTION | Medium | Medium | 0.041 |
| NQO1 | REDUCTION | Not detected | Not detected | 0.006 |

**Two findings.**

1. **DRAIN is hepatocyte-dominant at protein level too — 3 of 3 genes stain more
   strongly in hepatocytes than in cholangiocytes.** This independently corroborates
   the premise on which the entire D1 adjustment rests: that *MTARC1*, *MTARC2* and
   *POR* are hepatocyte-identity genes. That premise is the reason the identity
   adjustment was pre-registered at all, and it now has protein-level support.

2. **The single-cell result for two REDUCTION genes is contradicted by HPA's own
   protein data.** *CYB5R3* and *CYB5R1* stain at **Medium** in hepatocytes while
   their single-cell hepatocyte fractions are 0.041 and 0.131. *MTARC1* is starker:
   **Medium** hepatocyte protein against 1.0 nCPM — near-undetected — in the RNA
   atlas. Hepatocytes are large, fragile and RNA-rich, and are well known to be
   under-recovered and ambient-RNA-contaminated in droplet single-cell liver data.
   The RNA/protein discordance here is direct, within-source evidence of that.

**§6m-B's registered primary comparison — tumour versus normal staining — was NOT
performed.** HPA's current page structure does not expose the tumour staining table,
and §6m forbids substituting another source. Recorded as not performed, exactly as
§6j-S and §6m allowed.

---

## 3. What this does and does not change

**Does not change:**

- No verdict in §6g–§6l. The composite and module results, the identity and
  composition adjustments, and the §6l POR sensitivity all stand untouched.
- In particular it does not reopen the composition question. The C1 adjustment
  already answered that empirically **in the study's own cohorts**: non-parenchymal
  content explains none of the tumour rise (101% and 116% retained). A gene being
  expressed in non-parenchymal cells in *normal* liver is not the same claim as the
  tumour rise being *caused* by a shift in cell composition, and the latter was
  tested and rejected.

**Does change, and is now reported:**

- The DRAIN module's hepatocyte-identity premise is confirmed at two levels.
- The REDUCTION module's cell-type attribution is **not** confirmed by the registered
  rule. *NQO1* in particular is cholangiocyte- and endothelium-dominant in normal
  liver by RNA and undetected in hepatocytes by protein; it is the weakest member of
  the module on this evidence, and that belongs in the limitations.
- *AIFM2* and *NQO1* are **not detected in normal hepatocytes at protein level**. That
  is compatible with the REDUCTION arm being *induced* during transformation rather
  than constitutively present — which is what H2 reports — but this study cannot test
  that reading without the tumour staining that could not be obtained.

## 4. Honest ceiling

One atlas, one liver, 6,404 cells, ordinal antibody staining on a handful of donors.
This is orientation, not validation. The manuscript's principal limitation — no new
data and no orthogonal protein-level validation under this study's own control —
remains true and remains stated.
