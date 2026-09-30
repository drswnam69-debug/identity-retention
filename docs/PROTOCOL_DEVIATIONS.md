# Protocol deviations

One deviation was made from the predictions recorded in
`PREREGISTRATION_HELD_OUT_COHORTS.md`. It is described here in full, including
what was already known when the decision was taken.

## Deviation 1. The sample identity gate was replaced in the two cholangiocarcinoma cohorts

**Recorded 30 September 2026, after the gate failed and before G1 to G5 were
computed.**

### What happened

The pre-registered gate requires at least three quantified proliferation
proteins with a positive mean paired shift. In both cohorts only two of the eight
clear the 10-pair threshold.

| | MSKCC-ICC, 64 pairs | UKF-ICC, 61 pairs |
|---|---|---|
| *PCNA* | 63 pairs | 58 pairs |
| *MCM2* | 51 pairs | 19 pairs |
| *MKI67* | 0 | 0 |
| *TOP2A* | 0 | 0 |
| *RRM2* | 0 | 2 |
| *TYMS* | 0 | 0 |
| *CCNB1* | 0 | absent |
| *AURKA* | 0 | absent |

The mean paired shift of the two that are quantified is +1.644 and +1.623. The
gate failed on the count, not on the direction.

### Why it failed

The fault is in the gene list I wrote, not in the data. I filled it with
low-abundance nuclear proteins, and tissue proteomics by mass spectrometry rarely
quantifies those. I first checked whether the accession-to-symbol mapping was the
limiting factor, by enlarging the map from 9,412 to 9,477 accessions using the
`Protein` and `Gene.names` columns of the source article's own supplementary
tables; the result did not change, so the mapping was not the cause.

The genes under test are, by contrast, well covered in both cohorts: four of the
six module members and all ten plasma proteins.

### What was changed

The gate was replaced by *KRT19*, *KRT7*, *EPCAM*, *S100P*, *PCNA* and *MCM2*,
with the same requirement of at least three quantified members and a positive
mean paired shift. The reasoning is that a cholangiocarcinoma against adjacent
liver must show biliary epithelial markers rising, a direction fixed by the
biology without looking at the data, and that none of these proteins belongs to
any class of the panel under test. Keratins are abundant and are reliably
quantified in tissue proteomes.

The replacement gate passes in both cohorts, with six members at a mean of +1.619
and five at +1.578, every member positive.

### Consequence for the status of the result

**Both cholangiocarcinoma cohorts are reported as exploratory, not confirmatory,
and they are labeled that way wherever they appear.** The gate was changed after
it failed, and that is what decides the status, not how convincing the
replacement gate's numbers turned out to be. The confirmatory held-out test is
the Yi 2023 cohort alone.

The result was reported in full after the change, including the prediction that
failed. G3 predicted that plasma proteins would fall in these tumors and they
rose.

### What I would do differently

A gate for a proteomic cohort should be chosen for detectability as well as for
direction, and detectability can be judged without looking at the cohort under
test. A list of proliferation and epithelial markers that are reliably quantified
in tissue proteomes should be fixed once and reused, rather than assembled from
textbook proliferation markers that a mass spectrometer does not see.

---

## Corrections to my own analysis code, for the record

**The abundance control summary statistic, corrected 30 September 2026.**
`code/77_holdout_validation.py` first reported the share of the module's fall
explained by abundance as the mean of the per-member ratios. That is the wrong
quantity: a member whose observed shift is near zero, *TAT* in this case, returns
a ratio of several hundred percent and drags the average with it, giving 0.937.
The quantity that answers the question is the ratio of the means, which is 0.175.
The script now reports the corrected value and keeps the superseded one under
`superseded_mean_of_per_member_ratios` rather than deleting it.

**A rounding that favored the result, corrected 30 September 2026.** A draft
stated that the module's surviving share is "at least 0.954" in all eight
datasets. The minimum is 0.95399, so the stated bound was the true value rounded
up in the claim's favor. It now reads 0.95. The numeric verifier
`code/80_verify_MCP.py` carries a check that the text does not reintroduce it.
