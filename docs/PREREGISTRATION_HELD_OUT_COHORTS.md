# Predictions recorded before the two held-out cohorts were obtained

This file holds, in the order they were written, the predictions I recorded for
the two cohorts that were not used to define the metabolic module, together with
the amendments I made to them and the reason for each.

## What kind of record this is, and what it is not

These are my own dated records. They were written in a working document as the
analysis proceeded, and they are deposited here so that a reader can see exactly
what was predicted, with what thresholds and what decision rule, and compare it
with what was reported.

They are **not** a third-party registration. This archive was released on
30 September 2026, after both cohorts had been obtained, so the archive's own
timestamp does not establish that the predictions preceded the data. What can be
checked from the archive is weaker and still worth stating: `code/77_holdout_validation.py`
contains the numeric thresholds and the decision rule as constants, and its
`load_holdout` function, the only part that depends on the file layout of the
validation cohort, is marked in its own docstring as the one piece written after
the files arrived. Anyone re-reading that file can see which parts could not have
been tuned to the answer, because they do not mention the cohort at all.

I am stating this plainly rather than letting the word "pre-registered" carry
more weight than the evidence supports.

---

## 1. First held-out cohort, recorded 29 September 2026, before the files were obtained

### What was being tested

The claim under test was that a module of urea cycle and gluconeogenic enzymes
(*CPS1*, *OTC*, *ARG1*, *TAT*, *G6PC1*, *PCK1*) is a hepatocyte identity
covariate that means the same thing in a liver transcriptome and a liver
proteome. The module had been identified by looking at five datasets. Without a
held-out cohort the claim rests on the same data that produced it.

### Cohort

Yi X, Zhu J, Liu W, et al. Proteome landscapes of human hepatocellular carcinoma
and intrahepatic cholangiocarcinoma. *Mol Cell Proteomics* 2023;22:100604.
ProteomeXchange via iProX, PXD043265. Training cohort of 56 patients with paired
tumor and adjacent benign tissue, of whom 41 have hepatocellular carcinoma.
Wuhan, China. Data-independent acquisition on a Q Exactive HF.

Different investigators, city, acquisition mode and quantification software from
the Gao 2019 and Jiang 2019 proteomes used to define the module.

### Sample identity gate, which has to pass first

Of the proliferation proteins *MKI67*, *PCNA*, *TOP2A*, *RRM2*, *TYMS*, *MCM2*,
*CCNB1* and *AURKA*, at least three must be quantified in 10 or more pairs and
their mean paired shift must be positive in tumor. If it is not, the tumor and
adjacent labels or the pairing cannot be relied on, and the matrix is reported as
not testable and no contrast is computed.

None of these proteins belongs to the panel under test. The panel's own direction
cannot serve as this check, because the panel is the object of the study.

### Predictions

| | Prediction | What a failure means |
|---|---|---|
| H1 | The module's mean paired shift is below zero and its surviving share is at least 0.9 | The module claim is withdrawn |
| H2 | The ten plasma proteins have a mean paired shift above zero | The two-opposing-axes framing is cohort specific and is withdrawn |
| H3 | The module's mean sits below the 10th percentile of the quantified non-secreted proteins | The extremity claim is withdrawn |
| H4 | The 22-gene composite's Wilcoxon *P* is not significant, or its sign is opposite to the transcriptome sign | The cancellation claim is weakened |
| H5 | If the matrix is on an abundance scale, the module's fall is not explained by abundance-matched neighborhoods | The abundance explanation is not excluded |

**Decision rule.** H1 and H3 must both hold for the module claim to survive. H2
decides the framing. If H1 or H3 fails, this direction of work is closed.

### Amendment 1, 29 September 2026, after opening the matrix and before computing any contrast

Three things were fixed after seeing the file's structure and before computing
H1 to H5.

**A. Three quantification matrices exist, so the choice between them was fixed
in advance.** The supplement carries protein tables from three DIA packages built
from the same raw files. The primary matrix is the one with the highest
completeness, which is measurable without computing any contrast, and the others
are reported as sensitivity analyses. Completeness is 95.6% for EncyclopeDIA,
62.7% for Spectronaut and 55.4% for DIA-NN, so EncyclopeDIA is primary. It is
also the only one of the three that passes the sample identity gate.

**B. Gene symbol mapping.** The matrix is indexed by UniProt accession. UniProt's
web service was unreachable from the analysis environment, so accessions were
mapped through the `Protein IDs` and `Gene names` columns of the Jiang MaxQuant
`proteinGroups.txt` already deposited in this archive. The legacy symbol *G6PC*
maps to *G6PC1*. *G6PC1* is nevertheless absent from the primary matrix and is
reported as not assessable in this cohort; the module is scored on the five
members that are quantified, as the wording "those quantified" already provided
for. *ALB* is also absent from the primary matrix, so the plasma class is scored
on nine members.

**C. The target population was fixed before it was known whether it could be
identified.** The cohort mixes 41 hepatocellular carcinoma, 12 cholangiocarcinoma
and 3 mixed carcinoma patients. If the supplementary sample table gave a
per-patient diagnosis, the analysis would run on the 41 hepatocellular carcinoma
pairs. If it did not, the analysis would run on all 56 pairs and the population
would be reported as mixed primary liver cancer, with the note that mixing in
cholangiocarcinoma biases H1 in its own favor. Classifying patients from the data
using biliary markers was ruled out as a primary analysis. The table did give the
diagnosis, so the first branch applied.

---

## 2. Second held-out cohorts, recorded 30 September 2026, before the files were obtained

### Why this is a different test

Intrahepatic cholangiocarcinoma arises from biliary epithelium, and the adjacent
tissue is liver parenchyma. Hepatocyte identity should therefore collapse more
completely than in hepatocellular carcinoma. Cholangiocytes do not synthesize
albumin or fibrinogen.

**The plasma protein prediction therefore runs the opposite way from the first
cohort.** If plasma proteins rise in hepatocellular carcinoma proteomes because
the tumor cells are hepatocytes that still secrete, they should fall here.

### Cohorts

Werner T, Thiery J, Budau KL, et al. Proteomic characterization of intrahepatic
cholangiocarcinoma identifies risk-stratifying subgroups and EIF4A1 as a
therapeutic target. *Nat Commun* 2026;17:2741. Two independent cohorts, MSKCC-ICC
and UKF-ICC, from the United States and Germany, measured by DIA-PASEF and
DIA-NN. MassIVE MSV000095336.

### Predictions

The same sample identity gate has to pass first, judged separately in each cohort.

| | Prediction | Interpretation |
|---|---|---|
| G1 | The module's mean shift is below zero with a surviving share of at least 0.9, in both cohorts | The module reads hepatocyte loss across tumor types |
| G2 | The module falls further than in the three hepatocellular cohorts | Consistent with more complete hepatocyte loss |
| G3 | The plasma proteins have a mean shift below zero, the opposite of the first cohort | Supports residual tumor secretion as the cause of their rise in hepatocellular carcinoma |
| G4 | The gap between the plasma and metabolic classes is smaller than in the hepatocellular proteomes | The same reading, quantified |
| G5 | The module sits below the 10th percentile of the non-secreted proteins, in both cohorts | Extremity is independent of tumor type |

### What each outcome would mean, fixed in advance

* If G1 and G5 hold in both cohorts, the module works across tumor types.
* If G3 holds, the mechanistic reading of the plasma proteins' rise in
  hepatocellular carcinoma gains its first empirical support.
* **If G3 fails, that is, if plasma proteins rise here too, then their rise has
  nothing to do with hepatocyte origin.** Retention of plasma protein in tumor
  tissue becomes the better explanation, that reading carries over to
  hepatocellular carcinoma, and the residual secretion explanation is removed
  from the paper. This outcome is to be reported as it falls.
* If G1 fails, the module is specific to hepatocellular carcinoma, the claim is
  narrowed to that, and the first cohort's result is not retracted.
* G2 and G4 are directional only. A miss is reported and retracts nothing.

G3 failed, and the third bullet was applied.

---

## 3. What was run and what was not, at each stage

At the point Amendment 1 was written, only the premise checks had been run:
matrix diagnostics, missingness, the sample identity gate, and which of the genes
under test were quantified. H1 to H5 had not been computed.

At the point Amendment 3 was written (see `PROTOCOL_DEVIATIONS.md`), the same was
true of the cholangiocarcinoma cohorts: the gates and the coverage had been
examined and G1 to G5 had not been computed.
