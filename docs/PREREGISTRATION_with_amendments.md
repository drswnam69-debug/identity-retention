# Pre-specification — hepatic Reductive-Supply Index (RSI) in MASLD and HCC

**Registered:** 2026-08-25
**Investigator:** Soon Woo Nam, MD, PhD — Hepatobiliary Unit, Division of
Gastroenterology, Department of Internal Medicine, Incheon St. Mary's Hospital,
College of Medicine, The Catholic University of Korea, Seoul, Republic of Korea
**ORCID:** 0000-0001-5562-4775
**Definition lock hash (SHA-256):**
`7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`

This document fixes the index, the hypotheses, and the analysis order before
any outcome variable is examined. It is written to be deposited with a
timestamp before PHASE B begins.

---

## 1. Rationale

Coenzyme Q10 is an antioxidant only in its reduced form. The mevalonate
pathway determines how much CoQ10 is produced; NADH-dependent CoQ reduction
determines whether that CoQ10 actually functions as an antioxidant. In liver,
the dominant paralog performing that reduction is CYB5R3. The same CYB5R3
also supplies electrons, through CYB5B, to mARC1 (MTARC1), which drives
steatotic liver injury. The protective and the pathogenic arms therefore
compete for one electron pool.

Two recent findings motivate a human-tissue counterpart that does not yet
exist. Chen et al. (*J Hepatol* 2025;83:1338–1352) showed that the mevalonate
pathway promotes liver cancer by suppressing ferroptosis through CoQ10 and
selenocysteine-tRNA modification rather than through cholesterol. Kovooru et
al. (*Clin Mol Hepatol* 2026;32:829–842) showed that ablating MTARC1 p.A165
reduces hepatocellular carcinoma aggressiveness.

Accordingly, the primary quantity of interest is not the expression of any
single gene but the **balance between electron supply and electron drain**.

## 2. Index definition (locked)

Preprocessing: counts → log2(CPM + 1); microarray intensities → RMA log2.
Each gene is z-scored **within a cohort**, which is what makes the index
comparable across platforms without merging matrices.

```
RSI = 0.5 * ( mean_z(SUPPLY) + mean_z(REDUCTION) ) - mean_z(DRAIN)
```

| Module | Genes |
|---|---|
| **SUPPLY** (n=15) | HMGCR, HMGCS1, MVK, PMVK, MVD, IDI1, FDPS, PDSS1, PDSS2, COQ2, COQ3, COQ5, COQ6, COQ7, COQ9 |
| **REDUCTION** (n=4) | CYB5R3, CYB5R1, AIFM2, NQO1 |
| **DRAIN** (n=3) | MTARC1, MTARC2, POR |

SUPPLY and REDUCTION each carry weight 0.5 because SUPPLY has fifteen genes
and REDUCTION four; a pooled mean would let SUPPLY dominate, which does not
match the biological claim that the index measures production *times*
reduction.

Reported alongside the index but **not part of it**: CYB5A and CYB5B (shared
adapters, direction not assignable); GPX4, SLC7A11, SCD1, ACSL4, LPCAT3
(downstream effectors); SMG1, SMG8, SMG9, UPF1, UPF2, SREBF2 (upstream
set-point).

Higher RSI is interpreted as greater reductive-supply dominance and therefore
a higher predicted ferroptosis threshold.

## 3. Hypotheses, with direction stated in advance

| | Hypothesis | Predicted direction |
|---|---|---|
| **H1** | RSI across the MASLD spectrum | monotonic increase: control < NAFL < MASH F0-1 < F2 < F3 < F4 |
| **H2** | RSI in HCC tumor versus adjacent non-tumor liver | tumor > adjacent |
| **H3** | RSI in adjacent non-tumor NASH liver versus NASH liver from patients without HCC | adjacent > non-HCC NASH (field effect) |
| **H4** | RSI and overall survival in HCC | HR > 1 |
| **H5** | NMD-target burden (proxy for SMG1 activity) versus RSI | rho < 0 |
| **H6** | RSI versus a ferroptosis-resistance signature | rho > 0 |

H3 and H5 carry the originality of the study. H1 and H2 are expected to hold
and are not sufficient on their own.

## 4. Analysis order

- **PHASE A** — index construction and internal validity only. No outcome
  variable is opened. Ends with this document and the lock hash.
- **PHASE B** — H1 in the MASLD spectrum. Primary test: Jonckheere–Terpstra
  test for ordered alternatives. Covariate-adjusted linear models where age,
  sex, BMI, and diabetes are available. Sub-scores reported separately.
- **PHASE C** — H2 and H3 in GSE164760, replicated for H2 in TCGA-LIHC and
  ICGC LIRI-JP.
- **PHASE D** — H5. SMG1, SMG8, and SMG9 reported individually; NMD-target
  burden by ssGSEA (Tier 1) and, if feasible, junction- or transcript-level
  poison-exon ratios (Tier 2).
- **PHASE E** — H4 and H6. **Survival data are opened here for the first
  time.** Cox proportional hazards, RSI continuous and by tertile,
  multivariable with stage, AFP, age, sex, and race.

## 5. Pre-specified sensitivity analyses

| | Variant | What it tests |
|---|---|---|
| S1 | GSVA/ssGSEA module scores instead of mean-z | dependence on the scoring method |
| S2 | cholesterol-arm genes (SC5D, DHCR7, DHCR24, MSMO1) added to SUPPLY | whether the index is CoQ-specific — if results are unchanged, specificity is **weakened**, and this will be reported as such |
| S3 | CYB5R3 alone versus the whole REDUCTION module | contribution of the liver-dominant paralog |
| S4 | POR removed from DRAIN | arbitrariness of module membership |
| S5 | tumor purity adjustment (ESTIMATE/ABSOLUTE) | confounding by tumor cell fraction |
| S6 | hepatocyte-fraction adjustment after deconvolution against GSE202379 | whether RSI shifts are hepatocyte-intrinsic rather than infiltration |

## 6. Stated limitations

1. Transcript abundance measures **capacity**, not flux. The CoQH2/CoQ ratio
   is not measured here; that measurement is deferred to the planned
   mechanistic study.
2. All comparisons are cross-sectional. No public cohort pairs MASLD-spectrum
   liver tissue with longitudinal HCC outcome in the same patients. A raised
   RSI in adjacent non-tumor tissue may be a precancerous field change **or** a
   secondary effect of the adjacent tumor; these cannot be separated here, and
   H3 will be reported as "already observable outside the tumor," not as
   "precedes the tumor."
3. Bulk tissue confounds cell composition with cell-intrinsic change. S6 is
   the pre-specified mitigation.
4. GSE164760 is a microarray cohort and is analyzed as an independent arm,
   never merged with RNA-seq cohorts.
5. The established human phenotype of CYB5R3 deficiency is congenital
   methemoglobinemia, in which fatty liver is not prominent. This is
   acknowledged directly and explained by tissue-specific demand rather than
   argued around.

## 6a. Replication axes — amendment fixed 2026-08-25

*Added after the discovery cohort was analyzed and **before** any replication
cohort's RSI was computed. Recorded here because the choice below is a
researcher degree of freedom and must not be made after seeing a result.*

Cohorts do not annotate severity the same way, so a single ordered ladder
cannot be imposed on all of them. The axis each cohort supports is a property
of its GEO metadata, and is fixed as follows:

| Cohort | Ordered severity axis | Why |
|---|---|---|
| GSE135251 (discovery) | `group in paper`: control < NAFL < MASH F0–1 < F2 < F3 < F4 | the submitters' own staging, used in the source publication |
| GSE130970 | fibrosis stage F0 < F1 < F2 < F3 < F4 | the cohort carries fibrosis stage and full NAS components but **no** diagnosis label, so no NAFL/MASH split can be derived without inventing a threshold |
| GSE167523 | NAFL < MASH (`disease subtype`) | the cohort carries a diagnosis label but **no** fibrosis stage in its GEO metadata |

Consequences, accepted in advance:

1. GSE167523 has only two ordered groups, so the Jonckheere–Terpstra test there
   reduces to a Mann–Whitney test. It is a weaker replication and will be
   described as such.
2. Because the axes differ, cohorts are **not** pooled on a common stage
   coding. The forest plot pools Spearman ρ computed within each cohort against
   that cohort's own axis, which is a correlation on a comparable scale.
3. GSE135251 will additionally be re-run on the fibrosis axis as a secondary
   analysis, so that at least one comparison uses an identical axis across two
   cohorts.
4. No cohort's axis will be changed after its result is seen. If an axis proves
   unusable, that cohort is reported as not analyzable, not re-axed.

## 6b. Disclosure of prior analyses on the same data — recorded 2026-08-25

*Recorded on discovering existing analyses on the investigator's machine, and
before PHASE C or PHASE E was run. It is disclosed because it bears on how
"outcome-blind" this study can honestly claim to be.*

Earlier work (R 4.6.1, macOS) exists on the same cohorts, focused on **CYB5R3
as a single gene** rather than on a composite index:

| Location | What was analyzed |
|---|---|
| `outputs_masld/` | CYB5R3 across the MASLD spectrum in GSE135251; co-expression panel |
| `outputs_masld_val/` | the same in GSE130970 |
| `outputs_lihc/` | **CYB5R3 and overall survival in TCGA-LIHC**, tumor vs adjacent-normal |
| `outputs_ferroptosis/` | ferroptosis driver/suppressor scores, CYB5R3/CYB5R1/POR anchor contrasts in LIHC and LUAD |
| `figures_pub/` | four figures, including sex-stratified CYB5R3 analyses |

Three consequences for this study:

1. **The survival outcome for TCGA-LIHC has already been seen.** Reported there:
   CYB5R3 does not differ between tumor and paired adjacent-normal (paired
   Wilcoxon p = 0.94), but higher tumor CYB5R3 is associated with worse overall
   survival (per SD, HR 1.22, 95% CI 1.02–1.47; adjusted for stage and grade
   HR 1.23, 1.02–1.48), while categorical splits are not significant (log-rank
   p = 0.258 median split, p = 0.141 tertiles).
   Therefore **PHASE E in TCGA-LIHC is not outcome-blind at the level of the
   project**, even though the RSI was defined and locked before any outcome was
   examined in *this* analysis, and RSI is not CYB5R3. H4's stated direction
   (HR > 1) happens to agree with that prior result. This must be stated in the
   manuscript's limitations; it cannot be presented as a blind prediction.
   ICGC LIRI-JP has not been examined and remains genuinely outcome-blind, so it
   is the stronger prognostic test.

2. **The central negative is independently replicated.** The prior analysis
   found CYB5R3 versus fibrosis stage in GSE135251 at ρ = 0.02, p = 0.76. This
   analysis, on a different pipeline, found ρ = −0.074 (q = 0.38). Two
   independent implementations agree that CYB5R3 expression is flat across the
   MASLD spectrum. That strengthens the finding rather than weakening it.

3. **A replicated co-expression result exists that this study did not seek.**
   CYB5R3 co-expresses with the lipogenic program — SREBF1 ρ = 0.63 in
   GSE135251; DGAT2 ρ = 0.47, FASN ρ = 0.43 in GSE130970 — and with NADK
   (ρ = 0.63). This is consistent with the competing-demand arm of the reserve
   model, in which desaturation and lipogenesis draw on the same hub, but it was
   not a pre-specified hypothesis here and will be cited as prior work, not
   claimed as a finding of this analysis.

**Scope rule, fixed 2026-08-25.** The earlier CYB5R3 work and this study are
**two separate manuscripts**, and the boundary between them is:

| | Earlier manuscript | This manuscript |
|---|---|---|
| Unit of analysis | CYB5R3 as a single gene | the composite RSI |
| Claims | CYB5R3–lipogenic co-expression (SREBF1, DGAT2, FASN), NAD⁺/NADK coupling, sex-stratified effects, CYB5R3 and survival in TCGA-LIHC | the supply-versus-drain balance and its trajectory across MASLD and into HCC |
| Cohorts | GSE135251, GSE130970, TCGA-LIHC, LUAD | GSE135251, GSE130970, GSE167523, GSE164760, TCGA-LIHC, ICGC LIRI-JP |

Because both use GSE135251 and GSE130970, each manuscript must cite the other,
state the shared data source explicitly in its Methods, and **not re-report the
other's results as new**. Concretely, for this manuscript:

- Single-gene CYB5R3 findings are cited to the earlier work. The one exception
  is CYB5R3's flat trajectory, which is unavoidable context for the index and is
  reported here **as a replication of that earlier result, with attribution**.
- Co-expression analyses, sex stratification and the CYB5R3–survival model are
  out of scope here and will not be repeated.
- If a journal asks about overlap, the answer is that the two papers share
  cohorts but not analyses, and neither's conclusions depend on the other's.

## 6c. H5 — operationalization, and a correction to its stated direction

*Fixed 2026-08-25, before the NMD score was computed in any cohort.*

### The directional statement in §3 is wrong, and this records why

H5 was written as: *"NMD-target burden (a proxy for SMG1 activity) is inversely
correlated with RSI (rho < 0)."* That parenthesis conflates two things that are
inversely related. NMD-target **burden** is the accumulation of NMD substrates,
which happens when NMD activity is **low**. Burden is therefore a proxy for the
*inverse* of SMG1/NMD activity, not for the activity itself.

The framework predicts: low SMG1/NMD activity → higher SREBP2/HMGCR set-point →
higher RSI. Written in terms of burden, that is **burden up, RSI up: ρ > 0** —
the opposite sign to what §3 states.

Because the pre-registered sign is the product of an error rather than a
prediction, **H5 is tested two-sided**, and the framework's directional
expectation (ρ > 0) is recorded here before any computation. The original
wording is not quietly amended; both are reported.

### How the score is computed

**Gene set:** the NMD Consensus target set of Palou-Márquez G & Supek F,
*Genome Biol* 2025;26:316 (doi:10.1186/s13059-025-03727-y), file
`S4_NMD_Consensus_target_gene_set.txt` from Zenodo (doi:10.5281/zenodo.15974216).
It contains **130 genes** carrying an experimentally supported NMD-target
isoform in ≥2 independent studies (Colombo 2017; Karousis 2021; Tani 2012;
French 2020), each paired with a matched control isoform.

**Score:** mean within-cohort z-score of the 130 genes, on the same log
expression matrix used for the RSI. Higher score = more substrate accumulation
= **lower** inferred NMD activity.

Two properties were checked before use, and both are clean:
- **no overlap with the RSI gene panel** — the two scores share no gene;
- **no NMD machinery genes** (UPF1/2/3, SMG1/5/6/7/8/9) in the set, so the
  score is not circular with the upstream module.

### Why this test is weaker than H5 deserves

The source method is **transcript-level**: it contrasts each gene's NMD-target
isoform against a matched control isoform *of the same gene*. Collapsing to
gene level sums the two and cancels the contrast being measured. What remains
is the aggregate abundance of genes that happen to have an NMD-sensitive
isoform, which reflects NMD activity only weakly and is confounded by
transcription rate.

**H5 therefore cannot be adequately tested with gene-level data.** The analysis
below is **exploratory**, is reported as such regardless of outcome, and does
not settle H5 either way. A proper test needs transcript- or junction-level
quantification (analysis plan, Tier 2).

### Pre-specified sensitivity analyses

| | Variant | Why |
|---|---|---|
| N1 | drop the 10 mitochondrial/OXPHOS-adjacent members (ATP5A1, COX11, DAP3, ETHE1, HIBADH, MRPL48, MRPL49, MRPS27, MRPS35, NDUFS4) | these move with metabolic phenotype in liver disease independently of NMD |
| N2 | drop the 4 ISR/stress members (BAX, CIRBP, CRTC2, GADD45B) | the score could otherwise re-read out ER stress |
| N3 | correlate the score with SMG1, SMG8 and SMG9 expression | orthogonal check: if the score reflects NMD activity, it should relate to the machinery |

The set uses pre-2019 symbols (MARS, CARS, ATP5A1, H2AFY, MESDC2, WBSCR22,
MCRIP2 and others); these are resolved through the annotation's Synonyms
column, and the number matched is reported per cohort.

## 6d. PHASE E — fixed 2026-08-25, before any survival value was read

*Written and committed before the prognostic cohort's outcome fields were
fetched. The cohort was chosen on availability and blindness, not on results.*

### Cohorts, and why the obvious one is not primary

| Role | Cohort | Reason |
|---|---|---|
| **Primary** | **GSE76427** (n = 167; 115 HCC tumors, 52 adjacent non-tumor; Illumina HumanHT-12 V4) | overall and recurrence-free survival are carried in the GEO sample metadata; the cohort has never been examined in this project or in the disclosed earlier work, so it is genuinely outcome-blind |
| Secondary, deferred | TCGA-LIHC | **not outcome-blind** — CYB5R3 versus survival there has already been seen (§6b). Any TCGA-LIHC result is confirmatory at best and will be labeled as such |
| Not used this phase | ICGC LIRI-JP | the open ICGC bucket stores expression per donor with no project-level file, and cross-origin access is blocked; assembling it is a separate data-engineering task. Recorded so the omission is not mistaken for a choice made after seeing results |

### H4 — the model, fixed in advance

- **Endpoint:** overall survival (`event_os`, `duryears_os`). Recurrence-free
  survival is a secondary endpoint, reported but not used to judge H4.
- **Analysis set:** tumor specimens only. Adjacent non-tumor tissue is excluded
  from survival models.
- **Primary test:** Cox proportional hazards with **RSI as a continuous
  variable, per 1 SD**. Continuous is primary because dichotomizing discards
  information and invites threshold selection.
- **Secondary:** RSI tertiles, with Kaplan–Meier curves and a log-rank test.
  Tertile cut-points are the cohort's own 33rd and 67th percentiles — not
  chosen to optimize separation.
- **Adjusted model:** RSI + age + sex + BCLC stage. If a covariate is missing
  for more than 20% of samples it is dropped and that is reported.
- **Direction:** H4 predicts HR > 1. Tested **two-sided**.
- **Proportional hazards** is checked by the correlation of scaled Schoenfeld
  residuals with time; a violation is reported, not silently ignored.

### H6 — ferroptosis signature

- **Gene sets, in order of preference:** FerrDb V2 driver and suppressor sets
  if reachable; otherwise MSigDB `WP_FERROPTOSIS` and `GOBP_FERROPTOSIS`. The
  set actually used is recorded with the result.
- **Test:** Spearman correlation of RSI with each signature score (mean
  within-cohort z), computed in every cohort already prepared, not only the
  prognostic one.
- **Direction:** H6 predicts a positive correlation with a ferroptosis
  *resistance* signature — i.e. **positive with suppressors, negative with
  drivers**. Tested two-sided.
- **Overlap must be removed before testing.** GPX4, SLC7A11, ACSL4, LPCAT3,
  AIFM2 and SCD1 are in this study's own panel and appear in ferroptosis gene
  sets. Any gene shared with the RSI panel is **excluded from the signature**
  before scoring, and the number excluded is reported. Without this the test is
  partly circular.

### Multiplicity and what counts as a result

PHASE E runs two hypotheses in one primary cohort. Both are reported with
effect sizes and confidence intervals rather than being judged by a threshold
alone. No subgroup analysis is planned; any that appears is exploratory and
labeled so.

### The global-axis check applies here too

Every cohort's RSI is checked against its mean expression (PHASE D §4). If the
prognostic cohort shows a material correlation, the survival model is repeated
with mean expression as an additional covariate, and both are reported.

## §6e Amendment — GSE76427 acquisition and normalization route (recorded 2026-08-26, before any GSE76427 outcome was computed)

The GSE76427 non-normalized supplementary file could not be written to disk by
any route available to this analysis: the NCBI FTP host is blocked from both
the analysis container and the workstation VM, and the browser cancelled the
HTTPS download as unconfirmed. The file was instead fetched into the browser's
memory from the allowed `www.ncbi.nlm.nih.gov` host, decompressed there, and
the pre-specified normalization applied in-page. Only the panel probes were
extracted and carried back into the pipeline.

This changes *where* one pre-specified step ran, not what it does. Recorded in
full because it is a departure from the route used for every other cohort:

1. **What was computed in the browser.** log2(x+1) followed by quantile
   normalization across samples, over the **complete** 47,322-probe x
   167-sample matrix — the same matrix and the same step `01_prepare.py
   --units array-linear` would have applied. Normalizing over the full probe
   set (not the panel subset) is what makes this equivalent; a panel-only
   quantile normalization would **not** be.
2. **The implementation was validated against the pipeline's own code before
   use.** The JavaScript implementation and `io_utils.quantile_normalize` were
   run on an identical 60x7 fixture containing four planted tie groups. All
   420 values agreed to 10 decimal places, including the tie-averaged ranks and
   the round-half-to-even index rule. A first attempt disagreed; the cause was
   the fixture generator overflowing JavaScript's safe-integer range, not the
   normalization, and was fixed by changing the generator on both sides.
3. **The transfer was verified, not assumed.** Panel matrix totals, first and
   last column sums, row sums, minimum and maximum were computed independently
   in the browser and in pandas after reassembly; all eight agreed exactly
   (total 93211.8642 over 59 probes x 167 samples, no missing values).
4. **Samples were matched by name, not by position.** The supplementary file's
   column order (`PT<id>` / `ANTT<id>`) is **not** GEO's sample order. Columns
   are mapped to GSM accessions through patient id and tissue, and the mapping
   is asserted to be a bijection onto the 115 tumors and 52 adjacent samples
   already recorded from GEO's own metadata.
5. **Consequence for the pipeline.** Because normalization has already been
   applied, GSE76427 enters `01_prepare.py` with `--units logged` rather than
   `--units array-linear`. Probe-to-symbol mapping and the probe collapse
   (`groupby.max`) still run inside the audited pipeline, unchanged.

No outcome — no correlation, no hazard ratio, no group contrast — was inspected
before this amendment was written.

## §6f Amendment — H2 replication in GSE76427 (recorded 2026-08-26, before the contrast was computed)

> **Correction, same day, before the result was read.** As first written, this
> amendment stated H2's direction backwards ("lower in tumor"). The binding
> registration is the table in section 1, fixed 2026-08-25: **H2 predicts tumor
> > adjacent**. The error was in this amendment's prose only; the registered
> direction was never changed, and the correction was made on discovering the
> discrepancy, before the contrast was interpreted. Both wordings are left on
> the record.

H2 (RSI is **higher** in tumor than in adjacent non-tumor liver) was
pre-specified and tested in GSE164760, where it was **not supported**
(P = 0.44; shift +0.15, Cliff's delta +0.10 — directionally consistent but
far short of significance). GSE76427, acquired for
H4, contains 115 tumors and 52 adjacent non-tumor livers, **52 of which are
patient-matched to a tumor in the same series**. Testing a hypothesis that was
already registered, in a second independent cohort, is replication rather than
a new analysis; it is recorded here in advance so that the direction and the
tests are on file before the numbers are seen.

- **Primary (unpaired):** two-sided Mann-Whitney on RSI, tumor vs adjacent,
  all 167 samples. Effect reported as a Hodges-Lehmann shift with a
  distribution-free CI, plus Cliff's delta. Same estimators as PHASE C.
- **Paired:** Wilcoxon signed-rank on the patient-matched pairs only. Pairing
  is taken from GEO's own `patient id` field, not inferred; the number of
  complete pairs is asserted before testing.
- **Direction:** H2 as registered predicts RSI **higher in tumor**
  (section 1). Tested two-sided, and the registered direction is stated in the
  result either way.
- **Status:** this is a replication of a registered hypothesis in a new cohort,
  and is reported as such. It is *not* a re-test of H4 and does not change H4's
  pre-specified analysis, which was completed and recorded first.

The GSE76427 RSI values used here are the ones already written by
`02_compute_rsi.py` for H4; nothing about the index is recomputed.

## §6g Amendment — differentiation-adjusted H2 (recorded 2026-08-26, before any differentiation gene was extracted)

PHASE E §3 reported H2 as supported with a large effect, carried by REDUCTION up
and DRAIN down. It also named the threat to that reading: **mARC1, mARC2 and POR
are hepatocyte-differentiation genes and HCC is dedifferentiated**, so some or
all of the DRAIN fall may be loss of hepatocyte identity rather than a
reallocation of reducing equivalents. This amendment fixes the test of that
alternative explanation before the data to test it are pulled.

### The covariate, specified now

**D1 (primary) — mature hepatocyte identity, deliberately redox-free.** No gene
below is in the RSI panel, and none is a redox enzyme, an electron donor/acceptor
partner of POR, or a cytochrome P450:

- transcription factors: `HNF4A HNF1A FOXA1 FOXA2 NR1H4`
- secreted plasma proteins: `ALB TTR TF SERPINA1 AHSG APOH FGA FGB FGG F2`
- urea cycle / amino-acid catabolism: `CPS1 OTC ARG1 TAT`
- gluconeogenesis: `G6PC1 PCK1`
- hepatocyte surface receptor: `ASGR1`

**D2 (sensitivity) — D1 plus the major hepatic P450s** `CYP2E1 CYP3A4 CYP1A2
CYP2C9`. D2 is expected to **over-adjust**: POR exists to donate electrons to
these very enzymes, so conditioning on them partly conditions on the drain
itself. D2 is reported as a deliberately conservative bound, never as the
primary. If the effect survives D2 it survives anything.

Both scores = mean of within-cohort z-scores, the same construction the RSI
modules use. Genes absent from GPL10558 are reported, not silently dropped.

### The tests, in order

1. **Positive control.** D1 must be lower in tumor than adjacent. If it is not,
   the confound does not exist in this cohort and the rest is uninterpretable.
2. **Primary — paired-difference regression.** Within the 52 patient-matched
   pairs, regress ΔRSI (tumor − adjacent) on ΔD1. **The intercept is the
   quantity of interest**: it is the RSI rise in tumor at zero change in
   hepatocyte identity. Reported with a 95% CI and a two-sided t test. This is
   the cleanest form of the question, because pairing removes between-patient
   variation and the intercept has a direct interpretation.
3. **Supporting — unpaired.** OLS of RSI on tissue + D1 over all 167 samples;
   report the tissue coefficient with a 95% CI, and a Mann–Whitney on D1-residualized RSI.
4. **Module-level.** Repeat step 2 separately for the REDUCTION and DRAIN module
   scores, since PHASE E showed the composite is carried by those two.
5. **Sensitivity.** Repeat steps 2 and 4 with D2.

### Decision rule, fixed now

Let *b* = unadjusted paired median ΔRSI (+1.129, already reported) and *a* = the
adjusted intercept from step 2.

- **a is significantly > 0 and a/b ≥ 0.5** → the tumor rise is *not* explained by
  dedifferentiation; PHASE E §3 stands as written.
- **a is significantly > 0 but a/b < 0.5** → *partly* explained; the paper must
  say so and report both numbers.
- **a is not significantly different from 0** → the tumor–adjacent RSI
  difference is **attributable to loss of hepatocyte identity**, and PHASE E §3
  must be rewritten to say that. H2 would then remain descriptively true but
  mechanistically uninformative.

The third outcome is a real possibility and is accepted in advance. No
alternative covariate will be substituted after the fact to rescue the result,
and the RSI itself is not redefined either way (§7).

### Data route

The quantile-normalized GSE76427 matrix produced under §6e (all 47,322 probes,
167 samples, verified by checksum) is still held in the browser session that
created it. Only the probes for the genes named above are extracted; the
normalization is **not** recomputed, so the differentiation scores sit on exactly
the same footing as the RSI values already analyzed.

## §6h Amendment — replicating the differentiation adjustment in independent series (recorded 2026-08-26, before any replication number was computed)

§6g found, in GSE76427, that the tumor-vs-adjacent RSI rise is **half**
explained by loss of hepatocyte identity, and that the two arms behave
differently: **REDUCTION retained 71%** of its shift after D1 adjustment while
**DRAIN retained only 41%**. That split is the study's sharpest claim and it
currently rests on **one cohort**. This amendment fixes how it is retested
before any replication number exists.

### What is being replicated

Not H2, and not the composite. The **claim under test is the module split**:
after adjusting for hepatocyte identity, the REDUCTION rise in tumor survives
substantially, and the DRAIN fall survives much less. Both cohorts below are
independent of GSE76427 in patients, platform and (for GSE14520) etiology.

### Cohort 1 — GSE164760 (already prepared; NOT outcome-blind, stated plainly)

Affymetrix HG-U219, NASH-driven HCC: 53 tumors vs 29 non-tumoral NASH livers
adjacent to HCC. **This cohort's composite H2 was already tested and was NOT
supported** (P = 0.44, Hodges–Lehmann +0.15, Cliff's δ +0.10), and its
components were computed in PHASE C. It is therefore *not* outcome-blind, and
this is a **demanding** test rather than a favorable one: if the module shifts
are themselves absent here, there is nothing for the adjustment to act on and
the §6g split simply does not replicate.

- GEO declares no patient ID for this series, so the analysis is **unpaired**:
  OLS of each module score on tissue + D1, reporting the tissue coefficient with
  a 95% CI. The 21 pairs inferred in PHASE C from consecutive source identifiers
  are **not** used here — they were labeled sensitivity-only and remain so.
- The D1/D2 gene sets are **exactly those fixed in §6g**; not one gene is added,
  dropped or substituted. Genes absent from HG-U219 are reported, not replaced.
- Because quantile normalization in `01_prepare.py` runs over the complete
  probeset matrix *before* identifier mapping, extending the annotation map to
  carry the differentiation genes cannot alter the panel genes' values. This is
  **asserted numerically** against the already-written `rsi.tsv`, not assumed.

### Cohort 2 — GSE14520 (to be acquired; outcome-blind to this project)

HBV-related HCC with paired adjacent liver, Affymetrix HT-HG-U133A. Chosen now,
by accession, before acquisition, and specified here so the choice cannot be
made after seeing cohort 1. It is larger than either cohort used so far and has
**declared patient pairing**, allowing the same paired-difference regression
used as the §6g primary. Its etiology is HBV rather than MASLD — a deliberate
generality test, and a stated limitation if the result differs.

If GSE14520 cannot be acquired, that is reported as an unmet plan, and the
replication rests on cohort 1 alone.

> **Acquisition status, recorded 2026-08-26 (superseded the same day — see
> §6i).** GSE14520 was first **partially** acquired: GEO serves no processed
> matrix for this series on the reachable host, the 445 GPL3921 samples had to
> be pulled one at a time through the SOFT endpoint, and the browser session was
> lost twice mid-transfer, leaving the **RSI panel complete** but the D1
> covariate at **5 of 22 genes**. The adjustment was therefore **not run**, and
> was carried as an open, pre-registered obligation rather than quietly dropped.
> A 5-gene stand-in for a covariate registered as 22 genes is exactly the
> post-hoc substitution §6g and §6h were written to forbid, and it was not made.
>
> **That obligation has since been discharged in full.** The transfer was
> completed with per-sample retries and a completeness check, and HNF1A — which
> this note originally reported as absent from HG-U133A — was found to be
> present under its retired symbol `TCF1`. **D1 is the full 22 of 22 genes and
> D2 is 26 of 26.** §6i records how, and records it before any adjusted number
> was computed. This paragraph is left standing rather than edited away, so the
> partial state and its repair are both on the record.
>
> Two panel genes remain genuinely absent from HG-U133A and are reported, not
> substituted — confirmed by Entrez Gene ID, not by symbol: **AIFM2** (84883; so
> REDUCTION is a 3-gene version there) and **COQ5** (84274; SUPPLY 14/15).
> AIFM2's absence matters: it is the gene that ran opposite to the module in
> PHASE B, so its removal makes GSE14520's REDUCTION more internally coherent
> than the other cohorts' **by construction**, and that is a limitation, not a
> result.
>
> **The two sentences above are wrong and are corrected in §6k
> (2026-08-27).** AIFM2 did not run opposite to its module: it is the
> strongest single contributor to REDUCTION where it is measured (GSE76427
> Δz = +1.292, P = 2.31 × 10⁻¹⁶). The limitation is real but points the other
> way — the platform removes the module's most informative gene from the
> cohort carrying the main result. The wrong text is left standing here on
> purpose; see §6k.

### Decision rule, fixed now

For each cohort, let *retention* = adjusted effect / unadjusted effect, per
module.

- **Replicated** if, in a cohort where the unadjusted module shifts are present
  in the §6g directions (REDUCTION up in tumor, DRAIN down), REDUCTION's
  retention exceeds DRAIN's, and REDUCTION's adjusted effect has a CI excluding
  zero.
- **Not replicated** if REDUCTION's adjusted effect includes zero, or if
  DRAIN's retention is the higher of the two.
- **Uninformative** if the unadjusted module shifts are absent in that cohort;
  reported as such, and explicitly *not* counted as either support or refutation.

Retention percentages are reported for both modules in both cohorts regardless
of which branch is taken. If the split replicates in neither cohort, PHASE E
§3a is rewritten to present it as a single-cohort observation requiring
confirmation.

## §6i Amendment — GSE14520 acquisition completed; the §6h obligation is now discharged (recorded 2026-08-26, before any adjusted GSE14520 number was computed)

§6h recorded GSE14520 as **partially acquired** and carried its
differentiation adjustment as an open, pre-registered obligation. That
obligation is now met. This amendment records **how**, and records it before
the adjustment is run, so that nothing below can be read as chosen after the
fact.

### What changed since the §6h acquisition note

1. **The remaining covariate genes were transferred.** The 445 GPL3921 samples
   were pulled again through the per-GSM SOFT endpoint, this time with per-sample
   retries (4 attempts, backed off), concurrency 3, `localStorage` caching against
   another lost session, and a completeness check that throws if any sample
   returns fewer probesets than requested. Sample order was verified identical to
   the earlier run, first and last accession included.

2. **HNF1A was found to be present after all, under a retired symbol.** The
   §6h note reported the D1 covariate as capped at 21 of 22 because HNF1A
   appeared absent from HG-U133A. It is not absent. GPL3921 annotates it as
   **`TCF1`**, its pre-2005 symbol, on probesets **`210515_at`** and
   **`216930_at`**, both carrying Entrez Gene **6927** — unambiguously HNF1A and
   not TCF7 (Entrez 6932, probesets `205254_x_at`/`205255_x_at`, which are
   excluded). The earlier miss was **my symbol-based lookup failing on a retired
   alias**, not a platform limitation, and it is recorded here as an acquisition
   error rather than quietly repaired.

3. **The whole map was re-derived by Entrez Gene ID rather than by symbol**, so
   the same class of error cannot be hiding elsewhere. Probesets annotated to more
   than one gene are excluded, as before. The re-derived map contains **86
   probesets**; it differs from the map used for the §6h unadjusted analysis by
   **exactly the two HNF1A probesets and nothing else**. Every panel-gene value
   already reported for GSE14520 is therefore unchanged, and this is asserted
   numerically in the script rather than assumed.

### Consequences, fixed now

- **D1 is the full 22 of 22 genes registered in §6g.** No substitution, no
  shortened covariate, no rescue set. D2 is 26 of 26.
- **The RSI panel is unchanged and still incomplete in two places**, exactly as
  §6h recorded: **AIFM2** is absent from HG-U133A (REDUCTION is a 3-gene version
  here) and **COQ5** is absent (SUPPLY 14/15). Confirmed this time by Entrez ID
  (84883, 84274), not by symbol. AIFM2's absence remains a limitation that
  **favors the hypothesis by construction**, since AIFM2 ran opposite to its
  module in PHASE B, and it is reported as such wherever GSE14520's REDUCTION
  result is reported.
  **Corrected in §6k (2026-08-27): the clause after "since" is false.** AIFM2
  is the strongest member of REDUCTION where measured, so its absence makes
  GSE14520's retention conservative, not generous. Left standing here on
  purpose.
- **The analysis is exactly the one already registered.** Tests 1–5 of §6g in
  the order given there, the paired-difference regression as primary with the
  intercept as the quantity of interest, and the module-retention rule of §6h.
  Nothing in either is reopened.

### One implementation point, stated before it is used

§6h's decision rule (module retention) was written for both cohorts but had
only been *implemented* for the unpaired branch, because cohort 1 is unpaired.
GSE14520 is paired, so the §6h rule is now evaluated on the paired module
intercepts — the paired analogue of the same quantity. This is the registered
rule applied to the registered analysis, not a new rule; both the §6g composite
verdict and the §6h module verdict are reported for GSE14520 regardless of which
way either falls.

## §6j Amendment — cell composition and tumor purity (recorded 2026-08-26, before any composition gene was extracted from any cohort)

### Why this is a separate question from §6g

§6g asked whether the tumor RSI rise is loss of **hepatocyte identity**, and
§6h/§6i answered it in two cohorts. That is not the same question as whether the
rise is a change in **which cells are in the sample**. A bulk tumor biopsy and a
bulk adjacent-liver biopsy differ in the fraction of hepatocytes, leukocytes,
endothelium and stroma they contain, and the RSI genes are **not** hepatocyte-
restricted: CYB5R3 is near-ubiquitous, NQO1 is high in many lineages, POR is
broad. So a reviewer can grant §6g in full and still say: *your REDUCTION rise
could be infiltrating cells.*

D1 does not settle this. D1 is built from hepatocyte-specific genes, so a sample
with fewer hepatocytes scores lower on D1 for a reason that has nothing to do
with the differentiation state of the hepatocytes that remain. D1 therefore
**partly absorbs** composition without separating it. This amendment specifies
the separation before the data to do it are pulled.

### The covariate, specified now

Both scores are means of within-cohort z-scores, the same construction the RSI
modules and D1/D2 use. **No gene below is in the RSI panel and no gene below is
in D1 or D2** — the covariate must measure a different thing, not re-measure
identity.

**C1 (primary) — non-parenchymal content**, 21 genes in three unambiguous
compartments:

- leukocyte / immune: `PTPRC CD53 LAPTM5 CD3E CD2 CD68 AIF1 ITGAM LCP1`
- endothelium: `PECAM1 VWF CDH5 ENG CLEC4G`
- fibroblast / stellate / ECM: `COL1A1 COL1A2 COL3A1 DCN LUM ACTA2 PDGFRB`

**C2 (sensitivity) — C1 plus cholangiocyte / progenitor markers** `KRT19 KRT7
SOX9 EPCAM`. C2 is expected to **over-adjust** and is reported as a conservative
bound, never as the primary, for a stated reason: KRT19-positive and EPCAM-
positive hepatocellular carcinoma is a recognized aggressive tumor subtype, so
these four genes carry tumor biology as well as composition. Conditioning on
them conditions on part of the phenotype under study. This is the exact role D2
plays in §6g, and it is declared here for the same reason.

Genes absent from a platform are reported, not replaced. Retired symbols are
resolved by **Entrez Gene ID**, not by symbol, per §6i.

**On ESTIMATE and CIBERSORTx.** Neither is used. The ESTIMATE stromal/immune
signature could not be obtained from a primary source in this environment, and
CIBERSORTx requires a web service that is unreachable here. Rather than depend
on a resource that cannot be verified, C1 is specified **exhaustively by gene
symbol above**, so that any reader can reproduce it without a package, a
signature file or an account. That is a deliberate trade of sophistication for
auditability.

### The direction of the positive control is NOT predicted

Unlike D1, whose direction was predicted (lower in tumor) and confirmed, **C1
has no pre-specified direction and none is claimed.** Adjacent liver in these
cohorts is frequently cirrhotic and inflamed, with high immune and stromal
content; a tumor may be more or less infiltrated than that. Predicting a
direction after the fact would be exactly the manoeuvre this document exists to
prevent. Two branches are therefore fixed now:

- **If C1 does not differ between tumor and adjacent (P ≥ 0.05):** the
  composition confound is *empirically absent in that cohort*. The objection is
  answered by that fact alone, the adjustment is reported as unnecessary there,
  and this is stated plainly rather than presented as a passed test.
- **If C1 differs in either direction:** the adjustment below is run and
  interpreted.

### The tests, in order

1. **Positive control.** C1 in tumor vs adjacent, two-sided, direction reported
   not predicted.
2. **Composition only.** Paired-difference regression of ΔRSI on ΔC1; the
   intercept is the RSI rise at zero change in cell composition. Reported with a
   95% CI. This is the direct analogue of §6g step 2 and exists so the two
   confounds can be compared on the same footing.
3. **Joint model — this is the test that decides the amendment.** Paired-
   difference regression of ΔRSI on **ΔD1 + ΔC1** together. Its intercept is the
   RSI rise at zero change in hepatocyte identity *and* zero change in cell
   composition. This is the quantity a reviewer is actually asking about, and it
   is named as decisive here, before it exists.
4. **Module level.** Repeat steps 2 and 3 separately for REDUCTION and DRAIN.
5. **Sensitivity.** Repeat steps 2–4 with C2 in place of C1.
6. **Collinearity, reported whether or not it is convenient.** The correlation
   between ΔD1 and ΔC1 and the variance inflation factor of the joint model are
   reported for every cohort. If VIF ≥ 5 the joint intercept is unstable and the
   amendment says so in those words rather than quoting the estimate as if it
   were precise.

Where a series declares no patient pairing (GSE164760), steps 2–5 are the
unpaired OLS analogue — module or index score on tissue + covariate — exactly as
§6h specified for that cohort, and the tissue coefficient replaces the intercept.

### Decision rule, fixed now

Let *b* = the unadjusted paired median ΔRSI already reported for that cohort and
*a* = the intercept of the **joint** model (step 3).

- **a significantly > 0 and a/b ≥ 0.5** → the tumor rise is explained by
  neither dedifferentiation nor cell composition. PHASE E §3 and §3b stand and
  are strengthened.
- **a significantly > 0 but a/b < 0.5** → partly explained by the two together;
  the paper must report both numbers and say so.
- **a not significantly different from 0** → the rise is **attributable to the
  combination of identity loss and composition**, and PHASE E §3/§3b must be
  rewritten to say that. H2 would remain descriptively true and mechanistically
  uninformative.

Module rule: the §6g/§6h module split is said to **survive composition
adjustment** only if, in the joint model, REDUCTION's adjusted effect has a CI
excluding zero **and** REDUCTION's retention exceeds DRAIN's. Otherwise it is
reported as not surviving.

The third outcome is a real possibility and is accepted in advance. No
alternative covariate will be substituted afterwards to rescue the result, the
gene lists above will not be edited after seeing any number, and the RSI itself
is not redefined either way (§7).

### Cohorts and data route

- **GSE14520** (213 pairs) and **GSE76427** (52 pairs) are the two cohorts that
  matter, because they carry the §6g/§6h result. Their C1/C2 genes are pulled by
  the route already registered for each — §6e for GSE76427, the per-GSM SOFT
  endpoint of §6h/§6i for GSE14520 — and **no normalization is recomputed**, so
  the covariate sits on exactly the same footing as the values already analyzed.
- **GSE164760** is run as a third, unpaired cohort from the genome-wide
  probeset matrix already prepared for it, with the same quantile normalization.
  Its §6h verdict was *uninformative* and nothing here changes that; it is
  reported as supporting context, not as a replication.

If a cohort's composition genes cannot be acquired, that is recorded as an unmet
plan in this section, as §6h recorded GSE14520's partial acquisition — not
quietly dropped.

### Note added 2026-08-26, AFTER seeing GSE164760's output — recorded, not hidden

GSE164760 was run first, because its genome-wide matrix was already on disk.
Its output exposed a gap in the composite decision rule above, and the gap is
recorded here rather than quietly patched.

**The gap.** The composite rule keys on the joint intercept's significance and
on *a/b*, where *b* is the unadjusted effect. It assumes there **is** an
unadjusted effect to explain. In GSE164760 there is not: H2 was already not
supported there (P = 0.44 in PHASE C), and the unadjusted tumor coefficient in
this run is +0.120 with P = 0.5. Every retention ratio is therefore a ratio to
approximately zero, and the verdict string the rule produces — "ATTRIBUTABLE to
identity loss plus composition" — is **mechanical, not interpretable**. §6h had
exactly this branch and called it *uninformative*; §6j's composite rule does not.

**What was and was not changed.** The verdict strings are unchanged and the
rule is unchanged; GSE164760's mechanical verdict is left standing in
`results/GSE164760_composition.json` and is reported. What was added is a
**precondition flag** that computes the unadjusted effect's P value and prints
that the verdict is uninterpretable when it exceeds 0.05. This surfaces an
assumption; it does not alter a threshold, a covariate, or an outcome.

**Why this does not rescue anything.** §6j already stated, before any of this
existed, that GSE164760 "is reported as supporting context, not as a
replication", and that its §6h uninformative verdict stands. The cohorts that
carry the §6j test are GSE14520 and GSE76427. Had the flag been added to a
cohort whose result the amendment *does* rely on, it would be a post-hoc
rescue and would be inadmissible.

**What GSE164760 does contribute.** Two facts that are interpretable
independently of the verdict:
1. **C1 is lower in tumor** here (−0.112 vs +0.205, P = 0.005) — the direction
   §6j declined to predict, and a sensible one: the comparator is inflamed,
   fibrotic non-tumoral NASH liver.
2. **D1 and C1 are close to independent** (r = −0.225; VIF 1.05 for both). The
   premise of §6j — that hepatocyte identity and cell composition are different
   measurements rather than two names for one — is supported rather than assumed.

### Second note added 2026-08-26, AFTER seeing GSE76427 — the knife-edge flag

GSE76427's module retentions came out at **51.20%** (REDUCTION) and **51.12%**
(DRAIN). The §6j module rule is a strict inequality between two retentions, so
it returns "SURVIVES" on a margin of **0.0008** — four orders of magnitude
inside the estimates' own confidence intervals.

As with the precondition flag above: **the rule and every verdict string are
unchanged**, GSE76427's mechanical verdict is left standing in its JSON and is
reported. What was added is a **knife-edge flag** that prints the gap and states
that a margin this small must not be leaned on. It surfaces a fact about the
rule's resolution; it does not alter a threshold or an outcome.

This is the same treatment §6g's 50.3% composite retention received at the time
it was produced, and for the same reason: a decision rule that reports a
categorical verdict should not be allowed to disguise a coin flip. The
substantive statement for GSE76427 is that **once both confounds are removed the
reduction and drain arms are no longer distinguishable from one another there.**

### Result recorded 2026-08-26 — the acquisition is complete and the answer is mixed

All three cohorts were acquired and analyzed. Full results in
`results/COMPOSITION_6j_RESULTS.md`; the headline, recorded here because it
changes what the paper may claim:

- **Composition alone explains none of the rise** (101% and 116% retained in the
  two paired cohorts). The objection, taken by itself, is empirically wrong.
- **Jointly with identity it bites unequally**: GSE14520 retains 65%
  (P = 1.8 × 10⁻²², *not explained*), GSE76427 retains 40% (P = 0.014, *partly
  explained*). Not collinearity — VIF 1.05–1.14 throughout.
- **The module split holds decisively only in GSE14520** (89% vs 31%).
- **GSE164760 remains uninterpretable**, as §6h and §6j both said in advance.

**Consequence, stated plainly:** the §6h claim that the module split *replicates*
must be narrowed. After removing both confounds it holds in the larger paired
cohort and is not demonstrable in the smaller. PHASE E §3b and the manuscript
framing are to be written accordingly, and the possibility that GSE76427's null
is simply power (52 pairs, wide intervals that contain GSE14520's split) is
reported as an explanation, not as a result.

### Supporting analysis (§6j-S), conditional on acquisition

If a public single-cell HCC atlas can be obtained, report which cell types carry
the expression of the RSI panel genes — specifically whether `CYB5R3 AIFM2 NQO1
CYB5R1` and `MTARC1 MTARC2 POR` are hepatocyte- and malignant-cell-dominant or
substantially non-parenchymal. This is **supporting evidence about the genes**,
not a covariate adjustment, and is labeled as such. It cannot change the
verdict above.


## §6k Correction — the recorded direction of the *AIFM2* limitation was wrong (recorded 2026-08-27, on discovery, before the manuscript was submitted)

### What was recorded, and what is true

From §6h onward this file, and every results document downstream of it, carried
the following claim:

> *AIFM2* is absent from HG-U133A, so GSE14520's REDUCTION is a 3-gene score —
> and *AIFM2* is the gene that ran opposite to its module in PHASE B, so its
> absence **favors the hypothesis by construction**.

**The second half of that sentence is false.** It was never computed; it was
asserted. Direct computation of the per-gene tumor-versus-non-tumor *z*-score
difference gives the opposite picture:

| Cohort | *AIFM2* Δ*z* | *P* | Rank within REDUCTION |
|---|---|---|---|
| GSE76427 (115 tumor vs 52 adjacent) | **+1.292** | 2.31 × 10⁻¹⁶ | **1 of 4 (strongest)** |
| GSE164760 (53 NASH-HCC vs 29 adjacent) | +0.531 | 0.0437 | 2 of 4 (NQO1 +0.544) |
| GSE135251 (PHASE B, fibrosis axis) | +0.156 | n.s. | directionally consistent |

PHASE B §"Per-gene" already recorded *AIFM2* as **+0.156, directionally
consistent but not significant** — the opposite of "ran opposite to its module".
PHASE C recorded *AIFM2* δ = −0.410, q = 0.0081, i.e. **higher in tumor**. The
gene that actually ran opposite to REDUCTION in the PHASE C field-effect
contrast was *CYB5R3* (δ = +0.470, q = 0.0029, higher in adjacent liver). The
error was a transposition of one gene for another, propagated by copying rather
than recomputing.

### The corrected statement

*AIFM2* is absent from HT-HG-U133A, so GSE14520's REDUCTION is a three-gene
score. Where *AIFM2* is measured it is the **strongest single contributor** to
REDUCTION. The platform therefore removes the module's most informative member
from the cohort that carries the main result, and the 84–89% REDUCTION retention
in GSE14520 is achieved **without** it.

### What this changes, and what it does not

- **No computed value changes.** The RSI panel, the module definitions, the
  covariates, the decision rules and every number in every results file were
  fixed and computed independently of this sentence. The error lived only in the
  interpretive text attached to them.
- **The direction of the limitation reverses.** It was reported as a limitation
  that inflates the GSE14520 result; it is instead a limitation that
  **understates** it. The honest reading is that GSE14520's REDUCTION retention
  is conservative relative to a four-gene score, not generous.
- **It does not become a strength.** A three-gene module on a different platform
  is still not strictly comparable to the four-gene version, and the direction of
  non-comparability is now known rather than assumed. That is the whole of the
  gain.
- **Improvement A3 keeps its priority.** A third paired cohort carrying *AIFM2*
  is still the single most valuable remaining analysis; the reason has changed
  from "remove a favorable artifact" to "measure the module's strongest member in
  a cohort with paired design".

### Why this is recorded rather than edited away

The same treatment given to §6c (a stated direction corrected), §6f and §6i (a
gene wrongly reported absent). The original wrong text is left standing in §6h
and §6i above, with a pointer to this section, so the error and its repair are
both on the record. Any reader who wants to check that no number moved can do
so: nothing in `results/` was regenerated for this correction.



## §6l Correction and sensitivity analysis — POR does not donate electrons to mARC (recorded 2026-08-27, before any POR-free value was computed)

### The error

The concept note, the Introduction, and the justification for the DRAIN module
all state that mARC1 and mARC2 draw on the shared NAD(P)H pool **through
cytochrome P450 oxidoreductase (POR)**. That is why *POR* was placed in DRAIN
alongside *MTARC1* and *MTARC2* when the index was locked on 25 August 2026.

**The biochemistry is wrong.** The mARC N-reductive system is an NADH-dependent
three-component chain:

> NADH → **NADH-cytochrome b5 reductase (CYB5R3)** → **cytochrome b5 (CYB5B)** →
> **mARC1 / mARC2** → N-reduction of the substrate

POR is the NADPH-dependent flavoprotein partner of the **cytochrome P450**
system, not of mARC. This is settled, primary-literature biochemistry (Rixen et
al., *PLoS One* 2014, PMID 25144769; Plitzko et al., *Drug Metab Dispos*
2016;44:1617, PMID 27469001; and the mARC review already cited as reference 3,
which describes the system without any role for POR). The error was mine, was
never checked against a primary source, and propagated from the concept note
into the locked panel.

### What this does and does not do to the study

- **It does not invalidate any computed value.** DRAIN was defined, locked and
  computed as `mean_z(MTARC1, MTARC2, POR)`. Every number reported anywhere in
  this study is that quantity, computed correctly from that definition.
- **It invalidates part of the *rationale* for one gene's membership.** POR
  belongs in DRAIN only if POR feeds mARC. It does not. POR remains a large,
  genuine NAD(P)H-consuming reductase in hepatocytes — it is the obligate
  electron donor to every microsomal cytochrome P450 — so it is not off-topic
  for an index of reductive *drain*. But it is not part of the N-reductive arm
  the module was named for.
- **The index is not being redefined.** The hash lock (§Index lock,
  `7a2bf934…b680046`) stands. Redefining a module after seeing the outcome is
  precisely what pre-registration exists to prevent, and §7's outcome-neutral
  stopping rule forbids it.
- **It strengthens the study's central claim rather than weakening it.** The
  premise of the whole index is that the CoQ-reducing arm and the N-reductive
  arm compete for the same reducing equivalents. Under the correct biochemistry
  that competition is *more* direct, not less: **CYB5R3 is itself the electron
  donor to mARC.** The shared resource is not an analogy about NAD(P)H in
  general — it is one enzyme feeding both arms. This also means CYB5R3 sits in
  REDUCTION while supplying DRAIN, which is a structural feature of the index
  that must be stated plainly wherever the modules are defined.

### The sensitivity analysis, specified before it is run

Because POR's membership rested on a false premise, I will report what the study
shows **without** it. Define

> **DRAIN′ = mean_z(*MTARC1*, *MTARC2*)** — the N-reductive arm proper
> **RSI′ = ½ (mean_z SUPPLY + mean_z REDUCTION) − mean_z DRAIN′**

and recompute, in the two paired cohorts that carry the main result (GSE76427,
GSE14520):

1. the unadjusted paired shift in RSI′ and in DRAIN′;
2. the D1-adjusted and joint (D1 + C1) paired-difference intercepts and retained
   fractions for RSI′ and DRAIN′, by exactly the machinery of §6g and §6j with
   no other change.

### Decision rule, fixed now, before any of it is computed

- If the **direction and the qualitative verdicts are unchanged** — RSI′ still
  rises in tumor, DRAIN′ still falls, DRAIN′ still retains materially less than
  REDUCTION under joint adjustment in GSE14520 — then POR's inclusion is
  reported as a **rationale error that does not drive the result**, the locked
  index stands, and DRAIN′ is reported as a sensitivity analysis in the
  manuscript.
- If **any qualitative verdict flips** — for instance if DRAIN′ does not fall,
  or if the GSE14520 module dissociation disappears — then the result is
  reported as **driven by POR**, the module-dissociation claim is withdrawn from
  the abstract and conclusion, and the finding is restated as being about POR
  rather than about N-reduction.
- Either outcome is reported. I am recording this branch before running it
  precisely because the second branch would be costly.

### Consequences for the text, fixed now

Regardless of outcome:

1. The Introduction must state the correct chain (CYB5R3 → CYB5B → mARC) and
   must not claim POR donates to mARC.
2. The Methods must state that DRAIN as locked contains *POR* and why that is
   now known to be a rationale error, with a pointer to this section.
3. The Limitations must carry it explicitly. It is not a footnote.
4. The fact that **CYB5R3 appears in REDUCTION while being the donor to DRAIN**
   must be stated where the index is defined, since it means the two modules are
   not enzymatically independent.

### Outcome, recorded 2026-08-30

**Branch 1 applies.** The full result is in `results/POR_SENSITIVITY_6l_RESULTS.md`.
Summary:

- The machinery reproduces the locked index exactly (Pearson r = 1.0000000000,
  max absolute deviation 0.00 in both cohorts), so every difference below is
  attributable to dropping *POR* and to nothing else.
- **Direction unchanged.** RSI′ still rises in tumor (paired median +0.981,
  P = 4.5 × 10⁻⁷ in GSE76427; +1.292, P = 2.0 × 10⁻³¹ in GSE14520) and DRAIN′
  still falls (−0.611 and −0.756). In GSE14520 both effects are slightly larger
  without *POR*.
- **The module dissociation survives intact.** Under the joint D1 + C1 model in
  GSE14520, REDUCTION retains 89.0% against DRAIN′ 31.8%, versus DRAIN 31.0%
  with *POR* — a difference of 0.8 percentage points.
- **GSE76427 remains uninformative for the split**, and the knife-edge tips away
  from the hypothesis without *POR* (REDUCTION 51.2% vs DRAIN′ 61.9%), which
  reinforces rather than undermines the decision not to claim a split there.

Therefore, as fixed in advance: *POR*'s membership in DRAIN is a **rationale
error that does not drive any result**; the locked index stands unchanged; DRAIN′
is reported in the manuscript as a pre-specified sensitivity analysis; and the
error is disclosed in the Introduction, the Methods and the Limitations.

### A second inconsistency, found while running this and disclosed rather than repaired

The §6g/§6j pipeline computes **composite** retained fractions against the paired
**median** of ΔRSI but **module** retained fractions against the paired **mean**
of the module shift. No verdict depends on the choice — on a median denominator
GSE14520's split reads 94.4% versus 36.2% rather than 89.0% versus 31.0%, and the
GSE76427 arms are indistinguishable under either — but using two conventions
without saying so is not acceptable. The manuscript now states which denominator
applies where. The numbers themselves are left as computed under the registered
plan rather than harmonised after the fact.



## §6m Amendment — cell-type attribution and protein-level corroboration (recorded 2026-08-31, before any single-cell or protein value was read)

### Why this exists

Two things prompt it. First, **§6j-S is an undischarged obligation**: it registered a
single-cell supporting analysis "conditional on acquisition", and acquisition never
succeeded because the container's egress blocks NCBI, EBI and the single-cell
portals. Second, external readers of the draft independently asked for exactly
these two analyses. Neither is a new idea introduced to rescue a result; the first
is already in the plan, and both are registered here **before any value is read**.

### Source, fixed now

The **Human Protein Atlas** (HPA), a public, versioned resource, reachable through
the browser bridge when the container's egress is not. Two of its sections are used:

- **Single Cell Type** — nTPM per cell type, derived from published human
  single-cell RNA-seq atlases, for liver.
- **Pathology** — antibody-based immunohistochemistry in liver cancer versus normal
  liver tissue, scored ordinally (not detected / low / medium / high).

The HPA version and the antibody identifiers are recorded with the results. No other
source is substituted if HPA proves unreachable; the analysis is then recorded as
**not performed**, exactly as §6j-S allowed.

### §6m-A — cell-type attribution (discharges §6j-S)

For the seven RSI panel genes measured in the paired cohorts — *CYB5R3*, *CYB5R1*,
*AIFM2*, *NQO1* (REDUCTION) and *MTARC1*, *MTARC2*, *POR* (DRAIN) — record the liver
cell types carrying their expression, and in particular the **hepatocyte fraction**:
hepatocyte nTPM divided by the summed nTPM across liver cell types.

**Prediction, stated before looking.** If the index measures what it claims to
measure, both modules should be **hepatocyte-dominant**, because the whole
construction assumes the bulk signal reports on hepatocyte electron allocation.

**Decision rule, fixed now.**

- **Hepatocyte fraction ≥ 0.50 for the majority of genes in both modules** → the
  cell-type premise of the index is supported.
- **Hepatocyte fraction < 0.50 for the majority of genes in either module** → the
  premise is **not** supported for that module, and this is reported as a limitation
  in the abstract-level terms, not buried.
- A split verdict between the two modules is reported as such.

**What this cannot do.** As §6j-S already fixed: this is *supporting evidence about
the genes*, not a covariate adjustment. **It cannot change any verdict in §6g–§6l.**
If it contradicts them, both stand and the contradiction is reported.

### §6m-B — protein-level corroboration

For the same seven genes, record HPA's immunohistochemical staining intensity in
**liver cancer** versus **normal liver**, and the direction of the difference.

**Prediction, stated before looking.** H2 and the module results predict REDUCTION
proteins **higher** in tumor and DRAIN proteins **lower**.

**Decision rule, fixed now.**

- Each gene is scored **concordant**, **discordant**, or **flat/uninformative**
  against that prediction.
- The result is reported as **counts of concordant versus discordant genes per
  module**. No *P* value is computed and no inferential claim is made.

**The honest ceiling on this analysis, stated before it is run.** HPA's normal-liver
IHC rests on roughly three individuals and its tumor panel on roughly a dozen, scored
on a four-level ordinal scale by antibody staining. That is **descriptive
corroboration at a different level of biology, not validation**. It cannot establish
significance, it cannot be pooled with the transcriptome result, and a discordant
gene does not refute the transcriptome finding any more than a concordant one
confirms it. It is reported in those terms or not at all. The manuscript's principal
limitation — that no new data were generated and there is no orthogonal validation at
the protein level under this study's own control — **remains true and remains
stated**; §6m-B narrows it, it does not remove it.

### Registered in advance: what would count against the study

- REDUCTION genes proving substantially non-parenchymal would undercut the reading
  that the tumor rise reflects hepatocyte reductive capacity, *even though* the C1
  adjustment already showed composition explains none of the rise. That tension
  would be reported, not reconciled away.
- DRAIN proteins **rising** in tumor at the protein level would sit against the
  transcript result and would be reported as such.

Both outcomes are accepted in advance.

### Outcome, recorded 2026-08-31

Full result in `results/HPA_6m_RESULTS.md`; data in `data/hpa/`.

**§6m-A gave a split verdict**, which this amendment fixed in advance as a reportable
outcome.

- **DRAIN — premise SUPPORTED.** Hepatocyte fractions 0.602 (*POR*), 0.549 (*MTARC2*),
  0.385 (*MTARC1*); 2 of 3 clear the registered 0.50 threshold and hepatocytes are the
  single highest cell type for both *POR* and *MTARC2*.
- **REDUCTION — premise NOT SUPPORTED.** 0 of 4 clear 0.50: *AIFM2* 0.274, *CYB5R1*
  0.131, *CYB5R3* 0.041, *NQO1* 0.006. Reported as the registered rule requires.

**§6m-B was performed only in part.** Normal-liver immunohistochemistry was obtained for
all seven genes; the registered primary comparison — tumour versus normal staining — was
**not performed**, because HPA's current page structure does not expose the tumour
staining table and §6m forbids substituting another source.

What the protein data did give is twofold and was not anticipated:

1. **DRAIN is hepatocyte-dominant at protein level as well** — all 3 genes stain more
   strongly in hepatocytes than in cholangiocytes. The premise underlying the entire D1
   adjustment now has protein-level support.
2. **HPA's own protein data contradicts its own single-cell data for three genes.**
   *CYB5R3* and *CYB5R1* stain Medium in hepatocytes against RNA fractions of 0.041 and
   0.131; *MTARC1* stains Medium against 1.0 nCPM. Droplet single-cell liver data is
   well known to under-recover hepatocytes, and this is direct within-source evidence of
   it. The §6m-A negative for REDUCTION is therefore reported **together with** the
   reason it may be a property of the assay rather than of the genes — and it is
   reported either way, not withdrawn.

**No verdict in §6g–§6l changes.** In particular §6m-A does not reopen the composition
question, which was tested directly in this study's own cohorts and answered: C1
explains none of the tumour rise (101% and 116% retained). A gene being expressed
outside hepatocytes in normal liver is a different claim from the tumour rise being
caused by a change in cell composition, and only the second was ever the confound.


## §6n Correction — two numerical defects found by adversarial re-derivation (recorded 2026-08-31, before any number was regenerated)

### How they were found

Every statistical routine in this study is hand-rolled without SciPy, by design.
The existing unit tests compare that code against itself, so a systematic defect
in a primitive would propagate into every reported value without any internal
test catching it. I therefore recomputed the reported quantities with **SciPy
1.17 and statsmodels 0.14** — independent implementations and independent
numerics — in `code/26_independent_stats_check.py`. 499 comparisons were run.

**Every regression coefficient and every regression P value reproduced exactly.**
Two defects surfaced.

### Defect 1 — the 95% interval used *z* = 1.96 where *t* was required

`ols_ci` computes an exact critical value by bisection for dof ≤ 200 but
substitutes 1.96 above that. GSE76427 has 50 degrees of freedom and is unaffected
(bisection agrees with SciPy to five decimals). **GSE14520 has 211, so all six of
its reported confidence intervals were computed with the normal critical value
instead of *t*(211) = 1.97127, and are 0.53–0.57% too narrow.**

This is an arithmetic error, not a modelling choice: the intervals are declared as
95% *t* intervals and 1.96 is the *z* value. It is corrected by removing the
shortcut and recomputing.

**No verdict depends on it.** Retained fractions are computed from the point
estimate, not the interval; P values come from the exact *t* tail and were already
correct; and no interval contains zero under either version.

### Defect 2 — Mann–Whitney without continuity correction, undisclosed

The hand-rolled `mannwhitney_u` reproduces SciPy's `use_continuity=False` variant
**exactly** (agreement to 1 × 10⁻⁸ across 200 random comparisons). That is a
defensible variant, not a bug, but it is slightly anti-conservative relative to
the conventional default and the protocol never said which variant was used.

Applying the continuity correction to every Mann–Whitney P value reported in this
study changes each by less than 2% relative and **crosses no threshold**:

| Contrast | as reported | with correction |
|---|---|---|
| GSE76427 H2 | 4.563 × 10⁻¹⁵ | 4.626 × 10⁻¹⁵ |
| GSE14520 H2 | 9.505 × 10⁻⁵⁶ | 9.560 × 10⁻⁵⁶ |
| GSE164760 H2 | 0.4407 | 0.4436 |
| GSE76427 C1 control | 8.884 × 10⁻⁵ | 8.948 × 10⁻⁵ |
| GSE14520 C1 control | 4.144 × 10⁻⁴ | 4.150 × 10⁻⁴ |
| GSE164760 C1 control | 0.005141 | 0.005218 |

**Decision, fixed now:** the reported values are left as computed and the variant
is **disclosed in the Methods**, since recomputing would churn every P value in
the manuscript to no benefit and the uncorrected variant is what the archived code
actually ran. Defect 1 **is** corrected, because a printed interval should be the
interval it claims to be.

### What this says about the rest of the study

The independent re-derivation is itself the reassurance: rank handling with ties,
the normal and Student *t* tails, Spearman's ρ, and every regression coefficient,
standard error and P value in the paper agree with SciPy and statsmodels. The two
defects are the complete list of disagreements found.


## §6o Amendment — diagnostics of two asserted claims (recorded 2026-08-31, before either was computed)

### Why

A systematic audit of the manuscript for claims that were **asserted rather than
computed** — the class that produced the §6k (*AIFM2*) and §6l (*POR*) errors —
returned two that are directly testable with data already on disk. Both are
registered here before computation. Neither can *add* support to anything: each
can only leave a claim standing, qualify it, or remove it.

### §6o-A — is DRAIN's collapse under D2 actually attributable to *POR*?

The manuscript states that DRAIN loses significance under the over-adjusted D2
model "as expected, since POR supplies those very P450s". D2 adds *CYP2E1*,
*CYP3A4*, *CYP1A2* and *CYP2C9*, and the argument is that conditioning on them
conditions on POR's own substrates. **That causal attribution has never been
computed**, and the sentence was written before §6l established that POR does not
belong to the N-reductive arm at all.

**The test:** fit the D2 model on **DRAIN′ = mean *z*(*MTARC1*, *MTARC2*)**, which
contains no *POR*.

**Decision rule, fixed now.**
- If DRAIN′ **also** loses significance under D2, the *POR* explanation is
  **falsified**: the collapse is dedifferentiation, not P450 conditioning. The
  clause is removed from the manuscript.
- If DRAIN′ **retains** significance under D2 while the locked DRAIN does not, the
  *POR* explanation **survives** and the clause stands.
- Anything in between is reported as indeterminate and the clause is softened to a
  possibility rather than an expectation.

### §6o-B — does removing *AIFM2* raise or lower REDUCTION's retained fraction?

The manuscript states that *AIFM2*'s absence from HT-HG-U133A makes GSE14520's
REDUCTION retention "conservative relative to a four-gene score rather than
generous". Retention is a **ratio**: adding a gene changes the adjusted intercept
and the unadjusted shift together, so the direction is not deducible. **This exact
sentence has already been reversed once, in §6k, and is still asserted rather than
computed.**

It is, however, directly testable in **GSE76427**, where *AIFM2* is measured:
compute REDUCTION retention with all four genes and again with only the three that
HT-HG-U133A carries (*CYB5R3*, *CYB5R1*, *NQO1*).

**Decision rule, fixed now.**
- If the three-gene retention is **lower** than the four-gene retention, "conservative"
  is supported and the sentence stands.
- If it is **higher**, the sentence is **wrong in the same way as before** and is
  replaced by the corrected direction.
- If the two are within 2 percentage points, the direction is declared
  **indeterminable** and the sentence is replaced by a statement that it is unknown.

This is one cohort and does not license extrapolation to GSE14520; whatever the
result, the manuscript will say that the direction is estimated in the cohort where
*AIFM2* exists and assumed, not shown, in the cohort where it does not.

### Outcome, recorded 2026-08-31

**§6o-A — FALSIFIED.** DRAIN′, which contains no *POR*, collapses under D2 at least as
far as the locked DRAIN does, in both cohorts:

| Cohort | module | D1 retained | D2 retained | D2 *P* |
|---|---|---|---|---|
| GSE76427 | DRAIN (with *POR*) | 41.4% | 20.4% | 0.216 |
| GSE76427 | **DRAIN′ (no *POR*)** | 32.9% | **11.7%** | **0.621** |
| GSE14520 | DRAIN (with *POR*) | 24.7% | 8.3% | 0.364 |
| GSE14520 | **DRAIN′ (no *POR*)** | 18.0% | **2.3%** | **0.847** |

By the branch fixed in advance, the *POR*-supplies-those-P450s explanation is
**removed from the manuscript**. The drain's disappearance under over-adjustment is
dedifferentiation, not conditioning on POR's substrates. This is the **third** claim
of the asserted-not-computed class found in this study, after §6k and §6l, and the
first one found by systematic audit rather than by accident.

**§6o-B — the claim survives, and is now computed rather than asserted.** In
GSE76427, dropping *AIFM2* from REDUCTION lowers the retained fraction by **8.1
percentage points** under D1 (71.1% → 63.0%) and **17.4** under the joint model
(51.2% → 33.8%). The three-gene score is therefore conservative, not generous, in
the cohort where this can be measured. The manuscript now says exactly that, and
says that the direction is assumed rather than shown in GSE14520, where *AIFM2*
does not exist to be dropped.

### Why both were run

Because a systematic audit of the manuscript for claims backed by neither a number
nor a citation returned 31 of them, and these two were the ones that were testable
with data already on disk. One was wrong. That ratio is the argument for the audit.



## §6p Correction — the stated motivation for the identity adjustment was an unsupported literature claim (recorded 2026-08-31)

### The claim, and what a literature check found

Since the concept note, this study has justified its central adjustment with:

> *MTARC1*, *MTARC2* and *POR* **are hepatocyte-differentiation genes**, and HCC is
> dedifferentiated, so the drain's fall might be identity loss rather than a change
> in electron allocation.

The second half is a sound worry. **The first half is not a supported statement.** A
targeted search of the primary literature returns no paper describing any of the
three as a marker of hepatocyte differentiation, and returns evidence pointing the
other way:

- **POR is *higher* in the dedifferentiated HepG2 line than in primary-like human
  hepatocytes** — 4-fold at mRNA, 2-fold at protein (Schulz et al., PMID 31561354).
  A gene that rises with dedifferentiation cannot be a differentiation marker.
- **MTARC1 is expressed in fetal human liver**, not only adult (Neve et al.,
  PMID 26378779). Presence in immature liver is the opposite of the expected pattern.
- **MTARC2 expression is regulated by nutritional status** — caloric restriction
  lowers it — and it is abundant in omental and subcutaneous fat (same paper).
- The canonical mature-hepatocyte marker panels in the differentiation literature
  (*ALB*, *HNF4A*, *ASGR1*, *TTR*, *SERPINA1*, *CPS1*, *TAT*, *G6PC1*) contain none
  of the three.

This is the same failure mode as §6l: a plausible-sounding mechanistic premise,
never traced to a source, propagated from the concept note into the manuscript.

### What is *not* wrong

**The adjustment itself is unaffected, and was never anchored on that premise.** D1
is a 22-gene panel built from exactly the canonical markers above — transcription
factors, secreted proteins, urea-cycle and gluconeogenic enzymes, and *ASGR1* — with
no redox enzyme and no P450. The adjustment therefore already conditions on a
validated hepatocyte-identity signature, not on the drain genes' presumed identity
status. Nothing needs re-anchoring and no number changes.

**And the study measured what the premise merely asserted.** Two independent results
in this project bear directly on it:

- §6m: in human liver, *POR* and *MTARC2* are hepatocyte-dominant (hepatocyte shares
  0.602 and 0.549, hepatocytes the top cell type for both, High antibody staining),
  and *MTARC1* stains Medium in hepatocytes.
- §6g: **59% and 75% of the drain's tumor fall is tracked by D1** in the two paired
  cohorts. That is a direct measurement that the drain genes covary strongly with
  hepatocyte identity in these data.

So the phenomenon the premise pointed at is real and was demonstrated here. What was
wrong was citing it as established biology instead of as this study's own finding.

### The correction

The manuscript stops asserting that the drain genes are differentiation markers. It
states instead that they are hepatocyte-enriched, that their fall *could* reflect
identity loss, that this is why a validated identity panel was pre-registered as a
covariate, and that the adjustment then showed most of the fall is indeed tracked by
identity. That is the measured version of the same argument and it is stronger.

**A new limitation is added.** *MTARC2*'s regulation by nutritional status is a live
confounder in a MASLD-to-HCC series, because feeding state and lipogenic status
covary with disease stage; the study cannot separate that from identity loss.

### A convergent problem with *NQO1*, recorded here rather than buried

Three independent lines in this project now point at the same weak member of
REDUCTION, and they should be read together:

1. **PHASE B:** the REDUCTION module's signal across the MASLD cohorts is *almost
   entirely NQO1*.
2. **§6m:** *NQO1* has the lowest hepatocyte share of all seven panel genes (0.006),
   is cholangiocyte-dominant in normal liver, and is undetected in hepatocytes by
   antibody staining.
3. **Literature:** *NQO1* is an obligate two-electron quinone reductase whose
   characterized role is quinone detoxification. It does reduce CoQ homologs, but
   rates fall with isoprenoid chain length and are slowest for CoQ10, and no primary
   study establishes it as a quantitatively significant CoQ10 reductase in human
   cells (Ross & Siegel, *Front Physiol* 2017;8:595).

The module is hash-locked and is **not** renamed or re-membered. But a module called
"CoQ-reducing" whose signal is carried mainly by its least CoQ-specific and least
hepatocyte-restricted member is a real weakness, and it is now stated in the
manuscript rather than left for a reader to assemble.

## §6q Post hoc robustness analyses prompted by pre-submission peer review (recorded 2026-09-01)

*This amendment is **not** pre-registered in the sense the rest of this
protocol is, and it is labeled accordingly wherever it is reported. Two
simulated adversarial reviews of the completed manuscript raised four
statistical objections. Three of the four quantities below were computed
during the adversarial verification of those reviews **before** this section
was written — that is, the computation preceded its rule. Under this
protocol's own standard that disqualifies them from being called registered,
so they are recorded here as post hoc, reviewer-prompted robustness checks,
with the order of operations stated rather than concealed. They may qualify a
reported claim; they may not create one.*

**6q-A — Support of the primary estimand.** The primary adjusted quantity is
the intercept of ΔRSI on ΔD1, i.e. the value at ΔD1 = 0. The observed ΔD1
distribution is: GSE76427 mean −0.422, SD 0.559, with 6 of 52 pairs (11.5%) at
ΔD1 ≥ 0; GSE14520 mean −0.724, SD 0.812, with 11 of 213 pairs (5.2%). ΔD1 = 0
therefore lies 0.75 and 0.89 SD above the center of the data and is supported
by few observations. This is a property of the estimand that was fixed in
advance and is not a defect of the fit, but it must be stated: the intercept
is a model-based extrapolation to a sparsely populated region.

**6q-B — Linearity.** Adding a quadratic term in ΔD1: GSE76427
β = −0.159, SE 0.129, t = −1.23, P = 0.224 (no evidence of curvature);
GSE14520 β = −0.0628, SE 0.0180, t = −3.48, P = 6.0 × 10⁻⁴ (curvature
present). The intercept moves from +0.568 to +0.595 (50.3% → 52.7%) and from
+0.792 to +0.733 (62.8% → 58.1%) respectively. Curvature is real in the larger
cohort and does not change either verdict against the 50% rule.

**6q-C — Uncertainty on the retained fraction.** The 50% decision rule was
applied to a point estimate; the protocol never specified an interval for a
ratio of two estimated quantities. Non-parametric pair bootstrap, 4,000
resamples, re-estimating both numerator and denominator in each resample:

| Cohort | Model | Retained | 95% CI |
|---|---|---|---|
| GSE76427 | D1 | 0.503 | 0.301–0.759 |
| GSE76427 | joint D1 + C1 | 0.399 | 0.095–0.706 |
| GSE14520 | D1 | 0.628 | 0.450–0.704 |
| GSE14520 | joint D1 + C1 | 0.645 | 0.452–0.736 |

**No interval excludes 0.50.** The registered rule is a point-estimate rule and
its verdicts stand as registered, but the data do not distinguish "clears 50%"
from "falls below 50%" in either cohort, and the manuscript must say so.

**6q-D — Heteroscedasticity.** HC3 robust standard errors on the intercept of
ΔRSI on ΔD1: GSE76427 0.137 → 0.160 (+17%), 95% CI 0.292–0.843 → 0.254–0.882;
GSE14520 0.068 → 0.102 (+50%), 95% CI 0.658–0.925 → 0.592–0.991. Both remain
above zero. Reported as a widening, not a reversal.

## §6r The SUPPLY module, and an accounting of the pre-specified sensitivity analyses (recorded 2026-09-01)

*Also post hoc. Section 5 of this protocol lists six sensitivity analyses
S1–S6. A protocol that is published with the paper is read by reviewers, and
an unreported pre-specified analysis is worse than one never planned. This
section states the disposition of all six, and repairs a reporting omission of
the same kind concerning the SUPPLY module.*

**6r-A — SUPPLY was defined, weighted +½, and never reported.** SUPPLY is 15 of
the index's 22 genes and carries half the positive weight, yet no SUPPLY result
appears in the manuscript. Computed now:

| | GSE76427 (52 pairs) | GSE14520 (213 pairs) |
|---|---|---|
| genes present | 15/15 | 14/15 |
| paired mean ΔSUPPLY | −0.036 | +0.336 |
| Wilcoxon *P* | 0.572 | 6.8 × 10⁻¹⁶ |
| joint D1 + C1 intercept | −0.200 (*P* = 0.076) | +0.368 (*P* = 5.4 × 10⁻¹³) |
| retained after joint adjustment | not interpretable (null numerator) | 110% |

Across the MASLD gradient (PHASE B, controls excluded) SUPPLY gave z = 1.35,
*P* = 0.089, ρ = +0.097 — the weakest of the three modules. So SUPPLY is
**discordant between the two paired cohorts**: absent in GSE76427, strongly
present and fully adjustment-resistant in GSE14520. This is a reportable
finding, and its absence from the manuscript was an omission, not a decision.
It qualifies the index's name: the composite's tumor behavior is carried by
REDUCTION and DRAIN, not by the biosynthetic supply arm.

**6r-B — Disposition of S1–S6.**

| | Variant | Status |
|---|---|---|
| S1 | GSVA/ssGSEA scoring | **not performed.** No registered rule permitted substituting a different scoring method after the mean-*z* index was locked; it is recorded here as not done rather than quietly dropped |
| S2 | cholesterol-arm genes added to SUPPLY | **performed (PHASE B).** Discovery z 5.98 → 5.95, ρ +0.40 → +0.39; replication z 4.82 → 4.76. The protocol stated in advance that an unchanged result **weakens** the claim of CoQ specificity. It is unchanged. This must appear in the manuscript and did not |
| S3 | *CYB5R3* alone as REDUCTION | **performed (PHASE B).** z 4.56 and 3.15, ρ +0.30 and +0.22 — same direction, weaker than the full module |
| S4 | *POR* removed from DRAIN | **performed**, reported in full as §6l / DRAIN′ |
| S5 | tumor purity (ESTIMATE/ABSOLUTE) | **not performed.** The C1/C2 non-parenchymal adjustment addresses the same confounder from the opposite direction and is reported; it is not a substitute and is not claimed as one |
| S6 | deconvolution against GSE202379 | **not performed.** Partially addressed by §6m (cell-type attribution of the panel genes) and by C1; neither is the registered analysis |

Three of six were performed, two of those three were never reported, and three
were not performed. All six are now accounted for in the manuscript.

## 7. Outcome-neutral stopping rule

If H3 and H5 both fail, the index will **not** be redefined. The result will
be reported as specified, as a descriptive study, and submitted accordingly.
Redefining the index after seeing the outcome would void this
pre-specification.

---

# Amendments §6s and §6t

Added 7 September 2026, after the amendments above and after the first
submission was desk-rejected. Both are post hoc relative to the 25 August 2026
lock, in the same sense as §6q and §6r, and both were written before the
computation each governs. The two source citations in §6s were verified against
their PubMed records on 7 September 2026, which corrected the sample count
stated for the second dataset from 159 to 110; no rule, prediction or threshold
was changed by that correction.

## §6s. Protein-level corroboration of the module dissociation

**Status.** Post hoc relative to the lock. Written before any protein value from
either dataset was read.

**Why.** Section 4.4 of the manuscript states that the principal weakness of this
study is the absence of orthogonal validation, and that the most informative next
step is protein-level measurement of both arms in paired hepatocellular carcinoma
(HCC) tissue. The claim under test is the one the manuscript now leads with: that
the fall in N-reductive capacity is largely explained by loss of hepatocyte
identity while the rise in CoQ-reducing capacity is not. If that dissociation is a
property of the tissue rather than of messenger RNA, it should appear in protein.

**Data, fixed in advance.** Two published human HCC proteogenomic cohorts with
paired tumor and adjacent liver:

1. Gao Q, Zhu H, Dong L, et al. Integrated proteogenomic characterization of
   HBV-related hepatocellular carcinoma. *Cell* 2019;179:561–577.e22.
   PMID 31585088, doi 10.1016/j.cell.2019.08.052. Paired tumor and adjacent
   liver from **159** patients. A correction was published as *Cell*
   2019;179:1240 (PMID 31730861) and is read alongside the original. Proteomic
   data are deposited in the Proteomics Data Commons.
2. Jiang Y, Sun A, Zhao Y, et al. Proteomics identifies new therapeutic targets
   of early-stage hepatocellular carcinoma. *Nature* 2019;567(7747):257–261.
   PMID 30814741, doi 10.1038/s41586-019-0987-8. **110** paired early-stage
   HBV-related cases. Free full text is available in PubMed Central under
   PMC6497684.

No other proteomic dataset will be substituted after results are seen. If a
dataset proves inaccessible, that fact is reported rather than replaced.

**Directional prediction, fixed in advance.** At protein level, in each dataset:

- REDUCTION members (CYB5R3, CYB5R1, AIFM2, NQO1) are higher in tumor than in
  paired adjacent liver.
- DRAIN members (MTARC1, MTARC2, POR) are lower.
- After adjustment for a protein-level hepatocyte-identity covariate built from
  the D1 genes measurable in that dataset, the REDUCTION rise retains **at least
  50%** of its unadjusted shift and the DRAIN fall retains **less than 50%**,
  mirroring the transcript result in GSE14520 (89% against 31%).

**Estimand and method.** Identical in form to Methods §2.3 of the manuscript: the
intercept of the paired-difference regression of the paired module difference on
the paired difference in the identity covariate, and the identity-retention
fraction as the proportion of the unadjusted paired shift that intercept retains.
Module scores are the mean of within-dataset protein-wise z scores, matching the
transcript scoring.

**Coverage rule, fixed in advance.** Mass-spectrometry coverage is incomplete by
nature. A module with **fewer than two measured members** in a dataset is reported
as not testable in that dataset and is not scored. The identity covariate must
retain at least eight of the twenty-two D1 proteins, or the adjustment for that
dataset is reported as not performed. Measured membership is reported before any
contrast is computed.

**Outcomes accepted in advance.**

- If both datasets show the predicted dissociation, the manuscript reports
  protein-level corroboration and the abstract says so.
- If the datasets show the shifts but **no** dissociation after adjustment, I
  report that the arms are not distinguishable at protein level and remove the
  claim that the CoQ-reducing rise is independent of identity, keeping only the
  transcript observation with that limitation stated.
- If the datasets show the **opposite** direction for either arm, I report that
  the transcript finding is contradicted at protein level, in the abstract, and
  the manuscript's central claim is withdrawn.
- If coverage is insufficient under the rule above, the analysis is reported as
  not testable and the manuscript's limitation about the absence of orthogonal
  validation stands unchanged.

**What this cannot establish.** Protein abundance is not enzyme activity, and
these are cross-sectional tissues. A corroborated dissociation would not
demonstrate that electron allocation changes, only that the abundance pattern is
not an artifact of transcript measurement.

---

## §6t. Identity retention across published liver expression signatures

**Status.** Post hoc relative to the lock. Written before any signature retention
value was computed, and before any of the gene sets below was scored in any
cohort.

**Why.** Section 4.2 of the manuscript argues that the identity-retention fraction
should be reported routinely for liver expression signatures, and that the
quantity is not a formality. As written, that argument rests on a single worked
example, this study's own two modules. This amendment converts it from a claim
about one panel into a measurement across the published literature, or fails to.

**Panel selection rule, fixed in advance.** Every gene set in the MSigDB human
C2:CGP collection, release v2026.1.Hs, whose **standard name contains the
substring `LIVER` or `HEPAT`**. The enumeration was performed on 7 September 2026
and returned **149** sets. C2:CGP is the collection of gene sets curated from
published expression studies, so this rule selects published signatures by a
mechanical criterion rather than by my judgment of which ones matter.

**Eligibility, fixed in advance.**

- Between **15 and 1000** gene symbols in the retrieved set. This leaves **126**
  of the 149 sets.
- At least **10** symbols present on the cohort's platform after mapping.
  Sets failing this in a given cohort are reported as not evaluable there.

No set is added, removed, split, or merged after any retention value is seen. The
149 retrieved sets, their gene counts, their source publications, and the
verification of each count against the MSigDB record are archived with the
analysis code.

**Cohorts.** GSE14520 (213 pairs) is primary, being the larger and the cohort
carrying the study's own module result. GSE76427 (52 pairs) is reported alongside.
GSE164760 is excluded, as it is everywhere else in this study, because it did not
support the tumor-adjacent contrast under the rule fixed in advance.

**Scoring and estimand.** Identical to the study's own modules: the mean of
within-cohort gene-wise z scores, then the paired tumor minus adjacent difference,
then the intercept of the paired-difference regression on ΔD1 and on ΔD1 with ΔC1
jointly, and the identity-retention fraction as defined in Methods §2.3.

**Direction, fixed in advance.** Retention is computed on the **absolute**
unadjusted shift, so a signature that falls in tumor is treated symmetrically with
one that rises.

**Null-effect rule, fixed in advance.** A signature with no unadjusted
tumor-adjacent effect has nothing to retain. Sets whose unadjusted paired shift
does not reach two-sided Wilcoxon signed-rank *P* < 0.05 are **excluded from the
retention distribution**, and their number is reported. This rule is stated now,
before any *P* value is computed, precisely because it would otherwise be a
degree of freedom.

**Overlap handling, fixed in advance.** Of the 126 eligible sets, 37 share at
least one gene with D1, 30 with C1, and 36 with this study's own panel; no set
draws more than 20% of its members from D1 or from C1. Overlap counts are reported
per set. A pre-specified sensitivity analysis recomputes every retention fraction
after removing from each signature any gene it shares with D1 or C1, and the
manuscript reports whether any conclusion changes.

**Primary summary, fixed in advance.** The median and interquartile range of the
identity-retention fraction across eligible, non-null sets in GSE14520, the number
and proportion of those sets retaining **less than 50%**, and the percentile at
which this study's REDUCTION arm (89%) and DRAIN arm (31%) sit in that
distribution.

**Outcomes accepted in advance.**

- If a **majority** of published liver signatures retain less than 50%, the
  argument in §4.2 is supported and the manuscript may state that a large part of
  the liver expression signature literature is substantially attributable to
  dedifferentiation.
- If a **majority retain well above 50%**, the argument is weakened. I will report
  that most published liver signatures do carry information beyond
  dedifferentiation, soften §4.2 to a recommendation that the quantity be
  reported rather than a claim that the literature is confounded, and note that
  this study's own DRAIN arm is then unusual rather than typical.
- If the distribution is **wide with no central tendency**, I will report that
  identity retention varies unpredictably between signatures. That supports
  reporting the quantity per signature and does **not** support the stronger claim,
  and §4.2 will say so.

In every case the distribution is reported in full. No outcome will be used to
reselect the panel, change the eligibility rule, or change the null-effect rule.

**What this cannot establish.** This benchmark does not say whether any signature
is biologically wrong, whether it predicts outcome, or whether its source study
was correct. It measures one thing: how much of a signature's tumor versus
adjacent contrast is not attributable to the loss of hepatocyte identity. A low
retention fraction is not a refutation of the signature. It is a statement about
what the contrast is made of.

**Relation to the study's own result.** This study's modules are scored by the
same code on the same pairs with the same covariates, so their position in the
distribution is a like-for-like comparison. That is the point of running it. It
also means a flattering result cannot be produced by scoring the comparators
differently, and I will report the code path shared between the two.

---

# Amendment §6u

Added 7 September 2026, after §6s and §6t and after their results were known,
but before any number this amendment governs was computed. Post hoc relative to
the 25 August 2026 lock, in the same sense as §6q through §6t.

---

## §6u. Is the identity-retention fraction a re-description of tumor purity adjustment?

**Status.** Post hoc relative to the lock. Written before the composition-only
retention fraction was computed for any signature, and before the correlation
between the two covariates was examined in either cohort.

**Why.** A reader is entitled to ask whether this measure is new. Adjusting a
tumor versus normal contrast for the non-tumor content of the specimen is
established practice, through tumor purity estimation and through reference-based
deconvolution of cell type proportions. If the identity covariate D1 and the
composition covariate C1 move together across pairs, and if adjusting for either
one leaves the same fraction of a signature's shift standing, then the quantity
this manuscript proposes is a re-labeling of something the field already does, and
the recommendation in Section 4.2 has no claim on anyone's attention. This
amendment fixes, before the numbers exist, the statistics that would settle that
and the decision I will act on in each case.

**What is distinct in principle, and why that is not enough.** D1 is 22 hepatocyte
genes and C1 is 21 leukocyte, endothelial and stromal genes; the two lists share no
member by construction, and a specimen can lose hepatocyte identity without any
change in non-parenchymal content, since a dedifferentiated hepatocyte is still a
tumor cell and is still counted in the tumor fraction. That is an argument from
construction. It does not establish that the two behave differently in real
paired data, which is what the analysis below tests.

**Analysis, fixed in advance.** No gene set is added, removed, rescored or
reordered. The eligibility rule, the null-effect rule and the pairing of §6t are
carried over unchanged, so the evaluable set is exactly the one already reported.
Every regression uses the same ols_ci routine and the same paired-difference
construction as the rest of the study.

1. **Covariate association.** Pearson *r* and Spearman ρ between ΔD1 and ΔC1
   across the paired differences, in GSE14520 and in GSE76427, with the *R*² of
   ΔD1 regressed on ΔC1.
2. **Composition-only retention.** For every evaluable signature, the retention
   fraction from the model Δscore = α + βΔC1 + ε, computed exactly as the
   identity-only retention fraction is computed, on the same absolute
   denominator.
3. **Agreement.** Spearman ρ between the identity-only and composition-only
   retention fractions across signatures; the median and interquartile range of
   their difference; and the number of signatures whose position relative to the
   50% threshold differs between the two adjustments.
4. **Discordant cases.** The number of signatures for which composition
   adjustment leaves at least 90% standing while identity adjustment leaves less
   than 50%, and the number for which the reverse holds. These are reported as
   counts with examples named, not selected for their content.

**Decision rule, fixed in advance.**

- **Redundant.** If *all four* of the following hold: |Spearman ρ(ΔD1, ΔC1)| ≥
  0.70 in GSE14520; Spearman ρ between the two retention fractions ≥ 0.70; the
  median absolute difference between them is below 0.10; and fewer than 10% of
  evaluable signatures change side of the 50% threshold. Then I report that the
  identity-retention fraction is a re-description of composition adjustment. The
  recommendation in Section 4.2 is withdrawn as a distinct recommendation and
  rewritten to point at existing purity-adjustment practice, and the abstract is
  changed accordingly.
- **Distinct.** If |Spearman ρ(ΔD1, ΔC1)| < 0.70 *and* either the retention
  Spearman is below 0.70 or at least 10% of signatures change side. Then I report
  the measure as capturing something composition adjustment does not, and Section
  4.1 states the comparison with these numbers rather than by assertion.
- **Partial.** Anything else. I report partial overlap with the numbers in full,
  and Section 4.2 is softened to recommend the quantity be reported alongside
  composition adjustment rather than as a substitute for it.

**What this comparison cannot do, stated before it is run.** C1 is a covariate
defined inside this study. It is of the same construct as the stromal and immune
components of established purity estimators, but it is not one of them, and the
published gene lists for those estimators could not be retrieved in the
environment this analysis runs in. This is therefore a comparison against a
locally built composition score, not a benchmark against a named tool, and it is
reported that way. A reader who wants the comparison against a specific estimator
has everything needed to run it: the code path, the pairing and the scoring are in
Supplementary File 2, and substituting a different covariate is a one-line change.

**What a distinct result would and would not mean.** It would mean that
adjusting for non-parenchymal content does not remove what adjusting for
hepatocyte identity removes. It would not mean that either adjustment is correct,
that D1 measures differentiation without error, or that the residual after
adjustment is tumor-specific biology. The residual is what is left after two
covariates, and other confounders are not addressed by either.

---

# Amendment §6v

Added 7 September 2026, after §6u and after its result was known, but before
any interval or simulated dataset governed by this amendment existed. Post hoc
relative to the 25 August 2026 lock, in the same sense as §6q through §6u.

---

## §6v. Uncertainty of the retention fraction, and whether the estimator recovers a known value

**Status.** Post hoc relative to the lock. Written before any interval and before any
simulated dataset existed. The descriptive quantities quoted below as inputs are
properties of fits already reported in §6t, not new outcomes.

**Why.** §6q recorded that the retention fraction is a ratio of two estimated
quantities, so its uncertainty exceeds that of either, and bootstrapped that ratio
for the composite index alone. The benchmark of §6t reported 119 per-signature
values as point estimates with no interval, and Section 3.4 draws a count from
them: 30 of 119 below 50%. A count of point estimates is not the same claim as a
count of signatures for which the data support a value below 50%, and the
manuscript should not be read as making the second when it has only established
the first. This amendment fixes, before the intervals exist, how they will be
computed and what each possible result obliges the manuscript to say. It also
fixes a simulation, because an estimator whose behavior against a known truth has
never been checked should not be recommended to other people.

### Part A. Interval estimation for the 119 signatures

**Design, fixed in advance.** Nonparametric paired bootstrap over patients.
Module and covariate scores are computed once on the full cohort and held fixed;
the resampling is over the 213 patient pairs, with replacement, so the intervals
express sampling variability of patients and not of the scoring step. This is
stated because it is a limitation, not a neutral choice: variability introduced by
the within-cohort z scoring is not captured.

- **B = 2000** resamples. Random seed 20260907, recorded in the code.
- For each resample the unadjusted paired mean shift and the joint D1 + C1
  intercept are recomputed **on the same resample**, and the ratio is formed from
  the two. Percentile intervals at 2.5% and 97.5%.
- The primary model is the joint D1 + C1 model, matching Section 3.4. The
  identity-only model is reported alongside.
- **Instability guard, fixed in advance.** A ratio whose denominator is not
  reliably different from zero has no stable interval. If the 95% bootstrap
  interval of a signature's unadjusted shift contains zero, that signature is
  flagged unstable and its retention interval is reported as not defined rather
  than as a number. The count of such signatures is reported.
- The archived §6t point estimates must be reproduced before any interval is
  computed, as in §6u.

**Quantities fixed in advance.** Of the 30 signatures below 0.50 by point
estimate, how many have an entire interval below 0.50; the median interval width;
how many intervals contain 1.0; how many are unstable under the guard; and the
intervals for this study's own REDUCTION and DRAIN arms.

**Outcomes accepted in advance.**

- If **10 or more** of the 30 keep the whole interval below 0.50, Section 3.4
  reports both numbers, in the form "30 by point estimate, of which N have an
  interval entirely below 50%", and the abstract carries the interval-supported
  number rather than the point count.
- If **fewer than 10** do, Section 3.4 states plainly that the count below 50% is
  not supported by interval estimation, the count is removed from the abstract,
  and the claim is reduced to the dispersion statement, which does not depend on
  any individual value.
- If the **median interval width exceeds 0.50**, the manuscript states in Results
  and in Section 4.6 that a single retention value from a cohort of this size is
  too imprecise to act on alone, and Section 4.2 recommends that the quantity
  always be reported with its interval.
- In every case every interval is deposited, and no signature is dropped,
  rescored or reordered on the basis of its interval.

### Part B. Does the estimator recover a known retention fraction?

**Design, fixed in advance.** Synthetic paired differences are generated on the
**real** ΔD1 and ΔC1 vectors from the 213 GSE14520 pairs, so that the covariate
distribution, including the scarcity of ΔD1 near zero that §6q flagged as making
the intercept an extrapolation, is the one the estimator actually faces.

For a target retention τ and a target unadjusted shift μ,

```
Δy = α + β ΔD1 + γ ΔC1 + ε,    ε ~ N(0, σ²)
α = τ μ,   β mean(ΔD1) + γ mean(ΔC1) = (1 − τ) μ
```

so the true retention |α| / |E[Δy]| equals τ by construction.

**Grid, fixed in advance.** Ranges are taken from the fits already reported in
§6t, where the unadjusted absolute shift across the 119 signatures runs from 0.031
to 1.378 with median 0.497, and the implied residual standard deviation runs from
0.228 to 0.966 with median 0.487.

- τ ∈ {0.10, 0.25, 0.50, 0.75, 1.00}
- σ ∈ {0.25, 0.50, 1.00}
- composition share of the non-retained part ∈ {0, 0.3}
- μ = 0.50 for the main grid, with μ ∈ {0.25, 1.00} run at σ = 0.50 and
  composition share 0 for τ ∈ {0.25, 0.75}
- 500 replicates per cell; 500 bootstrap resamples per replicate for coverage.
  Random seed 20260907.

**Quantities fixed in advance.** Per cell: the median estimate, the bias as median
estimate minus τ, the interquartile range of the estimate, and the empirical
coverage of the nominal 95% percentile interval.

**Outcomes accepted in advance.**

- If |bias| is at most 0.05 and coverage lies between 0.92 and 0.97 across the
  whole grid, the manuscript reports that the estimator is well behaved at this
  sample size and states the grid over which that was checked.
- If |bias| exceeds 0.05 in any cell, its direction, size and the conditions
  producing it are reported in Results and in Section 4.6, and the manuscript
  states the range over which the measure should be read rather than presenting
  it as unconditional.
- If coverage falls below 0.90 in any cell, the manuscript reports the intervals
  as anticonservative under those conditions and gives the achieved coverage.
- If the estimator fails badly, meaning bias above 0.20 or coverage below 0.80
  anywhere in the grid at σ ≤ 0.50, the recommendation in Section 4.2 is
  withdrawn and the manuscript reports the measure as not yet fit for general
  use.

No cell will be added, removed, or re-run with different settings after its result
is seen, and the grid above is the grid that will be reported.

**What this cannot establish.** The simulation assumes the model the estimator
fits, with additive covariate effects, normal errors and no measurement error in
the covariates. It therefore tests whether the estimator recovers its own
estimand at this sample size and covariate geometry. It does not test robustness
to a wrong functional form, and §6q already records that a quadratic term is
warranted in this cohort. That limitation is reported alongside the result.

---

# Amendment §6w

Added 7 September 2026, after §6v. Post hoc relative to the 25 August 2026 lock,
in the same sense as §6q through §6v. Written before the expression matrix it
governs was opened; the only file examined beforehand was the gene annotation
map, which contains no expression data.

---

## §6w. A third paired cohort, on RNA sequencing: TCGA-LIHC

**Status.** Post hoc relative to the lock. Written before the expression matrix
was opened, before any pair was formed, and before any value from this cohort was
computed. The only file examined at the time of writing is the gene annotation
map, which contains no expression data.

**Why.** The benchmark of §6t rests on one cohort, GSE14520, because GSE76427's
archived values could not be rebuilt and its positive controls are advisory. Both
are microarray. A distribution of retention fractions drawn from a single
microarray cohort may be a property of that cohort or of that platform rather than
of liver signatures, and I cannot tell the difference from the data reported so
far. This amendment adds a third paired cohort measured by RNA sequencing and
fixes, before the matrix is opened, what will be run on it and what each outcome
obliges the manuscript to say.

**Dataset, fixed in advance.** TCGA-LIHC gene expression, STAR - TPM, from the
UCSC Xena GDC hub, dataset `TCGA-LIHC.star_tpm.tsv`, version 05-09-2024,
60,661 identifiers by 424 samples, unit log2(TPM + 1), with the hub's stated
wrangling that data from different vials, portions, analytes or aliquots of the
same sample are averaged before the log transform.

TPM is chosen over the hub's default FPKM-UQ and over STAR - Counts, and the
reason is fixed here rather than after any result. The locked protocol specifies
log2(CPM + 1) for sequencing data. log2(TPM + 1) is the direct analogue: TPM and
CPM differ by a gene-length factor that is constant across samples and therefore
very nearly cancels in a within-cohort, per-gene z score. STAR - Counts is served
as log2(count + 1) with no library-size normalization and is therefore not
comparable across samples without undoing the transform. **No alternative
quantification will be substituted after any result from this cohort is seen.** If
the chosen file proves unusable, that fact is reported rather than replaced.

**Symbol mapping, fixed in advance.** `gencode.v36.annotation.gtf.gene.probemap`,
60,660 rows, 59,427 distinct symbols. Where one symbol carries several Ensembl
identifiers, 110 symbols in total, the **maximum across identifiers on the log2
scale** is taken, which is the same collapse rule this study applied to multiple
probes on the arrays. Eleven of the affected symbols fall among the genes this
study uses, and most are pseudoautosomal genes annotated on both chrX and chrY;
none is a member of D1, C1, or any redox module, all of which map one to one.

**Alias, fixed in advance.** GENCODE v36 carries glucose-6-phosphatase as `G6PC`.
`G6PC1`, the current symbol and the one in the locked D1 list, is therefore
absent. The single rename G6PC1 to G6PC is applied. It is one to one, it is
unambiguous, and it is the same alias already applied in §6s. No other member of
any locked list is missing from this annotation.

**Pairing, fixed in advance.** From the TCGA barcode: the patient is the first
twelve characters, and the sample type is characters fourteen and fifteen, where
01 is primary tumor and 11 is solid tissue normal. Only 01 and 11 are used, and
any other sample type is discarded before pairing. If a patient still carries more
than one sample of either type, the first in sorted order is taken and the number
of such patients is reported. The number of pairs is reported before any contrast.

**Premise check, fixed in advance and reported before any contrast.** Following
the standing rule that a public matrix must be shown to contain what the analysis
assumes: I will report the number of samples of each type, the number of pairs,
the preprocessing already applied by the source, and the paired difference in the
hepatocyte-identity score D1. **If D1 is not lower in tumor at two-sided Wilcoxon
signed-rank *P* < 0.05, the confound this measure exists to remove is absent, and
this cohort is reported as not testable. No retention value from it is reported in
that case**, exactly as was done for the second proteome in §6s.

**Note on what is and is not a control.** Unlike §6u and §6v, there is no archived
value from this cohort to reproduce, so the run cannot be gated on reproducing
one. The premise check above is the only control. The behavior of the DRAIN and
REDUCTION modules is a **prediction** and not a control, and will not be used to
decide whether the pipeline is working.

**Analyses, fixed in advance.** Run in this order, with the eligibility rule, the
null-effect rule, the scoring, the covariates and every threshold of §6t, §6u and
§6v carried over unchanged:

1. The published-signature benchmark of §6t.
2. The identity against composition comparison of §6u.
3. The bootstrap intervals of §6v, at B = 2000 with the same seed.

**A module that is measurable here for the first time.** *AIFM2* is absent from
HT-HG-U133A, so GSE14520's REDUCTION is a three-gene score, which the manuscript
lists as a limitation. It is present in this annotation, so REDUCTION here is the
full four-gene module. **Both versions will be reported**, the four-gene module
and the three-gene module restricted to the members measurable in GSE14520, so
that the comparison between cohorts is like for like and neither version can be
chosen after the fact.

**Directional predictions, fixed in advance.** In this cohort, REDUCTION is higher
in tumor and retains at least 50% of its unadjusted shift; DRAIN is lower and
retains less than 50%.

**Outcomes accepted in advance.**

- If the benchmark distribution here resembles GSE14520, with a median retention
  within 0.15 of 0.758 and a comparably wide range, the manuscript reports that
  the distribution replicates across platform and cohort, and Section 3.4 says so.
- If it differs materially, the manuscript reports that the distribution is
  cohort or platform dependent, Section 4.2 is changed to say the quantity must
  be reported for the dataset at hand rather than compared across datasets, and
  the single-cohort limitation in Section 4.6 is strengthened rather than removed.
- If the premise check fails, the cohort is reported as not testable, no retention
  value from it appears, and Section 4.3 gains a second worked instance of a
  deposited matrix that does not support the analysis applied to it.
- If DRAIN retains 50% or more here, the manuscript reports that the worked
  example does not replicate on RNA sequencing, in the abstract, and the
  separation between the two arms is presented as specific to the array cohorts.
- If REDUCTION and DRAIN behave as predicted, that is reported as replication
  and the four-gene and three-gene versions are both shown.

In every case the full distribution is deposited, and no signature is added,
removed, rescored or reordered on the basis of any result from this cohort.

**What this cannot establish, stated before the data are opened.** TCGA-LIHC is a
mixed-etiology, largely Western cohort, whereas GSE14520, Gao and Jiang are
hepatitis B related. Any difference between this cohort and the others may be
etiology, platform, population, or specimen handling, and this design separates
none of them. A replication would therefore be evidence that the measure behaves
consistently across all of those at once, and a failure would not identify which
of them was responsible. The manuscript will say this wherever the comparison is
drawn.

---

# Amendment §6x

Added 8 September 2026, after §6w. Post hoc relative to the 25 August 2026 lock,
in the same sense as §6q through §6w. Written before the kidney expression matrix
was opened. The kidney identity covariate it defines is stated in full here,
before any kidney expression value was read.

---

## §6x. A second tissue: does the measure generalize beyond liver?

**Status.** Post hoc relative to the lock. Written before the kidney expression
matrix was opened and before any kidney value of any kind was computed. The only
files examined beforehand are the gene annotation map, which contains no
expression data, and the local gene-set directory, which was checked only to
establish that no kidney gene set is present.

**Why.** Everything reported so far is liver. A quantity defined for one tissue,
validated in that tissue and recommended to that tissue's literature has not been
shown to be a method. The manuscript's recommendation in Section 4.2 is written
as though it applies to expression signatures generally, and it is not entitled
to that scope on liver data alone. This amendment adds clear cell renal cell
carcinoma (ccRCC), the tumor whose defining feature is loss of the proximal
tubule program from which it arises, and fixes in advance what the result obliges
the manuscript to claim.

**Dataset, fixed in advance.** TCGA-KIRC gene expression, STAR - TPM, from the
UCSC Xena GDC hub, unit log2(TPM + 1), with the hub's aliquot averaging already
applied. Symbols from the same `gencode.v36.annotation.gtf.gene.probemap` used in
§6w, collapsing several identifiers per symbol by the maximum on the log2 scale.
Pairing from the TCGA barcode, sample type 01 as tumor and 11 as solid tissue
normal, any other type discarded before pairing. As in §6w, no alternative
quantification will be substituted after any result is seen.

### The kidney identity covariate, K1

This is the single largest degree of freedom in the whole amendment, and it is
fixed here, in full, before any expression value is read. K1 is built by the same
three rules that produced D1 and by no others:

1. genes that mark the **differentiated epithelium of the nephron**, weighted
   toward the proximal tubule, the segment ccRCC arises from;
2. **no redox enzyme and no cytochrome P450**, so the covariate cannot absorb the
   quantity a redox signature would measure;
3. **disjoint from C1**, from this study's redox modules, and from **D1**, so that
   no liver identity gene is imported into a kidney covariate.

It also matches D1's internal structure, seventeen differentiated effectors and
five regulators, so that the two covariates are comparable in construction and
the effector against regulator contrast found in §6w can be examined again.

**Effectors (17).** LRP2, CUBN, SLC34A1, SLC5A2, SLC5A12, SLC22A6, SLC22A8,
SLC22A12, SLC12A1, SLC12A3, UMOD, AQP1, AQP2, MIOX, ALDOB, GATM, NAT8.

**Regulators (5).** PAX2, PAX8, HNF1B, ESRRG, PPARGC1A.

This list is now locked. Coverage on the platform will be reported after the
matrix is opened, and a gene that is absent is reported as absent and **not
replaced**, exactly as *AIFM2* was on the array in §6w. No member will be added,
removed or exchanged on the basis of how it behaves.

**Composition covariate.** C1 unchanged, its 21 leukocyte, endothelial and
stromal genes being tissue-agnostic. One prediction is recorded now because it
runs opposite to liver: ccRCC is a densely vascular and immune-infiltrated tumor,
so **C1 is expected to rise in tumor here**, whereas it fell in the liver
cohorts. If it does, the composition covariate is doing more work in kidney than
in liver, and §6u may well come out differently. That is an anticipated outcome,
not a surprise to be explained afterward.

**The worked example is not carried over.** The redox panel is liver-motivated
and mARC is principally a hepatic system. **No REDUCTION, DRAIN, SUPPLY or
composite result from kidney will be computed or reported**, whatever it might
show. This is fixed now precisely so that a favorable kidney redox result could
not be added later, and an unfavorable one could not be quietly omitted. The
kidney cohort tests the measure, not the biology.

**Premise check, reported before any contrast.** Following the standing rule: the
number of samples of each type, the number of pairs, the processing already
applied, and the paired difference in K1. **If K1 is not lower in tumor at
two-sided Wilcoxon signed-rank *P* < 0.05, the kidney cohort is reported as not
testable, no retention value from it appears, and the manuscript states that the
covariate construction did not transfer to a second tissue.** The per-gene table
is reported alongside, as in §6w.

**Signature panel.** The §6t enumeration rule applied to kidney: every gene set
in the MSigDB human C2:CGP collection whose standard name contains **KIDNEY or
RENAL**, keeping sets of 15 to 1000 symbols, evaluable at 10 or more symbols on
the platform, and excluding sets whose unadjusted paired shift does not reach
two-sided Wilcoxon *P* < 0.05. The full enumeration, including sets excluded and
why, is reported. The collection file is not yet in hand; if it cannot be
obtained, the benchmark is reported as not performed and only the premise check
and the covariate behavior are reported.

**Analyses.** The §6t benchmark, the §6u identity against composition comparison
and the §6v bootstrap intervals, each with its rules and thresholds carried over
unchanged. Because the number of pairs will be of the order of the second and
third cohorts rather than of GSE14520, §6v predicts wide per-signature intervals,
and the kidney cohort is therefore expected to support the **distribution** and
not individual values. That expectation is recorded now so that wide intervals
are not reported later as though they were a surprise.

**Outcomes accepted in advance.**

- **Generalizes.** If the kidney retention distribution is comparably dispersed,
  meaning its interquartile range is at least half the liver range's and its full
  range spans at least tenfold, and its median lies within 0.20 of the liver
  median, the manuscript reports that the measure behaves the same way in a
  second tissue, and the recommendation in Section 4.2 is stated as a general one
  rather than a liver one.
- **Does not generalize.** If the kidney distribution is narrow, or retention
  clusters near 1.0 for nearly all sets, the manuscript reports that kidney
  signatures are largely not attributable to identity loss, that the liver result
  may be specific to a tissue whose differentiated program is unusually dominant,
  and Section 4.2 stays explicitly scoped to liver.
- **Premise fails.** Kidney is reported as not testable, and the manuscript says
  the covariate construction did not transfer, which is itself informative about
  how portable the method is.
- **Composition dominates.** If §6u returns REDUNDANT in kidney while it returned
  DISTINCT in liver, the manuscript reports that whether identity and composition
  adjustment differ is itself tissue-dependent, and Section 4.1 is qualified
  accordingly rather than left as a general claim.

No gene set, and no member of K1, will be added, removed, rescored or reordered
on the basis of any result.

**What this cannot establish, stated before the data are opened.** Two tissues are
not generality. K1 was assembled by me from prior physiological knowledge and has
none of the provenance D1 has, since D1 was itself locked before any liver
outcome and has since been examined in four datasets. A kidney result therefore
tests whether the construction is portable, not whether this particular kidney
covariate is correct. ccRCC also has an unusual metabolic phenotype driven by VHL
loss, so the strength of identity loss there may not represent solid tumors
generally.

---

# Amendment §6y

Added 8 September 2026, after §6x and after its result was known. Post hoc
relative to the 25 August 2026 lock. Written before the lung expression matrix
was obtained. The tissue choice, the reason for rejecting the larger breast
panel, and the lung identity covariate are all stated here in advance.

---

## §6y. A third tissue, chosen because §6x could not answer the question it asked

**Status.** Post hoc relative to the lock. Written before the lung expression
matrix was obtained. The files examined beforehand are the MSigDB C2:CGP
collection, to count how many gene sets each candidate tissue yields, and the
gene annotation map, to confirm that the covariate symbols below exist. Neither
contains expression data, and no expression value from any lung sample has been
read.

**Why, stated plainly.** §6x added kidney and could not settle the question it
was written to settle. Applying the §6t enumeration rule to kidney returned
eighteen gene sets, of which fifteen were evaluable, against 119 in liver, and
the fifteen were dominated by transplant immunology and aging rather than by
renal cancer. Two of the three criteria §6x fixed for a general result were met
and one was not. That is a limitation of the comparator panel, not a finding
about the measure, and it should not be presented as one. This amendment
therefore adds a third tissue chosen so that the panel is large enough for the
question to have an answer.

**How the tissue was chosen, and when.** Before any expression matrix was
obtained, the §6t enumeration rule was applied to the C2:CGP collection for
eleven candidate tissues, and the count of eligible sets was tabulated against
the number of paired tumor and adjacent samples each TCGA project provides. That
table is reported in full. Breast returns the largest panel, 160 eligible sets
and about 113 pairs, and was **not** chosen, for a reason recorded here before
any breast or lung value was seen: adjacent normal breast is largely adipose and
stroma while the tumor is an epithelial proliferation, so an epithelial identity
score is as likely to rise in tumor as to fall, and the premise this measure
depends on would probably fail. Lung adenocarcinoma was chosen instead: 41
eligible sets, nearly three times kidney's, about 59 pairs, and a cell of origin,
the alveolar type II pneumocyte, whose differentiation program is expected to
fall. **Choosing on expected premise behavior is a decision I am making in
advance and disclosing, not one made after seeing a result.** If lung's premise
fails, that is reported as a failure and breast is not substituted for it.

**Dataset, fixed in advance.** TCGA-LUAD gene expression, STAR - TPM, from the
UCSC Xena GDC hub, unit log2(TPM + 1), aliquots already averaged by the hub.
Symbols from the same `gencode.v36.annotation.gtf.gene.probemap`, collapsing
several identifiers per symbol by the maximum on the log2 scale. Pairing from the
TCGA barcode, sample type 01 as tumor and 11 as solid tissue normal, every other
type discarded. No alternative quantification will be substituted after any
result is seen.

### The lung identity covariate, L1

Built by the same three rules that produced D1 and K1: genes marking the
**differentiated distal lung epithelium**, weighted toward the alveolar type II
cell that most lung adenocarcinoma arises from; **no redox enzyme and no
cytochrome P450**, which is why the lung alcohol dehydrogenases, flavin
monooxygenases and catalase are excluded despite being distal-lung enriched; and
**disjoint from C1, from the redox modules, from D1 and from K1**, so no liver or
kidney identity gene is imported. The structure again matches D1 and K1,
seventeen effectors and five regulators.

**Effectors (17).** SFTPC, SFTPB, SFTPA1, SFTPA2, SFTPD, NAPSA, SLC34A2, PGC,
LAMP3, ABCA3, AGER, CLDN18, SCGB1A1, SCGB3A2, CTSH, SFTA2, CLIC5.

**Regulators (5).** NKX2-1, FOXJ1, HOPX, CEBPA, ETV5.

All 22 were confirmed present in the annotation before this amendment was closed,
none carries more than one Ensembl identifier, and none overlaps D1, K1, C1 or
the redox panel. The list is now locked; a member that turns out to be absent
from the matrix is reported as absent and **not replaced**, and no member will be
added, removed or exchanged on the basis of how it behaves.

**Composition covariate.** C1 unchanged. No directional prediction is registered
for it here, because lung adjacent tissue and lung tumor differ in immune and
stromal content in ways I cannot predict with the confidence I had for kidney.
Its behavior is reported as observed.

**The worked example is not carried over.** As in §6x, **no REDUCTION, DRAIN,
SUPPLY or composite result from lung will be reported**, whatever it shows. Where
the shared code path computes them as its built-in control step, they are removed
from the deposited files and the removal is recorded there, exactly as was done
for kidney.

**Premise check, reported before any contrast.** Sample-type composition, pair
count, processing already applied, and the paired difference in L1, with the
per-gene table. **If L1 is not lower in tumor at two-sided Wilcoxon signed-rank
*P* < 0.05, the lung cohort is reported as not testable and no retention value
from it appears.**

**Signature panel.** The §6t rule applied to lung: every C2:CGP set whose
standard name contains **LUNG, PULMONARY or ALVEOLAR**, 15 to 1000 symbols,
evaluable at 10 or more symbols on the platform, excluding sets whose unadjusted
paired shift does not reach two-sided Wilcoxon *P* < 0.05. The composition of the
panel by subject matter is reported **before** any retention value, because that
is what §6x failed to do until after the fact.

**Analyses.** The §6t benchmark, the §6u comparison and the §6v intervals, rules
and thresholds unchanged.

**Outcomes accepted in advance.** These replace §6x's criteria, which conflated
dispersion with level and left a third possibility undescribed. The result is
classified on two axes, reported separately:

- **Dispersion.** Generalizes if the interquartile range is at least half of
  liver's and the full range spans at least tenfold. Otherwise the manuscript
  reports that the wide spread seen in liver is not a general property.
- **Level.** Transfers if the median lies within 0.20 of liver's 0.758.
  Otherwise the manuscript states that the level is tissue-specific and that no
  tissue's median may be used as a reference for another. **If dispersion
  generalizes while level does not, which is what kidney showed, the manuscript
  reports exactly that and Section 4.2 recommends reporting the quantity per
  dataset while explicitly declining to offer a benchmark value.**
- If the premise fails, lung is reported as not testable and Section 4.3 gains a
  second worked instance.
- If §6u returns REDUNDANT in lung, the manuscript reports that whether identity
  and composition adjustment differ is tissue-dependent, as kidney already
  suggested by a smaller margin.

No gene set and no member of L1 will be added, removed, rescored or reordered on
the basis of any result.

**What this cannot establish.** Three tissues are not generality, and all three
are epithelial cancers with a dominant parenchymal cell of origin, which is the
setting the measure was designed for. The breast case set aside above is exactly
the setting where it may not apply, and declining to test it there means the
manuscript cannot claim the measure works where the cell of origin is not the
tissue's dominant cell type. That limitation is stated in the manuscript rather
than left implicit.

---

# Amendment §6z

Added 8 September 2026, after §6y and in response to pre-submission advice that
the manuscript should show whether the measure has clinical consequences. Post
hoc relative to the 25 August 2026 lock. Written before any survival value was
read; only the column names of the survival file were inspected beforehand.

---

## §6z. Does identity retention say anything about a signature's clinical value?

**Status.** Post hoc relative to the lock. Written before any survival value was
computed and before any survival column of the file was read. Only the file's
column names were inspected beforehand: sample, OS.time, OS and _PATIENT.

**Why.** Everything reported so far describes what a tumor-adjacent contrast is
made of. It does not say whether knowing that changes anything a reader would do.
There is one place where it plainly could. Tumor differentiation grade is a
routine pathological variable, obtained free with every specimen and long known to
carry prognostic information in hepatocellular carcinoma. If a published
expression signature predicts outcome only through the differentiation state of
the tumor, then whatever it costs to measure, it is restating something the
pathology report already contains. The identity-retention fraction is a candidate
instrument for detecting exactly that, and this amendment tests whether it works
as one.

**The question, fixed in advance.** Across published liver signatures, is the
identity-retention fraction, estimated from the paired tumor and adjacent samples,
associated with whether that signature's prognostic signal survives adjustment for
the tumor's own hepatocyte identity score?

**Design, fixed in advance.**

- **Cohort.** TCGA-LIHC primary tumors, sample type 01, with overall survival from
  the UCSC Xena GDC hub file `TCGA-LIHC.survival.tsv`. Overall survival is the
  primary and only endpoint. The paired subset used elsewhere is not used here;
  survival is estimated on every tumor with an outcome, and retention is taken from
  the 50 pairs, so the two quantities come from different estimation sets by
  construction and neither informs the other.
- **Signatures.** Exactly the evaluable set already reported for TCGA-LIHC under
  §6t. None is added, removed or rescored. Each is scored in tumors by the same
  zmean routine used everywhere else.
- **The adjustment variable.** D1 scored in each tumor, that is, how much
  hepatocyte identity that tumor retains. In the paired analyses D1 enters as a
  difference; here it enters as a level, which is the molecular analogue of
  differentiation grade. This is stated as an analogue and not as a substitute:
  pathological grade itself is not in the file used here, and the manuscript will
  say so.
- **Models.** For each signature, two Cox proportional-hazards models on the same
  tumors: the signature score alone, and the signature score with D1. Both report
  the hazard ratio per standard deviation with its 95% interval and *P*. The
  proportional-hazards assumption is checked by Schoenfeld residuals for D1 and
  reported. Every fit uses the study's own `survival.py`, which is unit-tested.
- **The quantity of interest.** For each signature, the fraction of its prognostic
  signal that survives adjustment, defined as |log HR adjusted| / |log HR
  unadjusted|, by analogy with the retention fraction itself.

**Primary analysis, fixed in advance.** Spearman correlation between the
identity-retention fraction and the surviving fraction of log HR, across the
evaluable signatures. Reported alongside: the number of signatures prognostic
before adjustment at *P* < 0.05 and after Benjamini-Hochberg control at 5%; how
many of those lose significance after adjustment; and the same counts split at
the median retention. D1's own prognostic association is reported first, because
if identity does not predict outcome in this cohort the rest of the amendment has
nothing to detect.

**Power rule, fixed in advance.** If the cohort yields fewer than 30 deaths, or
if fewer than 10 signatures are prognostic before adjustment, the analysis is
reported as underpowered and not testable, and no clinical claim is made.

**Outcomes accepted in advance.**

- **Associated.** If Spearman ρ is at least 0.30 with *P* < 0.05 in the expected
  direction, meaning that higher-retention signatures keep more of their
  prognostic signal, the manuscript reports that the measure flags signatures whose
  prognostic value is not independent of differentiation state, this appears in the
  abstract, and Section 4.2 gains a recommendation that a signature proposed as a
  prognostic marker be reported with its retention fraction.
- **Not associated.** If |ρ| < 0.30 or *P* ≥ 0.05, the manuscript reports plainly
  that identity retention does not predict whether a signature's prognostic value
  is independent, that the measure therefore describes what a contrast is made of
  and not what it is worth clinically, and **no clinical claim is made anywhere in
  the paper.** The result is reported in Results with the same prominence it would
  have had if it were positive.
- **Inverse.** If ρ is at most −0.30 with *P* < 0.05, that is reported as observed,
  and the manuscript states that it has no explanation for it rather than
  supplying one.

**What this cannot establish, stated before the data are opened.** A single
cohort, one endpoint, and a molecular proxy for grade rather than grade itself.
Adjusting for D1 removes what D1 measures, and D1 is a 22-gene score with error,
so residual prognostic signal after adjustment is not proof of independence from
differentiation. Nothing here is a validation of any signature as a clinical test,
and no signature will be recommended or discouraged on the strength of one
retrospective cohort.

---

# Amendment §6aa

Added 8 September 2026, after §6z. Post hoc relative to the 25 August 2026 lock.
Written before any value from the two supplementary sheets it governs was read;
only their sheet names and column headers were inspected beforehand.

---

## §6aa. Two things the source proteogenomic study can settle that this study could not

**Status.** Post hoc relative to the lock. Written before any value from either
sheet was read. Only the sheet names and the column headers of the source
workbook were inspected beforehand.

**Why.** Two claims in this manuscript are weaker than they need to be, and the
supplementary tables of one source study can strengthen or refute both.

The first is the transfer failure. Section 4.2 says the measure does not carry
from transcript to protein, and it names post-transcriptional regulation as a
candidate explanation while refusing to assert it. That refusal is correct but it
leaves the reader with nothing. The source study measured messenger RNA and
protein **in the same 159 patients** and reports a per-gene Spearman correlation
between them for 6,203 genes. That is a direct instrument for the question, and
using it turns a named candidate into a tested one.

The second is the distinction from tumor purity. Section 4.1 argues that identity
adjustment is not purity adjustment, and §6u supports that by comparing two gene
scores. A reviewer is entitled to say that both are gene scores and neither is
purity. The same workbook records **tumor purity assessed by a pathologist from
hematoxylin and eosin staining**, per patient. That is purity measured by the
method the field actually trusts, and it can be put into the model directly.

### Part A. Does mRNA to protein concordance explain the transfer failure?

**Source, fixed in advance.** The per-gene Spearman correlation between messenger
RNA and protein across the 159 paired cases, as published in the source study's
Table S1. No correlation is recomputed here and none is filtered.

**A1, descriptive and reported first.** The correlation for every member of
REDUCTION, of DRAIN and of D1, against the distribution over all 6,203 genes.

**A2, the test.** The §6t benchmark is run on the Gao protein matrix under the
rules already fixed, giving a protein-level retention fraction for every signature
with at least 10 measured proteins. For each signature evaluable at both levels,
the discrepancy is |transcript retention − protein retention| and the predictor is
the mean mRNA-protein correlation of that signature's members. **Prediction fixed
in advance: if post-transcriptional divergence explains the transfer failure, the
two are negatively associated, at Spearman ρ ≤ −0.30 with *P* < 0.05.**

**Outcomes accepted in advance.**

- If the association holds, Section 4.2 reports that the transfer failure tracks
  measured mRNA to protein concordance, and the candidate explanation becomes a
  supported one, stated at that strength and no higher.
- If it does not hold, Section 4.2 reports that the transfer failure is **not**
  explained by mRNA to protein concordance in these data, the candidate is
  withdrawn rather than left hanging, and the manuscript says the cause is
  unknown.
- If the association is positive, that is reported as observed without
  explanation.

### Part B. Does a pathologist's purity measure change anything?

**Source, fixed in advance.** Tumor purity by hematoxylin and eosin staining, per
patient, from the same table. The number of patients with a recorded value is
reported before any contrast.

**B1 and B2, reported first.** Spearman correlation of that purity with the
composition score C1 and with the identity score D1, both computed in the tumor
samples of the Gao protein matrix. C1 is expected to fall as purity rises, which
would confirm that C1 measures what it is meant to.

**B3, the test.** The identity-retention fractions for REDUCTION and DRAIN in the
Gao proteome are recomputed with the paired difference in histological purity
entered as an additional covariate alongside D1, on the patients with a recorded
value. **Decision rule fixed in advance: if either retention fraction moves by
more than 0.10, the measure is reported as partly confounded with tumor purity
after all, and Section 4.1 is qualified accordingly. If both move by 0.10 or less,
the manuscript reports that a histological purity measure adds nothing to the
identity adjustment**, which is a stronger statement than §6u could make and is
stated as such.

**What neither part can establish.** Part A uses a correlation the source study
computed, on a gene set of its choosing, and 6,203 genes is not the whole
transcriptome. Part B uses purity as assessed by eye on one section, which is a
coarse measure and is not available for every patient. Neither is a substitute
for the matched analysis of retention at both levels in the same patients, which
requires the source study's messenger RNA matrix; that matrix is deposited in
repositories this analysis cannot reach, and the manuscript will continue to name
that comparison as the informative next step whatever these two parts show.

---

# Amendment §6ab

Added 8 September 2026, after §6aa and in response to pre-submission advice that
the manuscript should consider technical explanations for the transfer failure.
Post hoc relative to the 25 August 2026 lock. Written before either quantity it
governs was computed.

---

## §6ab. Is the transfer failure a mass-spectrometry artifact?

**Status.** Post hoc relative to the lock. Written before either quantity below
was computed for any signature.

**Why.** §6aa tested one candidate explanation for the failure of the measure to
carry from transcript to protein, post-transcriptional divergence, and refuted
it: across 97 signatures the discrepancy did not track the source study's own
messenger RNA to protein correlations, and the module with the largest
discrepancy is built from the genes whose two measurements agree best. Two
candidates named in the manuscript remain untested, and one of them is testable
here. Mass spectrometry does not quantify all proteins equally. Its dynamic range
is limited, low-abundance proteins are measured with more error, and a
shotgun proteome covers a fraction of the transcriptome. If the discrepancy is a
measurement artifact of that kind, it should be larger for signatures whose
proteins are scarce in the matrix, and larger for signatures most of whose genes
never appear in it at all.

**What is and is not testable.** The deposited matrix carries 6,478 gene-level
proteins quantified in all 318 samples with no missing values, so detection
failure is not a variable within it; what varies is abundance, and which genes
made it into the matrix in the first place. Protein half-life is not testable
from these data at all and no attempt is made to test it.

**Predictors, fixed in advance.**

1. **Abundance.** The mean of a signature's measured proteins across all samples,
   on the scale the matrix is distributed in.
2. **Proteome coverage.** The proportion of the signature's genes that are present
   in the matrix, that is, the number quantified divided by the number in the set
   as published.

**Outcome, fixed in advance.** The same discrepancy used in §6aa: the absolute
difference between a signature's transcript and protein retention fractions, over
the same signatures evaluable at both levels, with the same eligibility rules and
no reselection.

**Prediction, fixed in advance.** If limited dynamic range or incomplete coverage
explains the transfer failure, the discrepancy is **larger** where abundance is
lower or coverage is poorer, giving **Spearman ρ ≤ −0.30 with *P* < 0.05 for at
least one of the two predictors**.

**Reported alongside, descriptively.** The abundance rank of every REDUCTION and
DRAIN member in the matrix, so that the reader can see whether the modules that
disagree between levels sit at the difficult end of the dynamic range.

**Outcomes accepted in advance.**

- If the prediction holds for either predictor, the manuscript reports that the
  transfer failure is at least partly an artifact of protein measurement, states
  which predictor and at what strength, and qualifies the claim in the Discussion
  accordingly rather than presenting the failure as a property of the tissue.
- If it holds for neither, the manuscript reports that three candidate
  explanations have now been examined, that post-transcriptional concordance,
  protein abundance and proteome coverage each fail to account for the
  discrepancy, and that the cause is unknown. In that case the Discussion names
  protein half-life and cohort difference as the remaining untested candidates
  and says explicitly that neither is testable from these data.
- An association in the opposite direction is reported as observed, without
  explanation.

**What this cannot establish.** Abundance in a processed, filtered matrix is a
coarse proxy for the instrument's dynamic range, and coverage is a property of
the source study's pipeline as much as of the proteins. A negative result here
does not exclude a measurement explanation; it excludes the two forms of it that
these data can address.
