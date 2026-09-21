# identity-retention

Code, gene sets, result files and the locked protocol for:

> **How much of a tumor expression signature survives adjustment for the dominant
> tumor-adjacent axis? A pre-registered measure in three tissues and its failure to
> transfer between a transcriptome and a proteome.** Soon Woo Nam. *Under review.*

The repository also contains **`retentionfrac`**, a one-function Python package
that computes the quantity the article defines.

Pre-registration locked **25 August 2026**, SHA-256
`7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`,
with 31 timestamped amendments (`docs/PREREGISTRATION_with_amendments.md`).
Most were written before the computation they govern; five state in their own
text that they are not pre-registered in any sense, and several record
predictions that failed.

---


## Correction, 20 September 2026, revised 21 September 2026

The analyses in this archive up to v1.2.2 treat one of the two deposited liver
proteomes as violating a hepatocyte-identity premise and the other as satisfying
it, on the direction of a composite panel score. That reasoning is wrong, and the
conclusion drawn from it, that one matrix is not testable, is withdrawn.

Re-derived from the raw deposits, both proteomes behave the same way. The panel
splits in two. In Jiang the twelve members that rise are its ten secreted plasma
proteins and the two transcription factors quantified there, and the six that
fall are the five quantified metabolic enzymes and ASGR1, with no exception on
either side; in Gao the same partition holds for 17 of the 18 quantified members,
the exception being TTR. The secreted and metabolic classes do not overlap in
either dataset. The composite therefore averages to near zero, -0.139 in one and
+0.106 in the other, and the sign of a near-zero composite is arbitrary.

That one raw Jiang file admits three normalizations and two member sets. Across
those six pipelines the composite stays positive but its paired difference runs
from +0.002 to +0.090 and its Wilcoxon P from 0.025 to 0.77, while the class
separation is complete in every one of them. The split is the stable quantity;
the composite built from it is not.

Four points in the first version of this notice, dated 20 September, were stated
more strongly than the data support, and are corrected here.

  * Blood was excluded on erythrocyte proteins alone, which report cellular
    blood rather than plasma exudation. The exudation control has now been run
    on immunoglobulins and JCHAIN, plasma proteins the liver does not make:
    they average -0.285 in Gao and +0.074 in Jiang, against +0.408 and +0.514
    for the panel's secreted members. The conclusion stands; the earlier control
    did not establish it.
  * Abundance was excluded in both proteomes, and was described as pointing away
    from the split. Neither is right. The Gao supplementary table is centered
    within each protein, with 50.8% of its values below zero, so it carries no
    abundance information and the control reported for it in v1.2.3 is
    withdrawn. In Jiang abundance points the same way as the split and is far
    too small to produce it: proteins sitting where the secreted members sit
    move by +0.078 and those sitting where the metabolic members sit by -0.113,
    a gap of 0.19 against an observed class difference of 1.41.
  * Missingness was not examined. The Jiang matrix is 41% unquantified, and
    unevenly, at 36.2% of values in tumor samples against 45.5% in adjacent.
    Eleven panel members are quantified in both arms of all 124 pairs; among
    those eleven the split is unchanged and the composite is +0.294, three times
    its value over all 18, because dropping seven members changed the class
    balance. Detection among the panel members follows the direction of the
    shift rather than opposing it.
  * The rising group includes HNF4A and HNF1A, which are nuclear. No mechanism
    that appeals to secretion covers them, and no mechanism is proposed here.

A constructive finding follows from the same computation and is new in this
revision. The panel is not incoherent; three of its four classes hold together
in every matrix. The metabolic enzymes have a surviving share of 1.000 in all
five datasets and keep their direction at both measurement levels, ranking above
the 94th percentile among the 114 comparison sets in each. The secreted proteins
also hold together but reverse direction between levels. Only the union of the
two cancels. An analyst needing one hepatocyte identity covariate that means the
same thing in a liver transcriptome and a liver proteome can use the metabolic
subset; that subset has five or six members, was chosen by looking at these five
matrices, and was not tested on a held-out proteome.

The earlier work's failure to transfer the retention measure between transcript
and protein is consistent with this: a covariate whose members cancel barely
moves. That is offered as an explanation the earlier work did not consider, not
as proof that the three explanations it did test were the wrong ones.

`code/70_panel_coherence.py` writes every number to
`results/PANEL_COHERENCE.json`, including the surviving shares, the confidence
intervals, the per-member detection counts, the control member lists and the
panel's percentile among the 114 comparison sets, and `code/71_figures.py` draws
the figures in `figures/panel_coherence/` from that file alone.
`code/72_verify_manuscript.py` reads the manuscript and that JSON file and fails
if any number printed in the text cannot be reproduced from the JSON. Nothing else in
the archive has been altered: the superseded analyses are left in place rather
than removed, and this notice is what marks them.

## What the measure is

For patient-matched tumor and adjacent tissue, regress the paired difference in
a signature score on the paired difference in a tissue-identity score built from
disjoint genes:

```
Delta signature = alpha + beta * Delta identity + epsilon
```

The **identity-retention fraction** is `|alpha|` divided by the unadjusted paired
shift: the share of the reported contrast that survives once the paired change in
a pre-specified tissue-identity score is taken out. A value near 1 means the
signature carries information beyond that axis. A value near 0 means it is
largely restating it.

The measure is named for the covariate used to compute it and not for a
mechanism. Under amendment 6ae, a covariate drawn at random and matched to the
identity score on the size of its own paired tumor shift removes as much as the
identity score does, so what the quantity measures is alignment with a large
paired tumor-adjacent axis. Loss of tissue identity is the pre-specified and
interpretable instance of such an axis used here, not a demonstrated cause.

## Install and use

```bash
pip install -e tool
```

```python
import pandas as pd
from retentionfrac import identity_retention

expr = pd.read_csv("symbols_by_sample.tsv", sep="\t", index_col=0)  # log scale
tumors   = ["T01", "T02", "T03"]   # patient-matched, same order
adjacent = ["N01", "N02", "N03"]

r = identity_retention(expr, tumors, adjacent, my_signature, D1, C1)
print(r.retention, r.ci, r.notes)
```

`identity_retention` checks the identity premise before it computes anything and
returns `retention=None` with the reason when the premise fails, when the
unadjusted shift is not reliably different from zero, or when the bootstrap
interval is too wide to place the signature against a threshold. See
`tool/README.md` for the reasoning and for how to build an identity covariate for
a tissue that is not liver.

## Layout

```
tool/          the retentionfrac package, its README and pyproject.toml
code/          the full analysis pipeline, numbered in run order
genesets/      module definitions, the tissue covariates, and the complete
               MSigDB C2:CGP collection from which every enumeration in the
               article can be reproduced
results/       every derived result file, including
               SupplementaryTable_S8_enumeration.csv, which lists all 223
               enumerated gene sets and the inclusion status of each
figures/       the article's figures at 300 dpi
docs/          the locked protocol with all 31 amendments
SHA256SUMS.txt checksums for every file in this repository
```

Key scripts:

| Script | What it does |
|---|---|
| `code/rsi_config.py` | locked module definitions; not edited after the lock date |
| `code/k1_config.py`, `code/l1_config.py` | the kidney and lung identity covariates, fixed in §6x and §6y |
| `code/30_signature_benchmark.py` | enumerates the comparator panel and computes retention for every set |
| `code/34_covariate_comparison.py` | identity versus composition adjustment (§6u) |
| `code/35_retention_intervals.py` | bootstrap intervals on the ratio (§6v Part A) |
| `code/36_estimator_simulation.py` | recovery of a known retention fraction (§6v Part B) |
| `code/37_tcga_prepare.py` | TCGA matrix preparation and barcode pairing |
| `code/39_tissue_premise.py` | the premise check, for any cohort and covariate |
| `code/40_prognostic_value.py` | retention against prognostic value (§6z) |
| `code/41_gao_mrna_protein.py`, `code/42_gao_he_purity.py` | proteogenomic tests (§6aa) |
| `code/43_ms_artifact.py` | mass-spectrometry artifact test (§6ab) |

## Data

Expression matrices are not redistributed here. Every one is public:

| Source | Accession | Where |
|---|---|---|
| GEO series | GSE135251, GSE130970, GSE167523, GSE76427, GSE164760, GSE14520 | https://www.ncbi.nlm.nih.gov/geo/ |
| TCGA (STAR-TPM, GENCODE v36 probemap) | TCGA-LIHC, TCGA-LUAD, TCGA-KIRC | https://xenabrowser.net/datapages/ (GDC hub) |
| Proteome | Gao et al. 2019, Jiang et al. 2019 | the supplementary tables of the source articles |

Two small derived tables from the first of those, the per-gene messenger RNA to
protein correlations and the histological tumor purity, are in `data/` so that
the proteome analyses can be rerun without recovering them by hand. They are the
source article's material, not mine, and the MIT license below does not extend
to them.
| Gene sets | MSigDB human C2:CGP v2026.1.Hs | included in `genesets/` |

`code/01_prepare.py` and `code/37_tcga_prepare.py` take the downloaded files and
write the tidy matrices the rest of the pipeline reads. Place downloads under
`data/` and run `run_all.sh`, or run the numbered scripts in order.

## Paths, and what is not here

Every script resolves the archive from its own location through `code/paths.py`,
so the repository runs wherever it is unpacked and from whatever directory. Set
`IR_ROOT` to point somewhere else. Earlier releases named absolute paths on the
machine the analysis was run on, which meant a third of the code, including
every figure script, failed for anyone but the author.

The manuscript, its protocol and the cover letter are not deposited before
publication. Three scripts read them, `50_consistency_check.py`,
`51_check_the_checker.py` and `52_anchored_claims.py`; set `IR_DOCS` to the
directory holding them. Everything else runs from the archive alone.

## Testing

    python3 code/test_stats.py            # the hand-written statistics
    python3 code/test_survival.py         # the survival code
    python3 code/test_pipeline.py         # the pipeline, on a synthetic cohort
    python3 code/test_retentionfrac.py    # the package, against this archive

`test_retentionfrac.py` needs `data/GSE14520_symbols.tsv.gz`; the other three
run from the archive alone.

An earlier `29_cross_document_check.py` targeted a manuscript from a previous
submission and read files that are not in this archive. Its checks live in
`50_consistency_check.py` now, and it has been removed rather than shipped as
code that cannot run.

One known gap: the ferroptosis set file that `08_phase_e.py` reads for
hypothesis H6 (`--ferroptosis`) is not in this archive and is no longer on the
author's disk, so that one test cannot be re-run from the deposit. Its inputs,
outputs and the genes removed for overlap are recorded in
`results/phase_e_h6.json`. Everything else in the pipeline runs from what is
here plus the public matrices named above.

## Before you reuse a public matrix

The study's own hardest finding is procedural. One deposited proteome showed
albumin *higher* in tumor than in the same patient's adjacent liver. Nothing
about the accession, the publication or the file format said so, and every
retention value computable from that matrix was meaningless. Check that the
tissue's differentiated character actually falls in the tumor samples you have
before you compute a contrast from them. `identity_retention` does this for you
and refuses rather than returning a number.

## Citation

See `CITATION.cff`. Please cite both the software (Zenodo DOI) and the article.

## License

MIT for the code, and see `LICENSE`. Two sets of files in this archive are
redistributed under their own terms rather than under the MIT license: the
MSigDB collection in `genesets/`, and the two derived tables in `data/` that
come from the supplementary material of Gao et al. 2019.
