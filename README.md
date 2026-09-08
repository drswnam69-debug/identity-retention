# identity-retention

Code, gene sets, result files and the locked protocol for:

> **How much of a tumor expression signature is dedifferentiation? A pre-registered
> measure in three tissues and its failure to transfer between transcriptome and
> proteome.** Soon Woo Nam. *Preprint / under review.*

The repository also contains **`retentionfrac`**, a one-function Python package
that computes the quantity the article defines.

Pre-registration locked **25 August 2026**, SHA-256
`7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`,
with 28 timestamped amendments (`docs/PREREGISTRATION_with_amendments.md`).
Every amendment was written before the computation it governs, and several
record predictions that failed.

---

## What the measure is

For patient-matched tumor and adjacent tissue, regress the paired difference in
a signature score on the paired difference in a tissue-identity score built from
disjoint genes:

```
Delta signature = alpha + beta * Delta identity + epsilon
```

The **identity-retention fraction** is `|alpha|` divided by the unadjusted paired
shift: the share of the reported contrast that survives once the tumor's loss of
differentiated character is taken out. A value near 1 means the signature is
saying something beyond dedifferentiation. A value near 0 means it is largely
restating it.

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
docs/          the locked protocol with all 28 amendments
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
| Gene sets | MSigDB human C2:CGP v2026.1.Hs | included in `genesets/` |

`code/01_prepare.py` and `code/37_tcga_prepare.py` take the downloaded files and
write the tidy matrices the rest of the pipeline reads. Place downloads under
`data/` and run `run_all.sh`, or run the numbered scripts in order.

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

MIT. See `LICENSE`.
