# retentionfrac

How much of a paired tumor-adjacent expression contrast is not dedifferentiation.

    pip install -e .
    from retentionfrac import identity_retention

## Why one function, and why it can refuse

The measure is one regression. What is easy to get wrong is not the arithmetic,
it is the assumption underneath it: that the tissue's differentiated character
actually falls in the tumor samples you have. In the study this comes from, one
publicly deposited proteome showed albumin *higher* in tumor than in the same
patient's adjacent liver. Every retention value computable from that matrix was
meaningless, and nothing about the accession, the publication or the file format
said so.

So `identity_retention` checks the premise first and returns `retention=None`
with the reason in `notes` when it fails. It does the same when the unadjusted
shift is not reliably different from zero, because a ratio with no denominator is
not a small number, it is not a number. It also warns when the bootstrap interval
is wider than 0.50, which in practice means the cohort is too small to place the
signature against any threshold. In the source study only a 213-pair cohort was
large enough; 50 to 72 pairs were not.

## Usage

```python
import pandas as pd
from retentionfrac import identity_retention

expr = pd.read_csv("symbols_by_sample.tsv", sep="\t", index_col=0)  # log scale
tumors   = ["T01", "T02", "T03", ...]     # patient-matched, same order
adjacent = ["N01", "N02", "N03", ...]

D1 = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG",
      "F2", "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1",
      "HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4"]
C1 = [...]   # optional non-parenchymal covariate

r = identity_retention(expr, tumors, adjacent, my_signature, D1, C1)
print(r)
if r.retention is None:
    print("not reportable:", r.notes)
```

## Building an identity covariate for a tissue that is not liver

Three rules, from the source study:

1. genes marking the **differentiated** cell the tumor arose from, not the tissue
   generally;
2. **no** enzyme belonging to the process you are measuring, so the covariate
   cannot absorb your signal;
3. **disjoint** from the signature and from any composition covariate.

Fix the list before you look at any result. In the source study each tissue's
covariate was written out gene by gene in a timestamped amendment before its
expression matrix was opened, and the study reports what happened when a
pre-registered prediction then failed.

## Composition adjustment is a different operation

Adjusting for non-parenchymal content is not this. In one liver cohort it removed
essentially nothing (median 1.03 of the shift left standing) where identity
adjustment left 0.70; in lung it removed *more* than identity adjustment did.
Pass `composition_genes` to include it alongside, not instead.
