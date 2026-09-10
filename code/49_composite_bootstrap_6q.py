#!/usr/bin/env python3
"""49_composite_bootstrap_6q.py -- archive the §6q-C composite bootstrap.

Table 2 quotes four bootstrap intervals on the composite index's retained
fraction. They were computed during the adversarial verification that prompted
§6q and were written into the amendment table, but never into a results file,
so no automatic check could reach them and they could not be reproduced from
the archive. This script recomputes them under the specification §6q-C fixed:
a non-parametric pair bootstrap with 4000 resamples, re-estimating numerator
and denominator on each resample, for the composite index under the D1-alone
and the joint D1 + C1 models in both paired cohorts.

This changes no reported result. It fills an archive gap.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/claude/rsi/code")
import rsi_config as cfg  # noqa: E402

RES = "/home/claude/rsi/results"
TUMOR = {"GSE14520": "HCC tumor",
         "GSE76427": "primary hepatocellular carcinoma tumor"}
D1 = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2",
      "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1", "HNF4A", "HNF1A",
      "FOXA1", "FOXA2", "NR1H4"]
C1 = ["PTPRC", "CD53", "LAPTM5", "CD3E", "CD2", "CD68", "AIF1", "ITGAM", "LCP1",
      "PECAM1", "VWF", "CDH5", "ENG", "CLEC4G", "COL1A1", "COL1A2", "COL3A1",
      "DCN", "LUM", "PDGFRB", "ACTA2"]
B = 4000
SEED = 20260901


def zmean(expr, genes):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g]
    sub = sub[(sub.std(axis=1) > 0).values]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0)


def load(cohort):
    expr = pd.read_csv(f"{RES}/../data/{cohort}_symbols.tsv.gz",
                       sep="\t", index_col=0)
    rsi = pd.read_csv(f"{RES}/{cohort}/rsi.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(f"{RES}/{cohort}/phenotype.tsv", sep="\t", index_col=0)
    t = ph[ph.tissue == TUMOR[cohort]]
    a = ph[ph.tissue != TUMOR[cohort]]
    pt = dict(zip(t.patient_id, t.index))
    pn = dict(zip(a.patient_id, a.index))
    ids = sorted(set(pt) & set(pn))
    pair = lambda ser: np.array([ser[pt[i]] - ser[pn[i]] for i in ids], float)
    return (pair(rsi["RSI"]), pair(zmean(expr, D1)), pair(zmean(expr, C1)))


def fraction(dv, X, denom):
    """Retained fraction: |intercept| over |unadjusted paired shift|.

    §6l fixes the composite denominator as the paired median, which is what
    Table 2's composite rows use, so the bootstrap uses the same denominator
    it is an interval for."""
    b, *_ = np.linalg.lstsq(X, dv, rcond=None)
    return abs(float(b[0])) / abs(float(denom))


# The point estimates §6q recorded, so that reproduction can be checked rather
# than assumed. GSE14520 reproduces; GSE76427 does not, for the reason the
# manuscript already reports about that cohort: its archived expression values
# were assembled gene by gene and the covariate cannot be rebuilt exactly.
RECORDED = {"GSE76427": {"D1": (0.503, 0.301, 0.759),
                         "joint_D1_C1": (0.399, 0.095, 0.706)},
            "GSE14520": {"D1": (0.628, 0.450, 0.704),
                         "joint_D1_C1": (0.645, 0.452, 0.736)}}


def run(cohort):
    dv, dd1, dc1 = load(cohort)
    n = len(dv)
    rng = np.random.default_rng(SEED)
    out = {"n_pairs": n, "resamples": B, "seed": SEED,
           "denominator": "paired median, as §6l fixes for the composite"}
    for tag, cols in (("D1", (dd1,)), ("joint_D1_C1", (dd1, dc1))):
        X = np.column_stack([np.ones(n), *cols])
        point = fraction(dv, X, np.median(dv))
        boot = np.empty(B)
        for b in range(B):
            k = rng.integers(0, n, n)
            d = np.median(dv[k])
            if abs(d) < 1e-9:
                boot[b] = np.nan
                continue
            boot[b] = fraction(dv[k], np.column_stack([np.ones(n)] +
                                                      [c[k] for c in cols]),
                               d)
        lo, hi = np.nanpercentile(boot, [2.5, 97.5])
        rp, rlo, rhi = RECORDED[cohort][tag]
        ok = abs(point - rp) < 0.01 and abs(lo - rlo) < 0.02 and abs(hi - rhi) < 0.02
        out[tag] = {"retained": round(point, 4),
                    "ci95": [round(float(lo), 4), round(float(hi), 4)],
                    "excludes_50pct": bool(lo > 0.5 or hi < 0.5),
                    "recorded_in_6q": [rp, rlo, rhi],
                    "reproduces_6q": bool(ok)}
        print(f"  {cohort:9s} {tag:12s} {point:.3f} ({lo:.3f} to {hi:.3f})   "
              f"§6q recorded {rp:.3f} ({rlo:.3f} to {rhi:.3f})   "
              f"{'reproduces' if ok else 'DOES NOT reproduce'}")
    return out


def main():
    print("=== §6q-C composite bootstrap, 4000 resamples ===")
    res = {"plan": "PREREGISTRATION 6q-C, post hoc relative to the lock",
           "note": "archive gap fill; the values were in the amendment table only"}
    for c in ("GSE76427", "GSE14520"):
        res[c] = run(c)
    with open(os.path.join(RES, "COMPOSITE_BOOTSTRAP_6q.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("\nwrote results/COMPOSITE_BOOTSTRAP_6q.json")


if __name__ == "__main__":
    main()
