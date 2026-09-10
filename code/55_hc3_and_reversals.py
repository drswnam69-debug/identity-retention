#!/usr/bin/env python3
"""55_hc3_and_reversals.py -- archive two quantities the manuscript reports
that no result file carried.

1. §6q-D quotes HC3 robust standard errors on the intercept of ΔRSI on ΔD1 in
   both paired cohorts. Like the §6q-C bootstrap before it, that computation
   lived only in the amendment table, so nothing could check it. It is
   recomputed here. The amendment's intervals used a normal critical value;
   §6n corrected exactly that shortcut elsewhere in this pipeline, so the
   intervals are recomputed with the t critical value on the model's residual
   degrees of freedom and both forms are recorded.

2. The retention fraction is |intercept| over |unadjusted paired shift|, so it
   is positive whether adjustment shrinks an effect or reverses its sign. The
   manuscript now says so, and needs the counts. They are computed here from
   the archived benchmarks rather than recomputed, so they cannot drift from
   the values the benchmark reported.

This changes no reported result. It fills two archive gaps.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results")
DATA = os.path.join(os.path.dirname(RES), "data")
TUMOR = {"GSE14520": "HCC tumor",
         "GSE76427": "primary hepatocellular carcinoma tumor"}
D1 = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2",
      "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1", "HNF4A", "HNF1A",
      "FOXA1", "FOXA2", "NR1H4"]
# What §6q-D wrote down, so reproduction is checked rather than assumed.
RECORDED = {"GSE76427": {"se_ols": 0.137, "se_hc3": 0.160},
            "GSE14520": {"se_ols": 0.068, "se_hc3": 0.102}}


def zmean(expr, genes):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g]
    sub = sub[(sub.std(axis=1) > 0).values]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0)


def load(cohort):
    expr = pd.read_csv(f"{DATA}/{cohort}_symbols.tsv.gz", sep="\t", index_col=0)
    rsi = pd.read_csv(f"{RES}/{cohort}/rsi.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(f"{RES}/{cohort}/phenotype.tsv", sep="\t", index_col=0)
    t = ph[ph.tissue == TUMOR[cohort]]
    a = ph[ph.tissue != TUMOR[cohort]]
    pt = dict(zip(t.patient_id, t.index))
    pn = dict(zip(a.patient_id, a.index))
    ids = sorted(set(pt) & set(pn))
    pair = lambda ser: np.array([ser[pt[i]] - ser[pn[i]] for i in ids], float)
    return pair(rsi["RSI"]), pair(zmean(expr, D1))


def hc3(cohort):
    """OLS and HC3 standard errors on the intercept of Δcomposite on ΔD1."""
    y, x = load(cohort)
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    XtXi = np.linalg.inv(X.T @ X)
    beta = XtXi @ X.T @ y
    resid = y - X @ beta
    df = n - X.shape[1]
    se_ols = float(np.sqrt(np.diag(XtXi * (resid @ resid) / df))[0])
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)          # leverages
    omega = (resid / (1.0 - h)) ** 2                  # HC3 weights
    cov = XtXi @ (X.T * omega) @ X @ XtXi
    se_hc3 = float(np.sqrt(np.diag(cov))[0])
    tcrit = float(stats.t.ppf(0.975, df))
    b0 = float(beta[0])
    rec = RECORDED[cohort]
    return {
        "n_pairs": n,
        "intercept": round(b0, 4),
        "se_ols": round(se_ols, 4),
        "se_hc3": round(se_hc3, 4),
        "se_increase_pct": round(100.0 * (se_hc3 / se_ols - 1.0), 1),
        "ci95_ols_t": [round(b0 - tcrit * se_ols, 4), round(b0 + tcrit * se_ols, 4)],
        "ci95_hc3_t": [round(b0 - tcrit * se_hc3, 4), round(b0 + tcrit * se_hc3, 4)],
        "ci95_hc3_normal_as_6q_wrote_it": [round(b0 - 1.96 * se_hc3, 4),
                                           round(b0 + 1.96 * se_hc3, 4)],
        "t_critical_value": round(tcrit, 4),
        "residual_df": df,
        "hc3_excludes_zero": bool((b0 - tcrit * se_hc3) * (b0 + tcrit * se_hc3) > 0),
        "recorded_in_6q": rec,
        "reproduces_6q": bool(abs(se_ols - rec["se_ols"]) < 0.005
                              and abs(se_hc3 - rec["se_hc3"]) < 0.005),
    }


def reversals():
    """Signatures whose adjusted intercept has the opposite sign to their
    unadjusted shift, per cohort, from the archived benchmarks."""
    out = {}
    for coh in ("GSE14520", "GSE76427", "TCGA_LIHC", "TCGA_LUAD", "TCGA_KIRC"):
        p = f"{RES}/SIGNATURE_BENCHMARK_{coh}.json"
        if not os.path.exists(p):
            continue
        sg = [s for s in json.load(open(p, encoding="utf-8"))["signatures"]
              if s.get("retention_joint") is not None
              and isinstance(s.get("joint_intercept"), dict)]
        rev = [s for s in sg
               if s["joint_intercept"]["beta"] * s["unadjusted_mean_delta"] < 0]
        firm = [s for s in rev
                if s["joint_intercept"]["ci95"][0] * s["joint_intercept"]["ci95"][1] > 0]
        out[coh] = {
            "n_evaluable": len(sg),
            "n_sign_reversal": len(rev),
            "n_sign_reversal_intercept_ci_excludes_zero": len(firm),
            "names_with_ci_excluding_zero": sorted(s["name"] for s in firm),
        }
    return out


def main():
    res = {"plan": "archive gap fill for §6q-D and for the sign of the intercept",
           "note": "changes no reported result"}
    print("=== §6q-D, HC3 standard errors on the composite intercept ===")
    res["hc3"] = {}
    for c in ("GSE76427", "GSE14520"):
        d = hc3(c)
        res["hc3"][c] = d
        print(f"  {c:9s} intercept {d['intercept']:+.3f}  SE {d['se_ols']:.3f} -> "
              f"{d['se_hc3']:.3f} (+{d['se_increase_pct']:.0f}%)  "
              f"HC3 95% CI (t) {d['ci95_hc3_t'][0]:.3f} to {d['ci95_hc3_t'][1]:.3f}  "
              f"{'reproduces §6q' if d['reproduces_6q'] else 'DOES NOT reproduce §6q'}")
    print("\n=== signatures whose intercept reverses the unadjusted sign ===")
    res["sign_reversals"] = reversals()
    for coh, d in res["sign_reversals"].items():
        print(f"  {coh:10s} {d['n_sign_reversal']:3d} of {d['n_evaluable']:3d} "
              f"({d['n_sign_reversal_intercept_ci_excludes_zero']} with an "
              f"intercept interval excluding zero)")
    out = f"{RES}/HC3_AND_REVERSALS.json"
    json.dump(res, open(out, "w", encoding="utf-8"), indent=2)
    print("\nwritten", out)


if __name__ == "__main__":
    main()
