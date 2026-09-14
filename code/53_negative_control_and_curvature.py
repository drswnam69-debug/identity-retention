#!/usr/bin/env python3
"""53_negative_control_and_curvature.py -- two post hoc checks a referee will ask for.

Written 9 September 2026, after the results it concerns already existed. Post
hoc relative to the 25 August 2026 lock and to every amendment. Nothing here is
pre-registered; it is recorded as §6ae.

A. MATCHED-RANDOM NEGATIVE CONTROL.
   The measure is named for hepatocyte identity, and the obvious objection is
   that D1 is simply the dominant tumor-versus-normal axis under another name.
   The test: build covariates from random gene sets matched to D1 on size, on
   mean expression, and on the magnitude of their own paired tumor fall, then
   recompute all 119 retention fractions against each. If the distribution
   obtained from a matched random covariate looks like the distribution
   obtained from D1, nothing about identity is being measured.

B. CURVATURE SENSITIVITY.
   §6q found a quadratic term in dD1 warranted in this cohort (P = 6.0e-4) but
   tested it only on the composite index. The intercept is an extrapolation to
   a sparsely observed corner, so the whole 119-signature distribution is
   recomputed with a quadratic term and the summaries compared.

Both read the cohort's own rsi.tsv and expression matrix and the archived
benchmark, and neither changes a reported value on its own.
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np
import pandas as pd

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_cands = [_here, _os.path.join(_here, "rsi", "code"), _os.path.join(_here, "code"),
          _os.path.join(_os.path.dirname(_here), "code")]
if _os.environ.get("IR_ROOT"):
    _cands.insert(0, _os.path.join(_os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        if _c not in _sys.path:
            _sys.path.insert(0, _c)
        break
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, \
    DATA as IR_DATA, CODE as IR_CODE, FIGURES as IR_FIGURES, DOCS as IR_DOCS
sys.path.insert(0, f"{IR_ROOT}/code")

RES = f"{IR_RESULTS}"
GS = f"{IR_GENESETS}"
SEED = 20260909
N_DRAWS = 200

D1 = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2",
      "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1", "HNF4A", "HNF1A",
      "FOXA1", "FOXA2", "NR1H4"]
C1 = ["PTPRC", "CD53", "LAPTM5", "CD3E", "CD2", "CD68", "AIF1", "ITGAM", "LCP1",
      "PECAM1", "VWF", "CDH5", "ENG", "CLEC4G", "COL1A1", "COL1A2", "COL3A1",
      "DCN", "LUM", "PDGFRB", "ACTA2"]


def load():
    expr = pd.read_csv(f"{RES}/../data/GSE14520_symbols.tsv.gz", sep="\t",
                       index_col=0)
    ph = pd.read_csv(f"{RES}/GSE14520/phenotype.tsv", sep="\t", index_col=0)
    t = ph[ph.tissue == "HCC tumor"]
    a = ph[ph.tissue != "HCC tumor"]
    pt = dict(zip(t.patient_id, t.index))
    pn = dict(zip(a.patient_id, a.index))
    ids = sorted(set(pt) & set(pn))
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1), axis=0)
    z = z[expr.std(axis=1) > 0]
    it = [pt[i] for i in ids]
    inn = [pn[i] for i in ids]
    dz = z[it].to_numpy() - z[inn].to_numpy()          # genes x pairs
    return expr, z, pd.Index(z.index), dz, len(ids)


def score(dz, index, genes):
    idx = [index.get_loc(g) for g in genes if g in index]
    return dz[idx].mean(axis=0) if idx else None


def retention(dy, X):
    b, *_ = np.linalg.lstsq(X, dy, rcond=None)
    return abs(float(b[0])) / abs(float(dy.mean())), float(b[0])


def gmt(path):
    out = {}
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            out[f[0]] = f[2:]
    return out


def main():
    expr, z, index, dz, n = load()
    bench = json.load(open(f"{RES}/SIGNATURE_BENCHMARK_GSE14520.json"))
    sets = gmt(f"{GS}/liver_cgp.gmt")
    ok = [s for s in bench["signatures"] if s["status"] == "ok"]
    print(f"GSE14520, {n} pairs, {len(ok)} evaluable signatures")

    d_d1 = score(dz, index, D1)
    d_c1 = score(dz, index, C1)
    ones = np.ones(n)
    X_joint = np.column_stack([ones, d_d1, d_c1])

    # ---- reproduce the archived joint retention as a control ---------------
    rep_err = []
    base = {}
    for s in ok:
        g = [x for x in sets.get(s["name"], []) if x in index]
        if len(g) < 10:
            continue
        dy = score(dz, index, g)
        r, _ = retention(dy, X_joint)
        base[s["name"]] = r
        rep_err.append(abs(r - s["retention_joint"]))
    print(f"  archived joint retention reproduced for {len(base)} sets; "
          f"max |difference| {max(rep_err):.4f}, median {np.median(rep_err):.4f}")

    # ---- A. matched random covariates --------------------------------------
    mean_expr = expr.mean(axis=1)
    d1_present = [g for g in D1 if g in index]
    d1_fall = float(abs(np.mean([dz[index.get_loc(g)].mean() for g in d1_present])))
    d1_meanexpr = float(np.mean([mean_expr[g] for g in d1_present]))
    panel = set(D1) | set(C1)
    pool = [g for g in index if g not in panel]
    pool_fall = np.array([dz[index.get_loc(g)].mean() for g in pool])
    pool_expr = np.array([mean_expr[g] for g in pool])
    # candidates that fall in tumor as D1 does, at a similar magnitude and level
    cand = [g for g, f, e in zip(pool, pool_fall, pool_expr)
            if f < 0 and abs(abs(f) - d1_fall) < 0.35
            and abs(e - d1_meanexpr) < 1.5]
    print(f"  matched pool: {len(cand)} genes fall in tumor within 0.35 z of "
          f"D1's mean fall and within 1.5 of its mean expression")
    rng = np.random.default_rng(SEED)
    meds, below, corr = [], [], []
    # Figure 12 plots one point per draw. Recording the draws here, inside the
    # run that writes the summary, is what keeps the figure and the archive
    # from describing two different realizations of the same pools.
    draws = []
    for _ in range(N_DRAWS):
        pick = list(rng.choice(cand, size=len(d1_present), replace=False))
        d_r = score(dz, index, pick)
        Xr = np.column_stack([ones, d_r, d_c1])
        vals = []
        for s in ok:
            g = [x for x in sets.get(s["name"], []) if x in index]
            if len(g) < 10:
                continue
            r, _ = retention(score(dz, index, g), Xr)
            vals.append(r)
        vals = np.array(vals)
        meds.append(float(np.median(vals)))
        below.append(int((vals < 0.5).sum()))
        draws.append(("matched", round(float(d_r.mean()), 6),
                      round(float(np.median(vals)), 6), int((vals < 0.5).sum())))
        common = [base[s["name"]] for s in ok
                  if s["name"] in base and len([x for x in sets.get(s["name"], [])
                                                if x in index]) >= 10]
        corr.append(float(np.corrcoef(vals, common)[0, 1]))
    real_med = float(np.median(list(base.values())))
    real_below = int(sum(1 for v in base.values() if v < 0.5))
    print("\n  A. matched-random negative control, "
          f"{N_DRAWS} draws of {len(d1_present)} genes")
    print(f"     D1            median {real_med:.3f}   below 50%: {real_below}")
    print(f"     random matched median {np.median(meds):.3f} "
          f"(2.5-97.5th {np.percentile(meds,2.5):.3f} to "
          f"{np.percentile(meds,97.5):.3f})   "
          f"below 50%: median {int(np.median(below))} "
          f"({int(np.percentile(below,2.5))} to {int(np.percentile(below,97.5))})")
    print(f"     Pearson r between the D1 and the random retention vectors: "
          f"median {np.median(corr):.3f}")
    p_med = float((np.array(meds) <= real_med).mean())
    print(f"     fraction of random covariates giving a median at or below "
          f"D1's: {p_med:.3f}")

    # C1 on its own, the comparator Figure 12b marks beside D1
    X_c1 = np.column_stack([ones, d_c1])
    c1_vals = []
    for s in ok:
        g = [x for x in sets.get(s["name"], []) if x in index]
        if len(g) < 10:
            continue
        c1_vals.append(retention(score(dz, index, g), X_c1)[0])
    c1_vals = np.array(c1_vals)
    c1_med = float(np.median(c1_vals))
    c1_below = int((c1_vals < 0.5).sum())

    # ---- A2. what governs what a covariate removes? -------------------------
    # If a matched falling covariate behaves like D1, the next question is
    # whether the gene list matters at all, or only the size of the covariate's
    # own paired shift. Three further pools answer that.
    fall_all = np.array([dz[index.get_loc(g)].mean() for g in pool])
    arms = {"unrestricted": pool,
            "no_tumor_shift": [g for g, f in zip(pool, fall_all) if abs(f) < 0.05],
            "rises_in_tumor": [g for g, f in zip(pool, fall_all) if f > 0.35]}
    arm_out = {}
    for label, cands in arms.items():
        m2, b2 = [], []
        for _ in range(100):
            pick = list(rng.choice(cands, size=len(d1_present), replace=False))
            Xr = np.column_stack([ones, score(dz, index, pick), d_c1])
            vals = []
            for s in ok:
                g = [x for x in sets.get(s["name"], []) if x in index]
                if len(g) < 10:
                    continue
                vals.append(retention(score(dz, index, g), Xr)[0])
            vals = np.array(vals)
            m2.append(float(np.median(vals))); b2.append(int((vals < 0.5).sum()))
            dr = score(dz, index, pick)
            draws.append((label, round(float(dr.mean()), 6),
                          round(float(np.median(vals)), 6), int((vals < 0.5).sum())))
        arm_out[label] = {"pool_size": len(cands),
                          "median_of_medians": round(float(np.median(m2)), 4),
                          "median_ci95": [round(float(np.percentile(m2, 2.5)), 4),
                                          round(float(np.percentile(m2, 97.5)), 4)],
                          "n_below_50_median": int(np.median(b2))}
        print(f"     {label:16s} pool {len(cands):5d}  median "
              f"{np.median(m2):.3f}  below 50%: {int(np.median(b2))}")

    # ---- B. curvature -------------------------------------------------------
    X_q = np.column_stack([ones, d_d1, d_c1, d_d1 ** 2])
    lin, quad = [], []
    for s in ok:
        g = [x for x in sets.get(s["name"], []) if x in index]
        if len(g) < 10:
            continue
        dy = score(dz, index, g)
        lin.append(retention(dy, X_joint)[0])
        quad.append(retention(dy, X_q)[0])
    lin, quad = np.array(lin), np.array(quad)
    print("\n  B. curvature sensitivity, quadratic term in dD1")
    print(f"     linear    median {np.median(lin):.3f}  IQR "
          f"{np.percentile(lin,25):.3f} to {np.percentile(lin,75):.3f}  "
          f"below 50%: {int((lin<0.5).sum())} of {len(lin)}")
    print(f"     quadratic median {np.median(quad):.3f}  IQR "
          f"{np.percentile(quad,25):.3f} to {np.percentile(quad,75):.3f}  "
          f"below 50%: {int((quad<0.5).sum())} of {len(quad)}")
    print(f"     Spearman between them: "
          f"{pd.Series(lin).corr(pd.Series(quad), method='spearman'):.3f}")

    out = {"plan": "PREREGISTRATION 6ae, post hoc; nothing here is pre-registered",
           "cohort": "GSE14520", "n_pairs": n, "n_sets": len(lin),
           "reproduction_of_archive": {"n": len(base),
                                       "max_abs_diff": round(max(rep_err), 4),
                                       "median_abs_diff": round(float(np.median(rep_err)), 4)},
           "negative_control": {
               "n_draws": N_DRAWS, "n_genes_per_draw": len(d1_present),
               "matched_pool_size": len(cand),
               "match_rule": "falls in tumor; |mean fall - D1 mean fall| < 0.35 z; "
                             "|mean expression - D1 mean expression| < 1.5",
               "D1_median": round(real_med, 4), "D1_n_below_50": real_below,
               "random_median_of_medians": round(float(np.median(meds)), 4),
               "random_median_ci95": [round(float(np.percentile(meds, 2.5)), 4),
                                      round(float(np.percentile(meds, 97.5)), 4)],
               "random_n_below_50_median": int(np.median(below)),
               "random_n_below_50_ci95": [int(np.percentile(below, 2.5)),
                                          int(np.percentile(below, 97.5))],
               "median_r_with_D1_retention": round(float(np.median(corr)), 4),
               "fraction_random_at_or_below_D1_median": round(p_med, 4),
               "other_pools": arm_out},
           "curvature": {
               "linear_median": round(float(np.median(lin)), 4),
               "linear_iqr": [round(float(np.percentile(lin, 25)), 4),
                              round(float(np.percentile(lin, 75)), 4)],
               "linear_n_below_50": int((lin < 0.5).sum()),
               "quadratic_median": round(float(np.median(quad)), 4),
               "quadratic_iqr": [round(float(np.percentile(quad, 25)), 4),
                                 round(float(np.percentile(quad, 75)), 4)],
               "quadratic_n_below_50": int((quad < 0.5).sum()),
               "spearman_linear_quadratic": round(
                   float(pd.Series(lin).corr(pd.Series(quad), method="spearman")), 4)}}
    with open(f"{RES}/NEGATIVE_CONTROL_6ae.json", "w") as fh:
        json.dump(out, fh, indent=1)
    # the per-draw record Figure 12 is drawn from, and the two study covariates
    # it marks, both written by the run that wrote the summary above
    with open(f"{RES}/NEGATIVE_CONTROL_6ae_draws.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pool", "covariate_mean_paired_shift",
                    "median_retention", "n_below_50"])
        w.writerows(draws)
    with open(f"{RES}/NEGATIVE_CONTROL_6ae_reference.json", "w") as fh:
        json.dump({"D1": {"shift": round(float(d_d1.mean()), 4),
                          "median": round(real_med, 4),
                          "n_below_50": real_below},
                   "C1_alone": {"shift": round(float(d_c1.mean()), 4),
                                "median": round(c1_med, 4),
                                "n_below_50": c1_below}}, fh, indent=1)
    print("\nwrote results/NEGATIVE_CONTROL_6ae.json")


if __name__ == "__main__":
    main()
