#!/usr/bin/env python3
"""24_por_sensitivity.py -- PREREGISTRATION 6l.

POR does not donate electrons to mARC; it entered DRAIN on a false premise.
The locked index is NOT changed. This recomputes, as a sensitivity analysis
specified in 6l before it was run:

    DRAIN' = mean_z(MTARC1, MTARC2)
    RSI'   = 0.5*(mean_z SUPPLY + mean_z REDUCTION) - mean_z DRAIN'

and puts RSI' and DRAIN' through exactly the 6g / 6j machinery, unchanged.
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stats_lite as sl
import rsi_config as cfg


def _load(name, fn):
    sp = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


_da = _load("diffadj", "12_differentiation_adjust.py")
_ca = _load("compadj", "18_composition_adjust.py")
ols_ci, D1 = _da.ols_ci, _da.D1
C1 = _ca.C1

DRAIN_PRIME = ["MTARC1", "MTARC2"]
assert "POR" not in DRAIN_PRIME
assert set(DRAIN_PRIME) < set(cfg.MODULE_DRAIN), \
    "6l requires DRAIN' to be a strict subset of the locked DRAIN"


def zmean(mat, genes):
    g = [x for x in genes if x in mat.index]
    miss = sorted(set(genes) - set(g))
    z = mat.loc[g].sub(mat.loc[g].mean(axis=1), axis=0) \
                  .div(mat.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0), len(g), miss


def paired_intercept(dy, dx_cols, names):
    X = np.column_stack([np.ones(len(dy))] + list(dx_cols))
    return ols_ci(np.asarray(dy, float), X, ["intercept"] + names)


def _w(d):
    d = np.asarray(d, float)
    d = d[d != 0]
    n = len(d)
    r = sl.rankdata(np.abs(d))
    wp = r[d > 0].sum()
    mu = n * (n + 1) / 4.0
    _, cnt = np.unique(np.abs(d), return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24.0 - ((cnt ** 3 - cnt).sum()) / 48.0
    z = (wp - mu) / np.sqrt(var)
    return float(wp), float(z), sl.norm_two_sided(z), n


# the covariate panels live outside expr_log for GSE76427; identical sources
# to those 12_/18_ were run with, so nothing here is renormalized or re-derived
SOURCES = {
    "GSE76427": {"d1": "data/GSE76427_d1_panel.tsv",
                 "c1": "data/GSE76427_composition_panel.tsv"},
    "GSE14520": {"d1": None,   # D1 is inside expr_log for this cohort
                 "c1": "data/GSE14520_composition_panel.tsv"},
}


def run(cohort):
    e = pd.read_csv(f"results/{cohort}/expr_log.tsv.gz", sep="\t", index_col=0)
    ph = pd.read_csv(f"results/{cohort}/phenotype.tsv", sep="\t", index_col=0)
    rsi = pd.read_csv(f"results/{cohort}/rsi.tsv", sep="\t", index_col=0)
    keep = [c for c in e.columns if c in ph.index and c in rsi.index]
    e, ph = e[keep], ph.loc[keep]

    src = SOURCES[cohort]
    d1src = e if src["d1"] is None else \
        pd.read_csv(src["d1"], sep="\t", index_col=0)[keep]
    c1src = pd.read_csv(src["c1"], sep="\t", index_col=0)[keep]

    sup, _, _ = zmean(e, cfg.MODULE_SUPPLY)
    red, _, _ = zmean(e, cfg.MODULE_REDUCTION)
    drn, n_d, miss_d = zmean(e, cfg.MODULE_DRAIN)
    drp, n_p, miss_p = zmean(e, DRAIN_PRIME)
    d1, nd1, _ = zmean(d1src, D1)
    c1, nc1, _ = zmean(c1src, C1)
    assert nd1 == 22, f"D1 must be 22/22 as in 6i, got {nd1}"

    rsi_prime = 0.5 * (sup + red) - drp
    rsi_locked = 0.5 * (sup + red) - drn

    # integrity check: the re-derived locked index must match the written one
    r = rsi.loc[keep, "RSI"].astype(float)
    corr = float(np.corrcoef(rsi_locked.values, r.values)[0, 1])
    maxdev = float(np.abs(rsi_locked.values - r.values).max())

    is_t = ~ph["tissue"].astype(str).str.lower().str.contains(
        "adjacent|non-tumor|nontumor|normal")
    pid = ph["patient_id"]
    rows = []
    for p, g in pd.DataFrame({"p": pid, "t": is_t}).groupby("p"):
        ti = g.index[g["t"]]
        ni = g.index[~g["t"]]
        if len(ti) == 1 and len(ni) == 1:
            rows.append((ti[0], ni[0]))
    T = [a for a, _ in rows]
    N = [b for _, b in rows]

    out = {"plan": "PREREGISTRATION 6l", "cohort": cohort, "n_pairs": len(rows),
           "drain_locked_genes": n_d, "drain_prime_genes": n_p,
           "D1_genes": nd1, "C1_genes": nc1,
           "drain_locked_missing": miss_d, "drain_prime_missing": miss_p,
           "locked_index_reproduces": {"pearson_r": round(corr, 10),
                                       "max_abs_dev": round(maxdev, 10)}}

    dD1 = (d1[T].values - d1[N].values)
    dC1 = (c1[T].values - c1[N].values)
    for tag, series in (("RSI_locked", rsi_locked), ("RSI_prime", rsi_prime),
                        ("DRAIN_locked", drn), ("DRAIN_prime", drp),
                        ("REDUCTION", red)):
        d = series[T].values - series[N].values
        med = float(np.median(d))
        w, z, p, n = _w(d)
        blk = {"paired_median": round(med, 4),
               "paired_mean": round(float(d.mean()), 4),
               "wilcoxon_p": p}
        f1 = paired_intercept(d, [dD1], ["D1"])
        blk["adjusted_D1"] = {
            "intercept": f1["intercept"],
            "retained_fraction": round(f1["intercept"]["beta"] / med, 4)}
        fj = paired_intercept(d, [dD1, dC1], ["D1", "C1"])
        blk["adjusted_joint"] = {
            "intercept": fj["intercept"],
            "retained_fraction": round(fj["intercept"]["beta"] / med, 4)}
        out[tag] = blk
    return out


if __name__ == "__main__":
    res = {}
    for c in ("GSE76427", "GSE14520"):
        res[c] = run(c)
        o = res[c]
        print(f"\n=== {c} — 6l POR sensitivity ({o['n_pairs']} pairs) ===")
        print(f"  locked index reproduces: r={o['locked_index_reproduces']['pearson_r']:.10f} "
              f"max|dev|={o['locked_index_reproduces']['max_abs_dev']:.2e}")
        print(f"  DRAIN genes present: locked {o['drain_locked_genes']}/3, "
              f"DRAIN' {o['drain_prime_genes']}/2  missing {o['drain_locked_missing']}")
        for tag in ("RSI_locked", "RSI_prime", "DRAIN_locked", "DRAIN_prime",
                    "REDUCTION"):
            b = o[tag]
            print(f"  {tag:12s} median {b['paired_median']:+.3f} "
                  f"(P={b['wilcoxon_p']:.2e})"
                  f"  D1 {b['adjusted_D1']['intercept']['beta']:+.3f}"
                  f" [{b['adjusted_D1']['retained_fraction']*100:5.1f}%]"
                  f"  joint {b['adjusted_joint']['intercept']['beta']:+.3f}"
                  f" [{b['adjusted_joint']['retained_fraction']*100:5.1f}%]"
                  f" P={b['adjusted_joint']['intercept']['p']:.2e}")
    json.dump(res, open("results/POR_sensitivity_6l.json", "w"), indent=1)
    print("\nwrote results/POR_sensitivity_6l.json")
