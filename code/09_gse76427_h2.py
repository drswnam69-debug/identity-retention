#!/usr/bin/env python3
"""09_gse76427_h2.py -- H2 replication in GSE76427 (PREREGISTRATION 6f).

Tumor vs adjacent non-tumor liver, unpaired and patient-paired. The index is
read as written by 02_compute_rsi.py; nothing here recomputes it.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_lite as sl


def hodges_lehmann(x, y):
    """Median of all pairwise differences x_i - y_j, with a Moses CI."""
    d = np.subtract.outer(np.asarray(x, float), np.asarray(y, float)).ravel()
    d.sort()
    est = float(np.median(d))
    n1, n2 = len(x), len(y)
    # normal approximation to the Wilcoxon rank-sum quantile
    mu = n1 * n2 / 2.0
    sd = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    k = int(np.floor(mu - 1.959964 * sd))
    if k < 0:
        return est, (float("nan"), float("nan"))
    lo = d[k]
    hi = d[len(d) - 1 - k]
    return est, (float(lo), float(hi))


def cliffs_delta(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def wilcoxon_signed_rank(d):
    """Two-sided Wilcoxon signed-rank, normal approx with tie correction."""
    d = np.asarray(d, float)
    d = d[d != 0]
    n = len(d)
    r = sl.rankdata(np.abs(d))
    w_plus = r[d > 0].sum()
    mu = n * (n + 1) / 4.0
    _, counts = np.unique(np.abs(d), return_counts=True)
    tie = (counts ** 3 - counts).sum()
    var = n * (n + 1) * (2 * n + 1) / 24.0 - tie / 48.0
    z = (w_plus - mu) / np.sqrt(var)
    return float(w_plus), float(z), sl.norm_two_sided(z), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="GSE76427")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()

    base = os.path.join(a.outdir, a.cohort)
    rsi = pd.read_csv(os.path.join(base, "rsi.tsv"), sep="\t", index_col=0)
    ph = pd.read_csv(os.path.join(base, "phenotype.tsv"), sep="\t", index_col=0)
    ph = ph.loc[rsi.index]

    is_tumor = ~ph["tissue"].astype(str).str.lower().str.contains("adjacent")
    t = rsi.loc[is_tumor.values, "RSI"].to_numpy(float)
    n = rsi.loc[~is_tumor.values, "RSI"].to_numpy(float)
    assert len(t) == 115 and len(n) == 52, (len(t), len(n))

    print(f"=== {a.cohort} — H2 replication (PREREG 6f) ===")
    print(f"  tumor n={len(t)}  adjacent n={len(n)}")
    print(f"  RSI mean  tumor {t.mean():+.3f}   adjacent {n.mean():+.3f}")

    u, z_u, p_u = sl.mannwhitney_u(t, n)
    hl, ci = hodges_lehmann(t, n)
    delta = cliffs_delta(t, n)
    print(f"  Mann-Whitney (two-sided): P = {p_u:.3g}")
    print(f"  Hodges-Lehmann shift (tumor - adjacent): "
          f"{hl:+.3f}  (95% CI {ci[0]:+.3f} to {ci[1]:+.3f})")
    print(f"  Cliff's delta: {delta:+.3f}")

    # patient-paired
    pid = ph["patient_id"].astype(str) if "patient_id" in ph.columns else None
    paired = None
    if pid is not None:
        df = pd.DataFrame({"pid": pid.values, "tumor": is_tumor.values,
                           "rsi": rsi["RSI"].to_numpy(float)},
                          index=rsi.index)
        tt = df[df.tumor].set_index("pid")["rsi"]
        nn = df[~df.tumor].set_index("pid")["rsi"]
        tt = tt[~tt.index.duplicated()]
        nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        d = (tt.loc[common] - nn.loc[common]).to_numpy(float)
        w, z, p_w, n_used = wilcoxon_signed_rank(d)
        print(f"  patient-matched pairs: {len(common)} "
              f"({n_used} non-zero differences)")
        print(f"  Wilcoxon signed-rank (two-sided): P = {p_w:.3g}   "
              f"median paired difference {np.median(d):+.3f}")
        paired = {"n_pairs": len(common), "n_nonzero": n_used,
                  "median_diff": round(float(np.median(d)), 4),
                  "W": w, "z": round(z, 4), "p": p_w}

    # Registered direction: PREREGISTRATION section 1 -> tumor > adjacent.
    reg = ("H2 as registered predicts RSI HIGHER in tumor; observed shift is "
           f"{'HIGHER' if hl > 0 else 'LOWER'} in tumor -> "
           f"{'CONSISTENT' if hl > 0 else 'OPPOSITE'} with the registered "
           "direction.")
    print(f"  {reg}")

    out = {"plan": "PREREGISTRATION 6f", "cohort": a.cohort,
           "n_tumor": len(t), "n_adjacent": len(n),
           "unpaired": {"U": u, "z": round(z_u, 4), "p": p_u, "hodges_lehmann": round(hl, 4),
                        "ci95": [round(ci[0], 4), round(ci[1], 4)],
                        "cliffs_delta": round(delta, 4)},
           "paired": paired, "registered_direction": reg}
    path = os.path.join(a.outdir, "gse76427_h2.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
