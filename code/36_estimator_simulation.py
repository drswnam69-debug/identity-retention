#!/usr/bin/env python3
"""36_estimator_simulation.py -- PREREGISTRATION 6v, Part B.

Does the identity-retention estimator recover a known retention fraction?

Synthetic paired differences are generated on the REAL dD1 and dC1 vectors from
the 213 GSE14520 pairs, so the covariate geometry the estimator faces, including
the scarcity of dD1 near zero that makes the intercept an extrapolation, is the
real one. The true retention is tau by construction.

The grid, the replicate counts, the seed and the outcomes accepted are fixed in
6v and are not modified here.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


_da = _load("diffadj", "12_differentiation_adjust.py")
_ca = _load("compadj", "18_composition_adjust.py")
_bm = _load("bench", "30_signature_benchmark.py")
D1, C1, zmean = _da.D1, _ca.C1, _bm.zmean

# --- 6v grid, fixed before any simulated dataset existed --------------------
TAUS = [0.10, 0.25, 0.50, 0.75, 1.00]
SIGMAS = [0.25, 0.50, 1.00]
COMP_SHARES = [0.0, 0.3]
MU_MAIN = 0.50
EXTRA_MUS = [0.25, 1.00]
EXTRA_TAUS = [0.25, 0.75]
N_REP = 500
N_BOOT = 500
SEED = 20260907


def batched_intercept(Y, X, idx):
    """Y (R, n) responses, X (n, p) design, idx (B, n). Returns (R, B)."""
    Xb = X[idx]                                        # (B, n, p)
    XtX = np.einsum("bij,bik->bjk", Xb, Xb)            # (B, p, p)
    Xtyi = np.einsum("bij,rbi->rbj", Xb, Y[:, idx])    # (R, B, p)
    beta = np.linalg.solve(XtX[None], Xtyi[..., None])[..., 0]
    return beta[..., 0]


def main() -> None:
    rsi = pd.read_csv("results/GSE14520/rsi.tsv", sep="\t", index_col=0)
    ph = pd.read_csv("results/GSE14520/phenotype.tsv", sep="\t",
                     index_col=0).loc[rsi.index]
    expr = pd.read_csv("data/GSE14520_symbols.tsv.gz", sep="\t",
                       index_col=0)[rsi.index]
    is_t = (~ph["tissue"].astype(str).str.lower()
            .str.contains("adjacent")).to_numpy()
    fr = pd.DataFrame({"pid": ph["patient_id"].astype(str).values, "t": is_t},
                      index=rsi.index)

    def paired(s):
        f = fr.assign(v=s.to_numpy())
        tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
        nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
        k = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[k] - nn.loc[k]).to_numpy(float)

    dD1 = paired(zmean(expr, D1)[0])
    dC1 = paired(zmean(expr, C1)[0])
    n = len(dD1)
    mD1, mC1 = float(dD1.mean()), float(dC1.mean())
    X = np.column_stack([np.ones(n), dD1, dC1])
    print("=== estimator simulation (PREREG 6v, part B) ===")
    print(f"  {n} real pairs; mean dD1 = {mD1:+.4f}, mean dC1 = {mC1:+.4f}")
    print(f"  {N_REP} replicates x {N_BOOT} bootstrap per cell, seed {SEED}")
    print(f"  dD1 >= 0 in {int((dD1 >= 0).sum())} of {n} pairs "
          f"(the intercept is an extrapolation)")

    cells = [(t, s, c, MU_MAIN) for t in TAUS for s in SIGMAS for c in COMP_SHARES]
    cells += [(t, 0.50, 0.0, m) for m in EXTRA_MUS for t in EXTRA_TAUS]
    print(f"  {len(cells)} cells\n")

    rng = np.random.default_rng(SEED)
    IDX = rng.integers(0, n, size=(N_BOOT, n))
    rows = []
    for k, (tau, sigma, cshare, mu) in enumerate(cells, 1):
        alpha = tau * mu
        beta = (1 - cshare) * (1 - tau) * mu / mD1
        gamma = cshare * (1 - tau) * mu / mC1
        mean_signal = alpha + beta * mD1 + gamma * mC1
        assert abs(mean_signal - mu) < 1e-9, mean_signal
        base = alpha + beta * dD1 + gamma * dC1              # (n,)
        eps = rng.normal(0.0, sigma, size=(N_REP, n))
        Y = base[None, :] + eps                              # (R, n)

        den = Y.mean(axis=1)                                 # (R,)
        XtX = X.T @ X
        num = np.abs(np.linalg.solve(XtX, X.T @ Y.T)[0])     # (R,)
        est = num / np.abs(den)

        den_b = Y[:, IDX].mean(axis=2)                       # (R, B)
        num_b = np.abs(batched_intercept(Y, X, IDX))         # (R, B)
        d_lo = np.percentile(den_b, 2.5, axis=1)
        d_hi = np.percentile(den_b, 97.5, axis=1)
        stable = ~((d_lo <= 0) & (0 <= d_hi))
        ratio_b = num_b / np.abs(den_b)
        lo = np.percentile(ratio_b, 2.5, axis=1)
        hi = np.percentile(ratio_b, 97.5, axis=1)
        cov = float(((lo <= tau) & (tau <= hi))[stable].mean())
        med = float(np.median(est[stable]))
        rows.append({"tau": tau, "sigma": sigma, "comp_share": cshare, "mu": mu,
                     "median_estimate": round(med, 4),
                     "bias": round(med - tau, 4),
                     "iqr_estimate": [round(float(np.percentile(est[stable], 25)), 4),
                                      round(float(np.percentile(est[stable], 75)), 4)],
                     "coverage_95": round(cov, 4),
                     "n_stable": int(stable.sum()), "n_rep": N_REP})
        print(f"  [{k:2d}/{len(cells)}] tau={tau:.2f} sigma={sigma:.2f} "
              f"comp={cshare:.1f} mu={mu:.2f} -> median {med:.3f} "
              f"bias {med-tau:+.3f}  coverage {cov:.3f}  stable {int(stable.sum())}")

    bias = np.array([r["bias"] for r in rows])
    covg = np.array([r["coverage_95"] for r in rows])
    verdict = ("well_behaved" if np.abs(bias).max() <= 0.05
               and covg.min() >= 0.92 and covg.max() <= 0.97 else "qualified")
    fails_badly = bool(any(r for r in rows
                           if r["sigma"] <= 0.50
                           and (abs(r["bias"]) > 0.20 or r["coverage_95"] < 0.80)))
    print(f"\n  max |bias| = {np.abs(bias).max():.4f} "
          f"(at {rows[int(np.abs(bias).argmax())]['tau']:.2f}/"
          f"{rows[int(np.abs(bias).argmax())]['sigma']:.2f})")
    print(f"  coverage range = {covg.min():.3f} to {covg.max():.3f}")
    print(f"  6v verdict: {verdict}; fails_badly = {fails_badly}")

    res = {"plan": "PREREGISTRATION 6v part B",
           "amendment_sha256_of_text_as_written":
               "4b434bd2dc2af88fe04eb90cb04ec2a6e6dff39ed4de772b8261e5be5aead827",
           "n_pairs": n, "mean_dD1": round(mD1, 4), "mean_dC1": round(mC1, 4),
           "n_pairs_dD1_ge_0": int((dD1 >= 0).sum()),
           "n_rep": N_REP, "n_boot": N_BOOT, "seed": SEED,
           "max_abs_bias": round(float(np.abs(bias).max()), 4),
           "coverage_min": round(float(covg.min()), 4),
           "coverage_max": round(float(covg.max()), 4),
           "verdict": verdict, "fails_badly": fails_badly,
           "cells": rows}
    with open("results/ESTIMATOR_SIMULATION_6v.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print("\n  wrote results/ESTIMATOR_SIMULATION_6v.json")


if __name__ == "__main__":
    main()
