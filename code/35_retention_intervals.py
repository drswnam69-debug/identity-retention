#!/usr/bin/env python3
"""35_retention_intervals.py -- PREREGISTRATION 6v, Part A.

Bootstrap intervals for the identity-retention fraction of every evaluable
signature. The resampling is over the 213 patient pairs, with the module and
covariate scores computed once on the full cohort and held fixed, so the
intervals express sampling variability of patients and not of the within-cohort
z scoring. That limitation is stated in 6v and repeated in the manuscript.

Numerator and denominator are recomputed on the SAME resample, which is the whole
point: the ratio's uncertainty is not the intercept's uncertainty.

Usage:
  python3 code/35_retention_intervals.py --cohort GSE14520 \
      --expr-source data/GSE14520_symbols.tsv.gz
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon as _wilcoxon

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
ols_ci, D1 = _da.ols_ci, _da.D1
C1 = _ca.C1
zmean = _bm.zmean
MIN_SET, MAX_SET, MIN_ON_PLATFORM, NULL_P = (
    _bm.MIN_SET, _bm.MAX_SET, _bm.MIN_ON_PLATFORM, _bm.NULL_P)
import rsi_config as cfg  # noqa: E402

OWN = {"REDUCTION": cfg.MODULE_REDUCTION, "DRAIN": cfg.MODULE_DRAIN}

# 6v settings, fixed before any interval existed.
B = 2000
SEED = 20260907
TOL = 0.002


def batched_intercept(y: np.ndarray, X: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Intercept of OLS of y on X, for each bootstrap index row of idx."""
    Xb = X[idx]                                   # (B, n, p)
    yb = y[idx]                                   # (B, n)
    XtX = np.einsum("bij,bik->bjk", Xb, Xb)
    Xty = np.einsum("bij,bi->bj", Xb, yb)
    beta = np.linalg.solve(XtX, Xty[..., None])[..., 0]   # (B, p)
    return beta[:, 0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--expr-source", required=True)
    ap.add_argument("--genesets", default="genesets/eligible.json")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--adjacent-match", default="adjacent")
    ap.add_argument("--identity", choices=["D1", "K1", "L1"], default="D1",
                    help="identity covariate. D1 is the locked liver score; "
                         "K1 is the kidney score fixed in 6x. Default D1, so "
                         "every previously reported run is unchanged.")
    a = ap.parse_args()

    global D1
    if a.identity == "K1":
        from k1_config import K1 as _K1
        D1 = _K1
        print(f"  identity covariate: K1 ({len(D1)} genes, PREREG 6x)")
    elif a.identity == "L1":
        from l1_config import L1 as _L1
        D1 = _L1
        print(f"  identity covariate: L1 ({len(D1)} genes, PREREG 6y)")
    else:
        print(f"  identity covariate: D1 ({len(D1)} genes, locked)")
    base = os.path.join(a.outdir, a.cohort)
    rsi = pd.read_csv(os.path.join(base, "rsi.tsv"), sep="\t", index_col=0)
    ph = pd.read_csv(os.path.join(base, "phenotype.tsv"), sep="\t",
                     index_col=0).loc[rsi.index]
    expr = pd.read_csv(a.expr_source, sep="\t", index_col=0)[rsi.index]
    is_t = (~ph["tissue"].astype(str).str.lower()
            .str.contains(a.adjacent_match.lower())).to_numpy()
    pid = (ph["patient_id"].astype(str) if "patient_id" in ph.columns
           else pd.Series(rsi.index, index=rsi.index).astype(str))
    fr = pd.DataFrame({"pid": pid.values, "t": is_t}, index=rsi.index)

    def paired(s):
        f = fr.assign(v=s.to_numpy())
        tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
        nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
        k = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[k] - nn.loc[k]).to_numpy(float)

    d_d1 = paired(zmean(expr, D1)[0])
    d_c1 = paired(zmean(expr, C1)[0])
    n = len(d_d1)
    print(f"=== {a.cohort} — retention intervals (PREREG 6v, part A) ===")
    print(f"  {n} pairs, B = {B}, seed {SEED}")

    X_D1 = np.column_stack([np.ones(n), d_d1])
    X_J = np.column_stack([np.ones(n), d_d1, d_c1])
    rng = np.random.default_rng(SEED)
    IDX = rng.integers(0, n, size=(B, n))          # one index set for every signature

    def evaluate(name, genes, min_on=MIN_ON_PLATFORM):
        s, n_on = zmean(expr, genes)
        rec = {"name": name, "n_genes": len(genes), "n_on_platform": n_on}
        if s is None or n_on < min_on:
            rec["status"] = "not_evaluable_on_platform"
            return rec
        dv = paired(s)
        try:
            w_p = float(_wilcoxon(dv, alternative="two-sided").pvalue)
        except Exception:
            w_p = float("nan")
        unadj = float(dv.mean())
        rec.update({"unadjusted_mean_delta": round(unadj, 4), "wilcoxon_p": w_p})
        if not np.isfinite(w_p) or w_p >= NULL_P:
            rec["status"] = "null_effect_excluded"
            return rec
        f1 = ols_ci(dv, X_D1, ["intercept", "dD1"])
        fj = ols_ci(dv, X_J, ["intercept", "dD1", "dC1"])
        rec.update({"status": "ok",
                    "retention_D1": round(abs(f1["intercept"]["beta"]) / abs(unadj), 4),
                    "retention_joint": round(abs(fj["intercept"]["beta"]) / abs(unadj), 4)})

        # --- the bootstrap ---------------------------------------------------
        den = dv[IDX].mean(axis=1)                       # (B,)
        d_lo, d_hi = np.percentile(den, [2.5, 97.5])
        unstable = bool(d_lo <= 0.0 <= d_hi)
        rec["unadjusted_boot_ci95"] = [round(float(d_lo), 4), round(float(d_hi), 4)]
        rec["ratio_unstable"] = unstable
        for tag, X in (("joint", X_J), ("D1", X_D1)):
            num = np.abs(batched_intercept(dv, X, IDX))
            if unstable:
                rec[f"ci95_retention_{tag}"] = None
                continue
            ratio = num / np.abs(den)
            lo, hi = np.percentile(ratio, [2.5, 97.5])
            rec[f"ci95_retention_{tag}"] = [round(float(lo), 4), round(float(hi), 4)]
            rec[f"boot_median_retention_{tag}"] = round(float(np.median(ratio)), 4)
        return rec

    # positive controls first, as in 6u
    arch = json.load(open(os.path.join(a.outdir,
                                       f"SIGNATURE_BENCHMARK_{a.cohort}.json")))
    controls = {}
    for mod, genes in OWN.items():
        rec = evaluate(f"OWN_{mod}", genes, min_on=2)
        exp = arch["own_modules"][mod]
        if exp.get("status") != "ok":
            print(f"  control {mod:<10} excluded in the archive "
                  f"({exp.get('status')}); no value to reproduce")
            controls[mod] = rec
            continue
        for key in ("retention_D1", "retention_joint"):
            if rec.get(key) is None or abs(rec[key] - exp[key]) > TOL:
                raise SystemExit(f"CONTROL FAILED {mod} {key}")
        controls[mod] = rec
        ci = rec["ci95_retention_joint"]
        print(f"  control {mod:<10} {rec['retention_joint']:.4f} "
              f"({ci[0]:.3f} to {ci[1]:.3f})   6t value reproduced")

    sets = json.load(open(a.genesets))
    sets = {k: v for k, v in sets.items() if MIN_SET <= len(v) <= MAX_SET}
    out = [evaluate(k, v) for k, v in sorted(sets.items())]
    ok = [r for r in out if r["status"] == "ok"]
    arch_map = {r["name"]: r for r in arch["signatures"] if r["status"] == "ok"}
    bad = [r["name"] for r in ok if r["name"] in arch_map
           and abs(r["retention_joint"] - arch_map[r["name"]]["retention_joint"]) > TOL]
    if bad:
        raise SystemExit(f"{len(bad)} signatures do not reproduce 6t: {bad[:5]}")
    print(f"  all {len(ok)} point estimates reproduce 6t exactly")

    # --- the quantities fixed in 6v -----------------------------------------
    unstable = [r for r in ok if r["ratio_unstable"]]
    stable = [r for r in ok if not r["ratio_unstable"]]
    below_pt = [r for r in stable if r["retention_joint"] < 0.5]
    below_ci = [r for r in below_pt if r["ci95_retention_joint"][1] < 0.5]
    widths = np.array([r["ci95_retention_joint"][1] - r["ci95_retention_joint"][0]
                       for r in stable], float)
    contains1 = [r for r in stable
                 if r["ci95_retention_joint"][0] <= 1.0 <= r["ci95_retention_joint"][1]]

    print(f"\n  unstable under the 6v guard: {len(unstable)}")
    print(f"  below 0.50 by point estimate: {len(below_pt)}")
    print(f"  of those, whole interval below 0.50: {len(below_ci)}")
    print(f"  median interval width: {np.median(widths):.3f} "
          f"(IQR {np.percentile(widths,25):.3f} to {np.percentile(widths,75):.3f})")
    print(f"  intervals containing 1.0: {len(contains1)} of {len(stable)}")

    res = {"plan": "PREREGISTRATION 6v part A",
           "amendment_sha256_of_text_as_written":
               "4b434bd2dc2af88fe04eb90cb04ec2a6e6dff39ed4de772b8261e5be5aead827",
           "cohort": a.cohort, "n_pairs": n, "B": B, "seed": SEED,
           "n_evaluable": len(ok),
           "n_unstable": len(unstable),
           "unstable_names": [r["name"] for r in unstable],
           "n_below_50_point": len(below_pt),
           "n_below_50_interval": len(below_ci),
           "below_50_interval_names": [r["name"] for r in below_ci],
           "median_ci_width": round(float(np.median(widths)), 4),
           "iqr_ci_width": [round(float(np.percentile(widths, 25)), 4),
                            round(float(np.percentile(widths, 75)), 4)],
           "n_ci_contains_1": len(contains1),
           "own_modules": controls,
           "signatures": ok}
    path = os.path.join(a.outdir, f"RETENTION_INTERVALS_6v_{a.cohort}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  wrote {path}")


if __name__ == "__main__":
    main()
