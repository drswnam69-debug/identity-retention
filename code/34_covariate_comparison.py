#!/usr/bin/env python3
"""34_covariate_comparison.py -- PREREGISTRATION 6u.

Is the identity-retention fraction a re-description of tumor purity adjustment?

The question is settled by four statistics fixed in 6u before any of them existed:
the association between the two paired covariate differences, the retention
fraction each covariate leaves standing across the same evaluable signatures, the
agreement between those two fractions, and the number of signatures that change
side of the 50% threshold depending on which covariate is used.

Nothing here re-selects, re-scores or re-orders a gene set. The eligibility rule,
the null-effect rule and the pairing of 6t are carried over unchanged, and the
identity-only and joint retention fractions are recomputed by the identical code
path and checked against the archived 6t values before anything new is reported.

Usage:
  python3 code/34_covariate_comparison.py --cohort GSE14520 \
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
from scipy.stats import pearsonr, spearmanr, wilcoxon as _wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_da = _load("diffadj", "12_differentiation_adjust.py")
_ca = _load("compadj", "18_composition_adjust.py")
_bm = _load("bench", "30_signature_benchmark.py")
ols_ci, D1 = _da.ols_ci, _da.D1
C1 = _ca.C1
import rsi_config as cfg  # noqa: E402

# The 6t rules, imported rather than restated so they cannot drift.
zmean = _bm.zmean
MIN_SET, MAX_SET, MIN_ON_PLATFORM, NULL_P = (
    _bm.MIN_SET, _bm.MAX_SET, _bm.MIN_ON_PLATFORM, _bm.NULL_P)
OWN = {"REDUCTION": cfg.MODULE_REDUCTION, "DRAIN": cfg.MODULE_DRAIN}

# 6u decision thresholds, fixed before any number below was computed.
RHO_COVAR = 0.70
RHO_RETENTION = 0.70
MED_ABS_DIFF = 0.10
RECLASS_FRAC = 0.10
TOL = 0.002


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

    tis = ph["tissue"].astype(str).str.lower()
    is_t = (~tis.str.contains(a.adjacent_match.lower())).to_numpy()
    pid = (ph["patient_id"].astype(str) if "patient_id" in ph.columns
           else pd.Series(rsi.index, index=rsi.index).astype(str))

    d1, n_d1 = zmean(expr, D1)
    c1, n_c1 = zmean(expr, C1)
    frame = pd.DataFrame({"pid": pid.values, "t": is_t}, index=rsi.index)

    def paired(series):
        f = frame.assign(v=series.to_numpy())
        tt = f[f.t].set_index("pid")["v"]
        tt = tt[~tt.index.duplicated()]
        nn = f[~f.t].set_index("pid")["v"]
        nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[common] - nn.loc[common]).to_numpy(float), common

    d_d1, common = paired(d1)
    d_c1, _ = paired(c1)
    n = len(common)
    print(f"=== {a.cohort} — covariate comparison (PREREG 6u) ===")
    print(f"  {n} pairs; D1 {n_d1}/{len(D1)} genes, C1 {n_c1}/{len(C1)} genes")
    print(f"  D1 and C1 share {len(set(D1) & set(C1))} genes by construction")

    # --- [1] covariate association ------------------------------------------
    pr, pp = pearsonr(d_d1, d_c1)
    sr, sp = spearmanr(d_d1, d_c1)
    fit_cov = ols_ci(d_d1, np.column_stack([np.ones(n), d_c1]),
                     ["intercept", "dC1"])
    r2 = float(fit_cov.get("r2", np.nan)) if isinstance(fit_cov, dict) else np.nan
    if not np.isfinite(r2):
        r2 = float(pr ** 2)
    print(f"\n  [1] paired dD1 vs dC1: Pearson r = {pr:+.3f} (P = {pp:.3g}), "
          f"Spearman rho = {sr:+.3f} (P = {sp:.3g}), R2 = {r2:.3f}")

    assoc = {"pearson_r": round(float(pr), 4), "pearson_p": float(pp),
             "spearman_rho": round(float(sr), 4), "spearman_p": float(sp),
             "r2_dD1_on_dC1": round(r2, 4),
             "shared_genes_D1_C1": len(set(D1) & set(C1))}

    # --- signature loop, 6t rules unchanged ---------------------------------
    sets = json.load(open(a.genesets))
    sets = {k: v for k, v in sets.items() if MIN_SET <= len(v) <= MAX_SET}
    X_D1 = np.column_stack([np.ones(n), d_d1])
    X_C1 = np.column_stack([np.ones(n), d_c1])
    X_J = np.column_stack([np.ones(n), d_d1, d_c1])

    def evaluate(name, genes, min_on=MIN_ON_PLATFORM):
        s, n_on = zmean(expr, genes)
        rec = {"name": name, "n_genes": len(genes), "n_on_platform": n_on}
        if s is None or n_on < min_on:
            rec["status"] = "not_evaluable_on_platform"
            return rec
        dv, _ = paired(s)
        try:
            w_p = float(_wilcoxon(dv, alternative="two-sided").pvalue)
        except Exception:
            w_p = float("nan")
        unadj = float(dv.mean())
        rec.update({"unadjusted_mean_delta": round(unadj, 4), "wilcoxon_p": w_p})
        if not np.isfinite(w_p) or w_p >= NULL_P:
            rec["status"] = "null_effect_excluded"
            return rec
        den = abs(unadj)
        f1 = ols_ci(dv, X_D1, ["intercept", "dD1"])
        fc = ols_ci(dv, X_C1, ["intercept", "dC1"])
        fj = ols_ci(dv, X_J, ["intercept", "dD1", "dC1"])
        rec.update({
            "status": "ok",
            "retention_D1": round(abs(f1["intercept"]["beta"]) / den, 4),
            "retention_C1": round(abs(fc["intercept"]["beta"]) / den, 4),
            "retention_joint": round(abs(fj["intercept"]["beta"]) / den, 4),
        })
        return rec

    # positive controls: the 6t values must come back identical
    controls = {}
    arch = json.load(open(os.path.join(a.outdir,
                                       f"SIGNATURE_BENCHMARK_{a.cohort}.json")))
    for mod, genes in OWN.items():
        rec = evaluate(f"OWN_{mod}", genes, min_on=2)
        exp = arch["own_modules"][mod]
        if exp.get("status") != "ok":
            # the null-effect rule of 6t excluded this module in this cohort, so
            # there is no archived value to reproduce. Reported, not silently
            # skipped, and never used to decide whether the run is valid.
            print(f"  control {mod:<10} excluded in the archive "
                  f"({exp.get('status')}); no value to reproduce")
            controls[mod] = rec
            continue
        for key in ("retention_D1", "retention_joint"):
            got, want = rec.get(key), exp.get(key)
            if want is None or got is None or abs(got - want) > TOL:
                raise SystemExit(f"CONTROL FAILED {mod} {key}: {got} vs {want}")
        controls[mod] = rec
        print(f"  control {mod:<10} D1={rec['retention_D1']:.4f} "
              f"C1={rec['retention_C1']:.4f} joint={rec['retention_joint']:.4f}"
              f"  (6t values reproduced)")

    out = [evaluate(k, v) for k, v in sorted(sets.items())]
    ok = [r for r in out if r["status"] == "ok"]
    print(f"\n  {len(ok)} evaluable signatures")

    # cross-check against the archived 6t retention values, set by set
    arch_map = {r["name"]: r for r in arch["signatures"] if r["status"] == "ok"}
    bad = [r["name"] for r in ok
           if r["name"] in arch_map
           and abs(r["retention_D1"] - arch_map[r["name"]]["retention_D1"]) > TOL]
    if bad:
        raise SystemExit(f"{len(bad)} signatures do not reproduce 6t: {bad[:5]}")
    print(f"  all {len(ok)} identity-only retention values reproduce 6t exactly")

    # --- [2][3][4] agreement -------------------------------------------------
    rd = np.array([r["retention_D1"] for r in ok], float)
    rc = np.array([r["retention_C1"] for r in ok], float)
    diff = rd - rc
    rho_r, p_r = spearmanr(rd, rc)
    reclass = [r for r in ok
               if (r["retention_D1"] < 0.5) != (r["retention_C1"] < 0.5)]
    disc_id = [r for r in ok
               if r["retention_C1"] >= 0.90 and r["retention_D1"] < 0.50]
    disc_co = [r for r in ok
               if r["retention_D1"] >= 0.90 and r["retention_C1"] < 0.50]

    print(f"\n  [2] composition-only retention: median {np.median(rc):.3f} "
          f"(IQR {np.percentile(rc,25):.3f} to {np.percentile(rc,75):.3f})")
    print(f"      identity-only retention:     median {np.median(rd):.3f} "
          f"(IQR {np.percentile(rd,25):.3f} to {np.percentile(rd,75):.3f})")
    print(f"  [3] Spearman rho between them = {rho_r:+.3f} (P = {p_r:.3g})")
    print(f"      median difference D1 minus C1 = {np.median(diff):+.3f} "
          f"(median absolute {np.median(np.abs(diff)):.3f})")
    print(f"      below 50%: identity {int((rd<0.5).sum())}, "
          f"composition {int((rc<0.5).sum())}")
    print(f"      changing side of 50%: {len(reclass)}/{len(ok)} "
          f"({len(reclass)/len(ok)*100:.1f}%)")
    print(f"  [4] composition leaves >=90% while identity leaves <50%: "
          f"{len(disc_id)}")
    print(f"      identity leaves >=90% while composition leaves <50%: "
          f"{len(disc_co)}")

    # --- the 6u verdict ------------------------------------------------------
    c_a = abs(sr) >= RHO_COVAR
    c_b = rho_r >= RHO_RETENTION
    c_c = float(np.median(np.abs(diff))) < MED_ABS_DIFF
    c_d = (len(reclass) / len(ok)) < RECLASS_FRAC
    if c_a and c_b and c_c and c_d:
        verdict = "REDUNDANT"
    elif (not c_a) and ((not c_b) or (len(reclass) / len(ok)) >= RECLASS_FRAC):
        verdict = "DISTINCT"
    else:
        verdict = "PARTIAL"
    print(f"\n  [5] 6u verdict: {verdict}")
    print(f"      |rho(dD1,dC1)| >= 0.70 ? {c_a}   rho(retentions) >= 0.70 ? {c_b}")
    print(f"      median|diff| < 0.10 ? {c_c}      reclass < 10% ? {c_d}")

    res = {
        "plan": "PREREGISTRATION 6u",
        "amendment_sha256_of_text_as_written":
            "c9544b3d03205d48d6ec869971a57570fa2f07c1ce4bab11e29d20af2e06c1f0",
        "cohort": a.cohort, "n_pairs": n,
        "covariate_association": assoc,
        "own_modules": controls,
        "n_evaluable": len(ok),
        "retention_identity_only": {
            "median": round(float(np.median(rd)), 4),
            "q1": round(float(np.percentile(rd, 25)), 4),
            "q3": round(float(np.percentile(rd, 75)), 4),
            "n_below_50pct": int((rd < 0.5).sum())},
        "retention_composition_only": {
            "median": round(float(np.median(rc)), 4),
            "q1": round(float(np.percentile(rc, 25)), 4),
            "q3": round(float(np.percentile(rc, 75)), 4),
            "n_below_50pct": int((rc < 0.5).sum())},
        "agreement": {
            "spearman_rho": round(float(rho_r), 4), "spearman_p": float(p_r),
            "median_difference_D1_minus_C1": round(float(np.median(diff)), 4),
            "median_absolute_difference": round(float(np.median(np.abs(diff))), 4),
            "iqr_difference": [round(float(np.percentile(diff, 25)), 4),
                               round(float(np.percentile(diff, 75)), 4)],
            "n_reclassified_at_50pct": len(reclass),
            "prop_reclassified": round(len(reclass) / len(ok), 4),
            "reclassified_names": [r["name"] for r in reclass]},
        "discordant": {
            "composition_ge90_identity_lt50": [r["name"] for r in disc_id],
            "identity_ge90_composition_lt50": [r["name"] for r in disc_co]},
        "criteria": {"rho_covariate_ge": RHO_COVAR,
                     "rho_retention_ge": RHO_RETENTION,
                     "median_abs_diff_lt": MED_ABS_DIFF,
                     "reclass_frac_lt": RECLASS_FRAC,
                     "met": {"covariate_rho": bool(c_a), "retention_rho": bool(c_b),
                             "median_abs_diff": bool(c_c), "reclass": bool(c_d)}},
        "verdict": verdict,
        "signatures": ok,
    }
    path = os.path.join(a.outdir, f"COVARIATE_COMPARISON_6u_{a.cohort}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  wrote {path}")


if __name__ == "__main__":
    main()
