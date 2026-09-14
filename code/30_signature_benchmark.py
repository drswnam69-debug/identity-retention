#!/usr/bin/env python3
"""30_signature_benchmark.py -- PREREGISTRATION 6t.

Identity retention across published liver expression signatures.

Every published signature is scored, paired, adjusted and summarized by the SAME
functions this study used for its own modules: zmean and the pairing helper from
12_differentiation_adjust.py, the covariates from 18_composition_adjust.py, and
ols_ci for every regression. The study's own REDUCTION and DRAIN modules are run
through the identical loop as positive controls and must reproduce the archived
values, or the run aborts.

Nothing in this script selects, filters or reorders gene sets on the basis of a
result. The eligibility rule is applied before any retention value exists, as
fixed in 6t.

Usage:
  python3 code/30_signature_benchmark.py \
      --cohort GSE14520 \
      --expr-source data/GSE14520_symbols.tsv.gz \
      --genesets genesets/eligible.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stats_lite as sl  # noqa: E402
from scipy.stats import wilcoxon as _wilcoxon  # 6t null-effect rule only


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_da = _load("diffadj", "12_differentiation_adjust.py")
_ca = _load("compadj", "18_composition_adjust.py")
ols_ci, D1 = _da.ols_ci, _da.D1
C1 = _ca.C1
import rsi_config as cfg  # noqa: E402

OWN = {"REDUCTION": cfg.MODULE_REDUCTION, "DRAIN": cfg.MODULE_DRAIN}

# Archived values these positive controls must reproduce (joint D1 + C1 model).
EXPECTED = {
    "GSE14520": {"REDUCTION": 0.8905, "DRAIN": 0.3103},
    "GSE76427": {"REDUCTION": 0.5120, "DRAIN": 0.5112},
}
TOL = 0.002

MIN_SET = 15
MAX_SET = 1000
MIN_ON_PLATFORM = 10
NULL_P = 0.05


def zmean(mat: pd.DataFrame, genes) -> tuple[pd.Series | None, int]:
    """Mean of within-cohort gene-wise z scores. Identical to 6g/6j scoring."""
    g = [x for x in genes if x in mat.index]
    if not g:
        return None, 0
    sub = mat.loc[g]
    sd = sub.std(axis=1)
    keep = sd > 0
    sub, g = sub[keep.values], list(np.asarray(g)[keep.values])
    if not g:
        return None, 0
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), len(g)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--expr-source", required=True,
                    help="symbol x sample matrix, already normalized by "
                         "01_prepare.py; must cover the signature genes")
    ap.add_argument("--genesets", default="genesets/eligible.json")
    ap.add_argument("--adjacent-match", default="adjacent")
    ap.add_argument("--identity", choices=["D1", "K1", "L1"], default="D1",
                    help="identity covariate. D1 is the locked liver score; "
                         "K1 is the kidney score fixed in 6x. Default D1, so "
                         "every previously reported run is unchanged.")
    ap.add_argument("--controls", choices=["enforce", "advisory"],
                    default="enforce",
                    help="advisory: report the positive controls but do not "
                         "abort. Only for a cohort whose archived expression "
                         "source cannot be rebuilt, which must be disclosed.")
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
    expr = pd.read_csv(a.expr_source, sep="\t", index_col=0)
    missing = [s for s in rsi.index if s not in expr.columns]
    if missing:
        raise SystemExit(f"expression source lacks {len(missing)} samples, "
                         f"first few {missing[:4]}")
    expr = expr[rsi.index]

    tis = ph["tissue"].astype(str).str.lower()
    is_t = (~tis.str.contains(a.adjacent_match.lower())).to_numpy()
    pid = (ph["patient_id"].astype(str) if "patient_id" in ph.columns
           else pd.Series(rsi.index, index=rsi.index).astype(str))

    d1, n_d1 = zmean(expr, D1)
    c1, n_c1 = zmean(expr, C1)
    if d1 is None or c1 is None:
        raise SystemExit("covariates could not be built from the expression source")
    print(f"=== {a.cohort} — signature identity retention (PREREG 6t) ===")
    print(f"  {int(is_t.sum())} tumor / {int((~is_t).sum())} adjacent")
    print(f"  D1 {n_d1}/{len(D1)} genes, C1 {n_c1}/{len(C1)} genes")

    frame = pd.DataFrame({"pid": pid.values, "t": is_t,
                          "d1": d1.to_numpy(), "c1": c1.to_numpy()},
                         index=rsi.index)

    def paired(series: pd.Series):
        f = frame.assign(v=series.to_numpy())
        tt = f[f.t].set_index("pid")["v"]
        tt = tt[~tt.index.duplicated()]
        nn = f[~f.t].set_index("pid")["v"]
        nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[common] - nn.loc[common]).to_numpy(float), common

    d_d1, common = paired(d1)
    d_c1, _ = paired(c1)
    n_pairs = len(common)
    print(f"  {n_pairs} pairs")
    if n_pairs < 10:
        raise SystemExit("too few pairs to run the benchmark")

    sets = json.load(open(a.genesets))
    sets = {k: v for k, v in sets.items() if MIN_SET <= len(v) <= MAX_SET}
    covar = set(D1) | set(C1)
    panel = set(cfg.MODULE_SUPPLY) | set(cfg.MODULE_REDUCTION) | set(cfg.MODULE_DRAIN)

    def evaluate(name, genes, kind, min_on=MIN_ON_PLATFORM):
        s, n_on = zmean(expr, genes)
        rec = {"name": name, "kind": kind, "n_genes": len(genes),
               "n_on_platform": n_on,
               "overlap_D1": len(set(genes) & set(D1)),
               "overlap_C1": len(set(genes) & set(C1)),
               "overlap_panel": len(set(genes) & panel)}
        if s is None or n_on < min_on:
            rec["status"] = "not_evaluable_on_platform"
            return rec
        dv, _ = paired(s)
        try:
            w_p = float(_wilcoxon(dv, alternative="two-sided").pvalue)
        except Exception:
            w_p = float("nan")
        unadj = float(dv.mean())
        rec.update({"unadjusted_mean_delta": round(unadj, 4),
                    "unadjusted_median_delta": round(float(np.median(dv)), 4),
                    "wilcoxon_p": w_p})
        if not np.isfinite(w_p) or w_p >= NULL_P:
            rec["status"] = "null_effect_excluded"
            return rec
        X1 = np.column_stack([np.ones(n_pairs), d_d1])
        XJ = np.column_stack([np.ones(n_pairs), d_d1, d_c1])
        f1 = ols_ci(dv, X1, ["intercept", "dD1"])
        fj = ols_ci(dv, XJ, ["intercept", "dD1", "dC1"])
        denom = abs(unadj)
        rec.update({
            "status": "ok",
            "D1_intercept": f1["intercept"],
            "joint_intercept": fj["intercept"],
            "retention_D1": round(abs(f1["intercept"]["beta"]) / denom, 4),
            "retention_joint": round(abs(fj["intercept"]["beta"]) / denom, 4),
        })
        return rec

    # --- positive controls, identical code path -----------------------------
    controls = {}
    for mod, genes in OWN.items():
        # the study's own modules are 3 and 4 genes; the 10-gene platform rule
        # in 6t governs published signatures, not these controls
        rec = evaluate(f"OWN_{mod}", genes, "own_module", min_on=2)
        controls[mod] = rec
        exp = EXPECTED.get(a.cohort, {}).get(mod)
        got = rec.get("retention_joint")
        print(f"  control {mod:<10} retention_joint={got}  expected={exp}")
        if exp is not None and abs((got if got is not None else -9) - exp) > TOL:
            if a.controls == "advisory":
                print(f"  !! control {mod} does not reproduce the archived "
                      f"value; reported as advisory, see the run notes")
            else:
                raise SystemExit(
                    f"POSITIVE CONTROL FAILED for {mod}: got {got}, "
                    f"archived {exp}. The benchmark is not comparable to the "
                    f"study's own numbers; fix the input before proceeding.")
    res_controls_mode = a.controls
    print(f"  positive controls: {a.controls}")

    # --- the published signatures -------------------------------------------
    out = []
    for i, (name, genes) in enumerate(sorted(sets.items()), 1):
        rec = evaluate(name, genes, "published_signature")
        trimmed = [g for g in genes if g not in covar]
        if len(trimmed) != len(genes):
            rec["sensitivity_no_covariate_genes"] = evaluate(
                name + "__trimmed", trimmed, "sensitivity")
        out.append(rec)
        if i % 20 == 0:
            print(f"    {i}/{len(sets)} sets scored")

    ok = [r for r in out if r["status"] == "ok"]
    res = {
        "plan": "PREREGISTRATION 6t",
        "cohort": a.cohort,
        "n_pairs": n_pairs,
        "genes": {"D1": n_d1, "C1": n_c1},
        "eligibility": {"min_set": MIN_SET, "max_set": MAX_SET,
                        "min_on_platform": MIN_ON_PLATFORM,
                        "null_effect_p": NULL_P},
        "n_sets_supplied": len(sets),
        "n_sets_evaluable": len(ok),
        "n_sets_not_on_platform": sum(
            1 for r in out if r["status"] == "not_evaluable_on_platform"),
        "n_sets_null_effect": sum(
            1 for r in out if r["status"] == "null_effect_excluded"),
        "controls_mode": res_controls_mode,
        "own_modules": controls,
        "signatures": out,
    }
    if ok:
        vals = np.array([r["retention_joint"] for r in ok], float)
        res["summary_joint"] = {
            "median": round(float(np.median(vals)), 4),
            "q1": round(float(np.percentile(vals, 25)), 4),
            "q3": round(float(np.percentile(vals, 75)), 4),
            "n_below_50pct": int((vals < 0.5).sum()),
            "prop_below_50pct": round(float((vals < 0.5).mean()), 4),
        }
        for mod in OWN:
            v = controls[mod].get("retention_joint")
            if v is not None:
                res["summary_joint"][f"percentile_of_{mod}"] = round(
                    float((vals < v).mean() * 100), 1)

    path = os.path.join(a.outdir, f"SIGNATURE_BENCHMARK_{a.cohort}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  wrote {path}")
    if "summary_joint" in res:
        s = res["summary_joint"]
        print(f"  median retention {s['median']:.3f} "
              f"(IQR {s['q1']:.3f} to {s['q3']:.3f}); "
              f"{s['n_below_50pct']}/{len(ok)} below 50%")


if __name__ == "__main__":
    main()
