#!/usr/bin/env python3
"""32_protein_validation.py -- PREREGISTRATION 6s.

Protein-level corroboration of the module dissociation in Gao et al., Cell
2019;179:561-577.e22 (Table S1, 6,478 gene-level proteins quantified across all
159 paired tumor and adjacent liver samples).

The estimand, the scoring, and the retention fraction are identical in form to
the transcript analysis, and are computed by the same ols_ci routine. The
directional prediction fixed in 6s is that REDUCTION retains at least 50% of its
unadjusted shift and DRAIN retains less than 50%.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

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
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


_da = _load("diffadj", "12_differentiation_adjust.py")
_ca = _load("compadj", "18_composition_adjust.py")
ols_ci, D1, D2_EXTRA = _da.ols_ci, _da.D1, _da.D2_EXTRA
C1 = _ca.C1
import rsi_config as cfg  # noqa: E402


# Symbols retired in the 2019 annotation used by the source table. Only
# unambiguous, one-to-one renames, applied so that a locked gene list written in
# current symbols does not silently lose a measured protein.
ALIAS = {"MTARC1": "MARC1", "MTARC2": "MARC2", "G6PC1": "G6PC"}

MIN_MODULE_MEMBERS = 2   # 6s coverage rule
MIN_D1_PROTEINS = 8      # 6s coverage rule


def resolve(genes, index):
    out = []
    for g in genes:
        if g in index:
            out.append(g)
        elif ALIAS.get(g) in index:
            out.append(ALIAS[g])
    return out


def zmean(mat: pd.DataFrame, genes):
    g = resolve(genes, set(mat.index))
    if not g:
        return None, 0
    sub = mat.loc[g]
    sd = sub.std(axis=1)
    sub = sub[(sd > 0).values]
    if sub.empty:
        return None, 0
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), sub.shape[0]


def main() -> None:
    src = sys.argv[1] if len(sys.argv) > 1 else f"{IR_DOCS}/gao_proteins.tsv"
    label = sys.argv[2] if len(sys.argv) > 2 else "Gao2019"
    df = pd.read_csv(src, sep="\t", index_col=0, low_memory=False)
    df.index = [str(i).strip() for i in df.index]
    df = df[[c for c in df.columns if str(c).strip()]]
    df = df.apply(pd.to_numeric, errors="coerce")

    tum = [c for c in df.columns if str(c).startswith("T")]
    non = [c for c in df.columns if str(c).startswith("N")]
    pt = {c[1:]: c for c in tum}
    pn = {c[1:]: c for c in non}
    pids = sorted(set(pt) & set(pn))
    print(f"=== {label} protein validation (PREREG 6s) ===")
    print(f"  {df.shape[0]} proteins, {len(tum)} tumor, {len(non)} adjacent, "
          f"{len(pids)} pairs")
    miss = int(df.isna().sum().sum())
    print(f"  missing values in the matrix: {miss}")
    lo, hi = float(np.nanmin(df.values)), float(np.nanmax(df.values))
    print(f"  value range {lo:.4g} to {hi:.4g}")
    if hi > 100:
        print("  values are not on a log scale; applying log2(x + 1) before "
              "z scoring, as for linear array data")
        df = np.log2(df.clip(lower=0) + 1.0)

    res = {"plan": "PREREGISTRATION 6s",
           "source": label,
           "n_proteins": int(df.shape[0]), "n_pairs": len(pids)}

    # --- coverage, reported before any contrast (6s) ------------------------
    cov = {}
    for name, genes in [("REDUCTION", cfg.MODULE_REDUCTION),
                        ("DRAIN", cfg.MODULE_DRAIN),
                        ("D1", D1), ("D2", D1 + D2_EXTRA), ("C1", C1)]:
        found = resolve(genes, set(df.index))
        cov[name] = {"n_found": len(found), "n_defined": len(genes),
                     "found": found,
                     "missing": [g for g in genes
                                 if g not in found and ALIAS.get(g) not in found]}
        print(f"  coverage {name:<10} {len(found)}/{len(genes)}")
    res["coverage"] = cov
    for m in ("REDUCTION", "DRAIN"):
        if cov[m]["n_found"] < MIN_MODULE_MEMBERS:
            raise SystemExit(f"{m} has fewer than {MIN_MODULE_MEMBERS} measured "
                             f"members; reported as not testable under 6s")
    if cov["D1"]["n_found"] < MIN_D1_PROTEINS:
        raise SystemExit("fewer than 8 D1 proteins measured; the adjustment is "
                         "reported as not performed under 6s")

    scores = {}
    for name, genes in [("REDUCTION", cfg.MODULE_REDUCTION),
                        ("DRAIN", cfg.MODULE_DRAIN),
                        ("D1", D1), ("D2", D1 + D2_EXTRA), ("C1", C1)]:
        s, n = zmean(df, genes)
        scores[name] = s

    def paired(series):
        return np.array([series[pt[p]] - series[pn[p]] for p in pids], float)

    d_d1 = paired(scores["D1"])
    d_d2 = paired(scores["D2"])
    d_c1 = paired(scores["C1"])

    # positive control on the covariate itself: identity must fall in tumor
    p_d1 = float(wilcoxon(d_d1, alternative="two-sided").pvalue)
    print(f"\n  [1] positive control, D1 paired difference mean {d_d1.mean():+.3f}"
          f"  P = {p_d1:.3g}")
    premise_ok = bool(d_d1.mean() < 0 and p_d1 < 0.05)
    res["positive_control_D1"] = {"mean_paired_delta": round(float(d_d1.mean()), 4),
                                  "p": p_d1,
                                  "direction": "LOWER in tumor" if d_d1.mean() < 0
                                  else "HIGHER in tumor",
                                  "identity_confound_present": premise_ok}
    if not premise_ok:
        print("      !! the hepatocyte identity score is NOT lower in tumor in "
              "this matrix. The confound the retention fraction is built to "
              "remove is absent, so every retention value below is "
              "UNINTERPRETABLE and is reported as not testable under 6s.")

    print("\n  [2] modules, tumor minus paired adjacent liver")
    out = {}
    for mod in ("REDUCTION", "DRAIN"):
        dv = paired(scores[mod])
        unadj = float(dv.mean())
        p_un = float(wilcoxon(dv, alternative="two-sided").pvalue)
        rec = {"n_members": cov[mod]["n_found"],
               "unadjusted_mean_delta": round(unadj, 4),
               "unadjusted_median_delta": round(float(np.median(dv)), 4),
               "unadjusted_p": p_un}
        n = len(pids)
        for tag, X in [("D1", np.column_stack([np.ones(n), d_d1])),
                       ("D2", np.column_stack([np.ones(n), d_d2])),
                       ("joint_D1_C1", np.column_stack([np.ones(n), d_d1, d_c1]))]:
            names = ["intercept"] + (["dD1", "dC1"] if "joint" in tag
                                     else [f"d{tag}"])
            fit = ols_ci(dv, X, names)
            rec[tag] = {"intercept": fit["intercept"],
                        "retention": round(abs(fit["intercept"]["beta"]) /
                                           abs(unadj), 4)}
        out[mod] = rec
        print(f"      {mod:<10} unadjusted {unadj:+.3f} (P = {p_un:.3g})")
        for tag in ("D1", "D2", "joint_D1_C1"):
            b = rec[tag]["intercept"]
            print(f"        {tag:<12} intercept {b['beta']:+.3f} "
                  f"({b['ci95'][0]:+.3f} to {b['ci95'][1]:+.3f}), "
                  f"P = {b['p']:.3g}, retains {rec[tag]['retention']*100:.1f}%")
    res["modules"] = out

    # --- the 6s verdict ------------------------------------------------------
    r_red = out["REDUCTION"]["D1"]["retention"]
    r_dra = out["DRAIN"]["D1"]["retention"]
    dir_ok = (out["REDUCTION"]["unadjusted_mean_delta"] > 0 and
              out["DRAIN"]["unadjusted_mean_delta"] < 0)
    pred_ok = dir_ok and r_red >= 0.5 and r_dra < 0.5
    res["verdict"] = {
        "directions_as_predicted": bool(dir_ok),
        "reduction_retention_D1": r_red,
        "drain_retention_D1": r_dra,
        "prediction_supported": bool(pred_ok),
    }
    print(f"\n  [3] 6s verdict: directions as predicted = {dir_ok}; "
          f"REDUCTION retains {r_red*100:.1f}% (predicted >= 50), "
          f"DRAIN retains {r_dra*100:.1f}% (predicted < 50) "
          f"-> prediction {'SUPPORTED' if pred_ok else 'NOT SUPPORTED'}")

    path = f"results/PROTEIN_VALIDATION_6s_{label}.json"
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"\n  wrote {path}")


if __name__ == "__main__":
    main()
