#!/usr/bin/env python3
"""12_differentiation_adjust.py -- PREREGISTRATION 6g.

Is the tumor-vs-adjacent RSI difference explained by loss of hepatocyte
identity? Index values are read as written by 02_compute_rsi.py; the
differentiation panel is read as written under 6e/6g. Nothing is renormalized.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_lite as sl
import rsi_config as cfg

D1 = ["HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4",
      "ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2",
      "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1"]
D2_EXTRA = ["CYP2E1", "CYP3A4", "CYP1A2", "CYP2C9"]


def zmean(mat: pd.DataFrame, genes) -> pd.Series:
    g = [x for x in genes if x in mat.index]
    missing = sorted(set(genes) - set(g))
    if missing:
        print(f"    !! absent from the platform: {missing}")
    z = mat.loc[g].sub(mat.loc[g].mean(axis=1), axis=0) \
                  .div(mat.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0), len(g)


def ols_ci(y, X, names):
    """OLS with 95% CIs and two-sided t tests."""
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    s2 = float(resid @ resid) / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    tcrit = _tcrit(dof)   # 6n: the 1.96 shortcut above dof=200 was a z, not a t
    out = {}
    for i, nm in enumerate(names):
        out[nm] = {"beta": round(float(beta[i]), 4),
                   "se": round(float(se[i]), 4),
                   "ci95": [round(float(beta[i] - tcrit * se[i]), 4),
                            round(float(beta[i] + tcrit * se[i]), 4)],
                   "t": round(float(t[i]), 3),
                   "p": float(sl._t_two_sided(float(t[i]), dof))}
    return out


def _tcrit(dof, lo=1.5, hi=6.0):
    for _ in range(200):
        mid = (lo + hi) / 2
        if sl._t_two_sided(mid, dof) > 0.05:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="GSE76427")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--diff-panel",
                    default="data/GSE76427_diff_panel_qn.tsv")
    ap.add_argument("--diff-map",
                    default="data/GPL10558_differentiation_probe_to_symbol.tsv")
    ap.add_argument("--diff-from-expr", action="store_true",
                    help="differentiation genes are already symbols in the "
                         "cohort's expr_log (no separate panel file)")
    ap.add_argument("--adjacent-match", default="adjacent",
                    help="substring identifying adjacent tissue in `tissue`")
    ap.add_argument("--tumor-match", default=None,
                    help="substring identifying tumor; if given, samples "
                         "matching neither label are dropped")
    ap.add_argument("--unpaired-only", action="store_true",
                    help="no declared patient pairing in this series")
    ap.add_argument("--tag", default=None, help="output file stem")
    a = ap.parse_args()

    base = os.path.join(a.outdir, a.cohort)
    rsi = pd.read_csv(os.path.join(base, "rsi.tsv"), sep="\t", index_col=0)
    ph = pd.read_csv(os.path.join(base, "phenotype.tsv"), sep="\t",
                     index_col=0).loc[rsi.index]
    expr = pd.read_csv(os.path.join(base, "expr_log.tsv.gz"), sep="\t",
                       index_col=0)[rsi.index]

    tis = ph["tissue"].astype(str).str.lower()
    if a.tumor_match:
        keep = tis.str.contains(a.adjacent_match.lower()) | \
               tis.str.contains(a.tumor_match.lower())
        n_drop = int((~keep).sum())
        if n_drop:
            print(f"  restricted to the two tissue groups; dropped {n_drop} "
                  "samples from other groups")
        rsi, ph, expr, tis = (rsi[keep.values], ph[keep.values],
                              expr.loc[:, keep.values], tis[keep.values])
    is_t_all = ~tis.str.contains(a.adjacent_match.lower())

    if a.diff_from_expr:
        # differentiation genes already carried as symbols in expr_log
        dp = expr
    else:
        dp = pd.read_csv(a.diff_panel, sep="\t", index_col=0)
        amap = pd.read_csv(a.diff_map, sep="\t", header=None,
                           names=["probe", "symbol"])
        m = dict(zip(amap["probe"], amap["symbol"]))
        dp.index = [m.get(i, "") for i in dp.index]
        dp = dp.loc[dp.index != ""].groupby(level=0).max()
        dp.columns = [c.strip() for c in dp.columns]
        key = [("PT" if t else "ANTT") + p.replace("HCC", "").strip()
               for t, p in zip(is_t_all.values, ph["patient_id"].astype(str))]
        assert set(key) == set(dp.columns), "differentiation panel columns mismatch"
        dp = dp[key]
        dp.columns = ph.index
        assert list(dp.columns) == list(rsi.index)

    print(f"=== {a.cohort} — differentiation-adjusted H2 (PREREG 6g) ===")
    d1, n1 = zmean(dp, D1)
    d2, n2 = zmean(dp, D1 + D2_EXTRA)
    print(f"  D1 hepatocyte identity: {n1}/{len(D1)} genes")
    print(f"  D2 (+ hepatic P450s):   {n2}/{len(D1 + D2_EXTRA)} genes")

    is_t = is_t_all.to_numpy()
    r = rsi["RSI"].to_numpy(float)

    # 1. positive control
    _, _, p_d1 = sl.mannwhitney_u(d1.to_numpy()[is_t], d1.to_numpy()[~is_t])
    print(f"\n  [1] positive control — D1 in tumor {d1[is_t].mean():+.3f} "
          f"vs adjacent {d1[~is_t].mean():+.3f}   P = {p_d1:.3g}")
    if not (d1[is_t].mean() < d1[~is_t].mean() and p_d1 < 0.05):
        print("      !! D1 is NOT lower in tumor — the confound is absent here;"
              " steps 2-5 are uninterpretable as specified.")

    # 2. primary: paired-difference regression (only where pairing is declared)
    pid = ph["patient_id"].astype(str) if "patient_id" in ph.columns \
        else pd.Series(ph.index, index=ph.index)
    df = pd.DataFrame({"pid": pid.values, "t": is_t, "rsi": r,
                       "d1": d1.to_numpy(), "d2": d2.to_numpy()},
                      index=rsi.index)
    for mod, genes in (("reduction", cfg.MODULE_REDUCTION),
                       ("drain", cfg.MODULE_DRAIN)):
        s, _ = zmean(expr, genes)
        df[mod] = s.to_numpy()

    def pairs(col):
        tt = df[df.t].set_index("pid")[col]; tt = tt[~tt.index.duplicated()]
        nn = df[~df.t].set_index("pid")[col]; nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[common] - nn.loc[common]).to_numpy(float), common

    if a.unpaired_only:
        print("\n  [2] PRIMARY paired test SKIPPED — this series declares no "
              "patient ID (PREREG 6h); the unpaired model below is primary here.")
        d_rsi = d_d1 = d_d2 = None
        common = []
    else:
        d_rsi, common = pairs("rsi")
        d_d1, _ = pairs("d1")
        d_d2, _ = pairs("d2")
        print(f"\n  [2] PRIMARY — paired-difference regression, "
              f"{len(common)} pairs")
        print(f"      unadjusted mean ΔRSI = {d_rsi.mean():+.3f}  "
              f"(median {np.median(d_rsi):+.3f})")

    res = {"plan": "PREREGISTRATION 6g", "cohort": a.cohort,
           "n_pairs": len(common), "genes": {"D1": n1, "D2": n2},
           "positive_control": {"d1_tumor": round(float(d1[is_t].mean()), 4),
                                "d1_adjacent": round(float(d1[~is_t].mean()), 4),
                                "p": p_d1}}

    b_unadj = float(np.median(d_rsi)) if d_rsi is not None else float("nan")
    for tag, dd in ((("D1", d_d1), ("D2", d_d2)) if d_rsi is not None else ()):
        X = np.column_stack([np.ones_like(dd), dd])
        fit = ols_ci(d_rsi, X, ["intercept", "slope"])
        a_int = fit["intercept"]["beta"]
        ratio = a_int / b_unadj
        print(f"      adjusted for {tag}: intercept = {a_int:+.3f} "
              f"({fit['intercept']['ci95'][0]:+.3f} to "
              f"{fit['intercept']['ci95'][1]:+.3f})  "
              f"P = {fit['intercept']['p']:.3g}   "
              f"retains {ratio*100:.0f}% of the unadjusted shift")
        res[f"paired_adjusted_{tag}"] = {**fit, "retained_fraction":
                                         round(float(ratio), 4)}

    # 3. supporting unpaired model
    print(f"\n  [3] supporting — unpaired OLS on all {len(r)} samples")
    X = np.column_stack([np.ones(len(r)), is_t.astype(float), d1.to_numpy()])
    fit = ols_ci(r, X, ["intercept", "tumor", "D1"])
    print(f"      tumor coefficient = {fit['tumor']['beta']:+.3f} "
          f"({fit['tumor']['ci95'][0]:+.3f} to {fit['tumor']['ci95'][1]:+.3f})"
          f"  P = {fit['tumor']['p']:.3g}")
    res["unpaired_adjusted_D1"] = fit

    # 4. module level
    print("\n  [4] module level — the claim under replication (PREREG 6h)")
    tt = is_t_all.to_numpy().astype(float)
    for mod in ("reduction", "drain"):
        if d_rsi is not None:
            dm, _ = pairs(mod)
            Xm = np.column_stack([np.ones_like(d_d1), d_d1])
            f2 = ols_ci(dm, Xm, ["intercept", "slope"])
            un, ad = float(dm.mean()), f2["intercept"]["beta"]
            ci = f2["intercept"]["ci95"]; pv = f2["intercept"]["p"]
            kind = "paired"
        else:
            y = df[mod].to_numpy(float)
            f_un = ols_ci(y, np.column_stack([np.ones_like(tt), tt]),
                          ["intercept", "tumor"])
            f2 = ols_ci(y, np.column_stack([np.ones_like(tt), tt,
                                            df["d1"].to_numpy(float)]),
                        ["intercept", "tumor", "D1"])
            un, ad = f_un["tumor"]["beta"], f2["tumor"]["beta"]
            ci = f2["tumor"]["ci95"]; pv = f2["tumor"]["p"]
            f2 = {"unadjusted": f_un, "adjusted": f2}
            kind = "unpaired"
        ret = abs(ad / un) if un else float("nan")
        print(f"      {mod:<10} [{kind}] unadjusted {un:+.3f} -> adjusted "
              f"{ad:+.3f} ({ci[0]:+.3f} to {ci[1]:+.3f})  P = {pv:.3g}   "
              f"retains {ret*100:.0f}%")
        res[f"module_{mod}"] = {"kind": kind,
                                "unadjusted_mean_delta": round(float(un), 4),
                                "adjusted": round(float(ad), 4),
                                "ci95": ci, "p": pv,
                                "retained_fraction": round(float(ret), 4),
                                "fit": f2}

    # 5. sensitivity: the same module models with D2 (PREREG 6g step 5,
    #    "repeat steps 2 and 4 with D2"). D2 is expected to OVER-adjust.
    print("\n  [5] sensitivity — the same module models adjusted for D2")
    for mod in ("reduction", "drain"):
        if d_rsi is not None:
            dm, _ = pairs(mod)
            f3 = ols_ci(dm, np.column_stack([np.ones_like(d_d2), d_d2]),
                        ["intercept", "slope"])
            un, ad = float(dm.mean()), f3["intercept"]["beta"]
            ci, pv = f3["intercept"]["ci95"], f3["intercept"]["p"]
        else:
            y = df[mod].to_numpy(float)
            f_un = ols_ci(y, np.column_stack([np.ones_like(tt), tt]),
                          ["intercept", "tumor"])
            f3 = ols_ci(y, np.column_stack([np.ones_like(tt), tt,
                                            df["d2"].to_numpy(float)]),
                        ["intercept", "tumor", "D2"])
            un, ad = f_un["tumor"]["beta"], f3["tumor"]["beta"]
            ci, pv = f3["tumor"]["ci95"], f3["tumor"]["p"]
        ret2 = abs(ad / un) if un else float("nan")
        print(f"      {mod:<10} D2-adjusted {ad:+.3f} "
              f"({ci[0]:+.3f} to {ci[1]:+.3f})  P = {pv:.3g}   "
              f"retains {ret2*100:.0f}%")
        res[f"module_{mod}_D2"] = {"adjusted": round(float(ad), 4),
                                   "ci95": ci, "p": pv,
                                   "retained_fraction": round(float(ret2), 4),
                                   "fit": f3}

    # verdict per the registered decision rule
    if a.unpaired_only:
        red, dra = res["module_reduction"], res["module_drain"]
        present = red["unadjusted_mean_delta"] > 0 > dra["unadjusted_mean_delta"]
        if not present:
            verdict = ("UNINFORMATIVE — the unadjusted module shifts are not "
                       "present in the registered directions in this cohort "
                       "(PREREG 6h); counts as neither support nor refutation")
        elif (red["ci95"][0] > 0 or red["ci95"][1] < 0) and \
             red["retained_fraction"] > dra["retained_fraction"]:
            verdict = "REPLICATED — reduction survives adjustment better than drain"
        else:
            verdict = "NOT REPLICATED under the 6h rule"
        res["verdict"] = verdict
        print(f"\n  VERDICT (registered rule 6h): {verdict}")
        path = os.path.join(a.outdir,
                            f"{a.tag or a.cohort}_differentiation.json")
        with open(path, "w") as fh:
            json.dump(res, fh, indent=2, default=float)
        print(f"\n  wrote -> {path}")
        return 0

    # PREREG 6h module rule, evaluated on the paired module intercepts.
    # This is the registered 6h rule applied to the registered 6g analysis
    # (see 6i); it is reported alongside the 6g composite verdict, not
    # instead of it.
    red, dra = res["module_reduction"], res["module_drain"]
    present = red["unadjusted_mean_delta"] > 0 > dra["unadjusted_mean_delta"]
    if not present:
        v6h = ("UNINFORMATIVE — the unadjusted module shifts are not present "
               "in the registered directions in this cohort (PREREG 6h)")
    elif (red["ci95"][0] > 0 or red["ci95"][1] < 0) and \
            red["retained_fraction"] > dra["retained_fraction"]:
        v6h = ("REPLICATED — reduction survives adjustment better than drain "
               "(PREREG 6h)")
    else:
        v6h = "NOT REPLICATED under the 6h rule"
    res["verdict_6h_module_split"] = v6h
    print(f"\n  VERDICT (registered rule 6h, module split): {v6h}")

    a1 = res["paired_adjusted_D1"]["intercept"]
    frac = res["paired_adjusted_D1"]["retained_fraction"]
    if a1["p"] < 0.05 and frac >= 0.5:
        verdict = ("NOT explained by dedifferentiation — PHASE E section 3 "
                   "stands as written")
    elif a1["p"] < 0.05:
        verdict = ("PARTLY explained by dedifferentiation — both numbers must "
                   "be reported")
    else:
        verdict = ("ATTRIBUTABLE to loss of hepatocyte identity — PHASE E "
                   "section 3 must be rewritten")
    res["verdict"] = verdict
    print(f"\n  VERDICT (registered rule 6g): {verdict}")

    path = os.path.join(a.outdir,
                        f"{a.tag or a.cohort}_differentiation.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(f"\n  wrote -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
