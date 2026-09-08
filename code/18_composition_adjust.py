#!/usr/bin/env python3
"""18_composition_adjust.py -- PREREGISTRATION 6j.

Is the tumor-vs-adjacent RSI rise a change in WHICH CELLS are in the sample,
rather than a change in what the hepatocytes are doing? 6g answered the
hepatocyte-identity question; this answers the cell-composition question, and
the JOINT model (dRSI ~ dD1 + dC1) is the test 6j names as decisive.

Index values are read as written by 02_compute_rsi.py. Nothing is renormalized.
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
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "diffadj", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "12_differentiation_adjust.py"))
_da = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_da)
ols_ci, D1, D2_EXTRA = _da.ols_ci, _da.D1, _da.D2_EXTRA

# --- the covariate, exactly as fixed in 6j ---------------------------------
C1_IMMUNE = ["PTPRC", "CD53", "LAPTM5", "CD3E", "CD2", "CD68", "AIF1",
             "ITGAM", "LCP1"]
C1_ENDO = ["PECAM1", "VWF", "CDH5", "ENG", "CLEC4G"]
C1_STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PDGFRB"]
C1 = C1_IMMUNE + C1_ENDO + C1_STROMA
C2_EXTRA = ["KRT19", "KRT7", "SOX9", "EPCAM"]

assert not (set(C1) | set(C2_EXTRA)) & (
    set(cfg.MODULE_SUPPLY) | set(cfg.MODULE_REDUCTION) | set(cfg.MODULE_DRAIN)), \
    "6j requires the composition covariate to share no gene with the RSI panel"
assert not (set(C1) | set(C2_EXTRA)) & (set(D1) | set(D2_EXTRA)), \
    "6j requires the composition covariate to share no gene with D1/D2"


def zmean(mat, genes, label=""):
    g = [x for x in genes if x in mat.index]
    missing = sorted(set(genes) - set(g))
    if missing:
        print(f"    !! absent from the platform ({label}): {missing}")
    if not g:
        raise SystemExit(f"    !! none of the {label} genes are present — "
                         "the covariate cannot be built; recorded as unmet")
    z = mat.loc[g].sub(mat.loc[g].mean(axis=1), axis=0) \
                  .div(mat.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0), len(g), missing


def vif(X):
    """VIF for each non-intercept column of a design matrix."""
    out = []
    for j in range(1, X.shape[1]):
        y = X[:, j]
        others = np.delete(X, j, axis=1)
        beta, *_ = np.linalg.lstsq(others, y, rcond=None)
        resid = y - others @ beta
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1 - float(resid @ resid) / ss_tot if ss_tot > 0 else 0.0
        out.append(float("inf") if r2 >= 1 else 1.0 / (1.0 - r2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--comp-source", default=None,
                    help="TSV of composition genes x samples; default is to "
                         "read them from the cohort's own expr_log")
    ap.add_argument("--adjacent-match", default="adjacent")
    ap.add_argument("--tumor-match", default=None)
    ap.add_argument("--unpaired-only", action="store_true")
    ap.add_argument("--unadjusted-median", type=float, default=None,
                    help="b in the 6j decision rule; defaults to this "
                         "cohort's own paired median")
    ap.add_argument("--d1-source", default=None,
                    help="TSV of the D1 identity genes x samples, for cohorts "
                         "whose expr_log carries only the RSI panel")
    ap.add_argument("--tag", default=None)
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
        if (~keep).sum():
            print(f"  restricted to the two tissue groups; dropped "
                  f"{int((~keep).sum())} samples")
        rsi, ph, expr, tis = (rsi[keep.values], ph[keep.values],
                              expr.loc[:, keep.values], tis[keep.values])
    is_t = (~tis.str.contains(a.adjacent_match.lower())).to_numpy()

    comp = expr
    if a.comp_source:
        comp = pd.read_csv(a.comp_source, sep="\t", index_col=0)
        missing_cols = [c for c in rsi.index if c not in comp.columns]
        assert not missing_cols, f"composition source lacks {missing_cols[:4]}"
        comp = comp[rsi.index]

    print(f"=== {a.cohort} — composition-adjusted H2 (PREREG 6j) ===")
    print(f"  {int(is_t.sum())} tumor / {int((~is_t).sum())} adjacent")
    c1, n1, m1 = zmean(comp, C1, "C1")
    c2, n2, m2 = zmean(comp, C1 + C2_EXTRA, "C2")
    d1src = expr
    if a.d1_source:
        d1src = pd.read_csv(a.d1_source, sep="\t", index_col=0)
        assert list(d1src.columns) == list(rsi.index) or \
            set(rsi.index) <= set(d1src.columns), "D1 source lacks samples"
        d1src = d1src[rsi.index]
    d1, nd1, _ = zmean(d1src, D1, "D1")
    print(f"  C1 non-parenchymal content: {n1}/{len(C1)} genes")
    print(f"  C2 (+ cholangiocyte/progenitor): {n2}/{len(C1 + C2_EXTRA)} genes")
    print(f"  D1 hepatocyte identity: {nd1}/{len(D1)} genes")

    r = rsi["RSI"].to_numpy(float)
    res = {"plan": "PREREGISTRATION 6j", "cohort": a.cohort,
           "n_tumor": int(is_t.sum()), "n_adjacent": int((~is_t).sum()),
           "genes": {"C1": n1, "C1_missing": m1, "C2": n2, "C2_missing": m2,
                     "D1": nd1}}

    # ---- 1. positive control, direction reported not predicted -------------
    cv = c1.to_numpy()
    _, _, p_c1 = sl.mannwhitney_u(cv[is_t], cv[~is_t])
    d_c1_mean = float(cv[is_t].mean() - cv[~is_t].mean())
    direction = "HIGHER in tumor" if d_c1_mean > 0 else "LOWER in tumor"
    print(f"\n  [1] positive control — C1 tumor {cv[is_t].mean():+.3f} vs "
          f"adjacent {cv[~is_t].mean():+.3f}  ({direction})  P = {p_c1:.3g}")
    present = p_c1 < 0.05
    res["positive_control"] = {"c1_tumor": round(float(cv[is_t].mean()), 4),
                               "c1_adjacent": round(float(cv[~is_t].mean()), 4),
                               "delta": round(d_c1_mean, 4),
                               "direction": direction, "p": p_c1,
                               "confound_present": bool(present)}
    if not present:
        print("      -> the composition difference is ABSENT in this cohort "
              "(6j branch 1).\n"
              "         The objection is answered by that fact; the adjustment "
              "below is reported\n         as unnecessary here, not as a passed test.")

    # ---- module scores -----------------------------------------------------
    mods = {}
    for name, genes in (("reduction", cfg.MODULE_REDUCTION),
                        ("drain", cfg.MODULE_DRAIN)):
        s, _, _ = zmean(expr, genes, name)
        mods[name] = s.to_numpy()

    pid = ph["patient_id"].astype(str) if "patient_id" in ph.columns \
        else pd.Series(ph.index, index=ph.index)
    df = pd.DataFrame({"pid": pid.values, "t": is_t, "rsi": r,
                       "d1": d1.to_numpy(), "c1": cv, "c2": c2.to_numpy(),
                       "reduction": mods["reduction"], "drain": mods["drain"]},
                      index=rsi.index)

    def pairs(col):
        tt = df[df.t].set_index("pid")[col]; tt = tt[~tt.index.duplicated()]
        nn = df[~df.t].set_index("pid")[col]; nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[common] - nn.loc[common]).to_numpy(float), common

    paired = not a.unpaired_only
    if paired:
        d_rsi, common = pairs("rsi")
        d_d1, _ = pairs("d1"); d_c1, _ = pairs("c1"); d_c2, _ = pairs("c2")
        b = a.unadjusted_median if a.unadjusted_median is not None \
            else float(np.median(d_rsi))
        print(f"\n  paired over {len(common)} pairs; unadjusted mean ΔRSI "
              f"{d_rsi.mean():+.3f} (median {np.median(d_rsi):+.3f})")
        res["n_pairs"] = len(common)
        res["unadjusted_paired"] = {"mean": round(float(d_rsi.mean()), 4),
                                    "median": round(float(np.median(d_rsi)), 4)}
    else:
        print("\n  this series declares no patient ID (PREREG 6h/6j): the "
              "unpaired OLS is used and the\n  tissue coefficient replaces the "
              "intercept throughout.")
        res["n_pairs"] = 0
        tt = is_t.astype(float)
        f_un = ols_ci(r, np.column_stack([np.ones_like(tt), tt]),
                      ["intercept", "tumor"])
        b = f_un["tumor"]["beta"]
        res["unadjusted_unpaired_tumor_coef"] = f_un["tumor"]
        print(f"  unadjusted tumor coefficient {b:+.3f}")

    def fit(outcome, covs, names):
        """covs: list of covariate arrays (already paired-differenced or raw)."""
        if paired:
            y, _ = pairs(outcome)
            X = np.column_stack([np.ones(len(y))] + covs)
            return ols_ci(y, X, ["intercept"] + names), "intercept", y
        y = df[outcome].to_numpy(float)
        tt = is_t.astype(float)
        X = np.column_stack([np.ones(len(y)), tt] + covs)
        return ols_ci(y, X, ["intercept", "tumor"] + names), "tumor", y

    cov = {"C1": (d_c1 if paired else df["c1"].to_numpy(float)),
           "C2": (d_c2 if paired else df["c2"].to_numpy(float)),
           "D1": (d_d1 if paired else df["d1"].to_numpy(float))}

    def report(label, outcome, covs, names, unadj):
        f, key, _ = fit(outcome, [cov[c] for c in covs], names)
        est = f[key]["beta"]; ci = f[key]["ci95"]; pv = f[key]["p"]
        ret = abs(est / unadj) if unadj else float("nan")
        print(f"      {label:<34} {est:+.3f} ({ci[0]:+.3f} to {ci[1]:+.3f})"
              f"  P = {pv:.3g}   retains {ret*100:.0f}%")
        return {"estimate": round(float(est), 4), "ci95": ci, "p": pv,
                "retained_fraction": round(float(ret), 4), "fit": f}

    # ---- 2. composition only ----------------------------------------------
    print("\n  [2] composition only — ΔRSI adjusted for ΔC1")
    res["composite_C1"] = report("adjusted for C1", "rsi", ["C1"], ["C1"], b)

    # ---- 3. JOINT model — the decisive test --------------------------------
    print("\n  [3] JOINT — ΔRSI adjusted for ΔD1 + ΔC1  (6j decisive test)")
    res["composite_joint"] = report("adjusted for D1 + C1", "rsi",
                                    ["D1", "C1"], ["D1", "C1"], b)

    # ---- 4. module level ---------------------------------------------------
    print("\n  [4] module level")
    for mod in ("reduction", "drain"):
        if paired:
            un = float(pairs(mod)[0].mean())
        else:
            tt = is_t.astype(float)
            un = ols_ci(df[mod].to_numpy(float),
                        np.column_stack([np.ones_like(tt), tt]),
                        ["intercept", "tumor"])["tumor"]["beta"]
        res[f"module_{mod}_unadjusted"] = round(un, 4)
        res[f"module_{mod}_C1"] = report(f"{mod} · C1", mod, ["C1"], ["C1"], un)
        res[f"module_{mod}_joint"] = report(f"{mod} · D1 + C1", mod,
                                            ["D1", "C1"], ["D1", "C1"], un)

    # ---- 5. C2 sensitivity -------------------------------------------------
    print("\n  [5] sensitivity — C2 in place of C1 (deliberate over-adjustment)")
    res["composite_joint_C2"] = report("ΔRSI · D1 + C2", "rsi",
                                       ["D1", "C2"], ["D1", "C2"], b)
    for mod in ("reduction", "drain"):
        un = res[f"module_{mod}_unadjusted"]
        res[f"module_{mod}_joint_C2"] = report(f"{mod} · D1 + C2", mod,
                                               ["D1", "C2"], ["D1", "C2"], un)

    # ---- 6. collinearity, reported either way ------------------------------
    dd1, dc1 = cov["D1"], cov["C1"]
    rho = float(np.corrcoef(dd1, dc1)[0, 1])
    X = np.column_stack([np.ones(len(dd1)), dd1, dc1])
    v = vif(X)
    print(f"\n  [6] collinearity — corr(ΔD1, ΔC1) = {rho:+.3f}   "
          f"VIF D1 {v[0]:.2f}, C1 {v[1]:.2f}")
    unstable = max(v) >= 5
    if unstable:
        print("      !! VIF >= 5 — the joint intercept is UNSTABLE and the "
              "estimate above must not\n         be quoted as if it were precise.")
    res["collinearity"] = {"corr_dD1_dC1": round(rho, 4),
                           "vif_D1": round(v[0], 3), "vif_C1": round(v[1], 3),
                           "unstable": bool(unstable)}

    # ---- precondition flag, added 2026-08-26 AFTER seeing GSE164760 --------
    # Not a change to the rule and not a verdict override: the verdict strings
    # below are exactly those 6j fixed. This only SURFACES a precondition the
    # composite rule assumes -- that there is an unadjusted effect to explain.
    # 6h already had this branch ("uninformative"); 6j's composite rule does
    # not, and GSE164760 is the case that exposes it. Recorded in 6j.
    if paired:
        n_ = len(d_rsi)
        t_ = d_rsi.mean() / (d_rsi.std(ddof=1) / np.sqrt(n_))
        p_un = float(sl._t_two_sided(float(t_), n_ - 1))
    else:
        p_un = float(res["unadjusted_unpaired_tumor_coef"]["p"])
    no_effect = p_un >= 0.05
    res["unadjusted_effect_p"] = p_un
    res["no_unadjusted_effect"] = bool(no_effect)
    if no_effect:
        print(f"\n  !! PRECONDITION: the UNADJUSTED tumor effect is not "
              f"distinguishable from zero here (P = {p_un:.3g}).\n"
              "     There is no effect for an adjustment to explain, so every "
              "retention ratio below\n     is a ratio to ~0 and the composite "
              "verdict string is mechanical, not interpretable.\n"
              "     6h had already declared this cohort UNINFORMATIVE by a rule "
              "fixed in advance, and\n     6j states it is supporting context, "
              "not a replication.")

    # ---- verdict, by the rule fixed in 6j ----------------------------------
    j = res["composite_joint"]
    frac = j["retained_fraction"]
    sig = (j["ci95"][0] > 0) and j["p"] < 0.05
    if not present:
        verdict = ("CONFOUND ABSENT — C1 does not differ between tumor and "
                   "adjacent in this cohort; the composition objection does not "
                   "arise here (6j branch 1)")
    elif sig and frac >= 0.5:
        verdict = ("NOT explained by identity loss or composition — PHASE E "
                   "§3/§3b stand and are strengthened")
    elif sig:
        verdict = ("PARTLY explained by the two together — both numbers must "
                   "be reported")
    else:
        verdict = ("ATTRIBUTABLE to identity loss plus composition — PHASE E "
                   "§3/§3b must be rewritten")
    res["verdict_composite"] = verdict
    print(f"\n  VERDICT (6j composite rule): {verdict}")

    red, dra = res["module_reduction_joint"], res["module_drain_joint"]
    red_ci_excl = red["ci95"][0] > 0 or red["ci95"][1] < 0
    if red_ci_excl and red["retained_fraction"] > dra["retained_fraction"]:
        mv = ("SURVIVES composition adjustment — reduction retains more than "
              "drain in the joint model")
    else:
        mv = "DOES NOT survive composition adjustment under the 6j module rule"
    res["verdict_module_split"] = mv
    print(f"  VERDICT (6j module rule):    {mv}")

    # Knife-edge flag, added 2026-08-26 AFTER seeing GSE76427 — recorded in 6j.
    # Like the precondition flag above this SURFACES a fact; it changes no rule
    # and no verdict string. The 6j module rule is a strict inequality between
    # two retentions, so it can be decided by a margin far smaller than the
    # estimates' own uncertainty. When it is, that is said plainly.
    gap = red["retained_fraction"] - abs(dra["retained_fraction"])
    res["module_retention_gap"] = round(float(gap), 5)
    if present and abs(gap) < 0.05:
        print(f"      !! KNIFE-EDGE: the two retentions differ by {gap:+.4f} "
              f"({red['retained_fraction']*100:.2f}% vs "
              f"{abs(dra['retained_fraction'])*100:.2f}%).\n"
              "         The rule is a strict inequality and returns a verdict, "
              "but a margin this small is\n         far inside the estimates' "
              "own uncertainty. It must not be leaned on.")
        res["module_verdict_knife_edge"] = True

    if unstable:
        res["verdict_caveat"] = ("joint model collinear (VIF >= 5); the "
                                 "estimates are unstable and the verdict is "
                                 "reported with that caveat")

    path = os.path.join(a.outdir, f"{a.tag or a.cohort}_composition.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(f"\n  wrote -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
