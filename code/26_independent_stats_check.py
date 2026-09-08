#!/usr/bin/env python3
"""26_independent_stats_check.py -- adversarial re-derivation.

The whole analysis pipeline is hand-rolled without SciPy, by design. That means
a systematic bug in stats_lite or in ols_ci would propagate to every P value and
every confidence interval in the manuscript without any internal test catching
it: the unit tests compare the code against itself.

This script recomputes the reported quantities with SciPy and statsmodels --
independent implementations, independent numerics -- and reports any
disagreement. It is a check, not part of the analysis: nothing here is fed back
into the index or into any result.
"""
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps
import statsmodels.api as sm

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
ols_ci, D1, D2_EXTRA = _da.ols_ci, _da.D1, _da.D2_EXTRA
C1 = _ca.C1

SOURCES = {
    "GSE76427": {"d1": "data/GSE76427_d1_panel.tsv",
                 "c1": "data/GSE76427_composition_panel.tsv"},
    "GSE14520": {"d1": None,
                 "c1": "data/GSE14520_composition_panel.tsv"},
}

PROBLEMS = []
CHECKS = [0]


def cmp(label, mine, theirs, rtol=1e-6, atol=1e-9):
    CHECKS[0] += 1
    ok = np.isclose(mine, theirs, rtol=rtol, atol=atol)
    if not ok:
        PROBLEMS.append((label, float(mine), float(theirs)))
        print(f"  MISMATCH {label}: hand-rolled {mine:.10g} vs reference {theirs:.10g}")
    return ok


# ---------------------------------------------------------------- primitives
def check_primitives():
    print("\n=== 1. statistical primitives against SciPy, on random data ===")
    rng = np.random.default_rng(20260831)
    for i in range(200):
        n1, n2 = rng.integers(8, 60), rng.integers(8, 60)
        x = rng.normal(size=n1)
        y = rng.normal(loc=0.4, size=n2)
        _, _, p_mine = sl.mannwhitney_u(x, y)
        p_ref = sps.mannwhitneyu(x, y, alternative="two-sided",
                                 method="asymptotic").pvalue
        cmp(f"mannwhitney[{i}]", p_mine, p_ref, rtol=1e-4)
    # rank handling with ties, which is where a hand-rolled version usually breaks
    for i in range(100):
        x = rng.integers(0, 5, size=40).astype(float)
        r_mine = sl.rankdata(x)
        r_ref = sps.rankdata(x)
        cmp(f"rankdata_ties[{i}]", np.abs(r_mine - r_ref).max(), 0.0, atol=1e-9)
    # two-sided normal tail
    for z in [0.0, 0.5, 1.0, 1.96, 2.5, 3.0, 4.0, 5.0, 6.0]:
        cmp(f"norm_two_sided(z={z})", sl.norm_two_sided(z),
            2 * sps.norm.sf(abs(z)), rtol=1e-8, atol=1e-15)
    # Student t tail, used for every CI and every regression P value
    for dof in [10, 25, 50, 51, 211, 212, 500]:
        for t in [0.5, 1.5, 2.0, 3.0, 5.0, 10.0]:
            cmp(f"t_two_sided(t={t},dof={dof})", sl._t_two_sided(t, dof),
                2 * sps.t.sf(abs(t), dof), rtol=1e-6, atol=1e-14)
    # Spearman
    for i in range(100):
        n = rng.integers(20, 200)
        a = rng.normal(size=n)
        b = 0.4 * a + rng.normal(size=n)
        rho_mine = sl.spearman(a, b)[0] if isinstance(sl.spearman(a, b), tuple) \
            else sl.spearman(a, b)
        cmp(f"spearman[{i}]", rho_mine, sps.spearmanr(a, b).statistic, rtol=1e-6)


# ------------------------------------------------------------ t critical value
def check_tcrit():
    print("\n=== 2. the t critical value behind every confidence interval ===")
    # after 6n, ols_ci uses the bisection at every dof
    for dof in [50, 51, 55, 100, 199, 200, 201, 211]:
        exact = sps.t.ppf(0.975, dof)
        used = _da._tcrit(dof)   # 6n: shortcut removed; this is now ols_ci's real path
        CHECKS[0] += 1
        err = abs(used - exact)
        flag = "" if err < 5e-3 else "   <-- CI width off by more than 0.5%"
        print(f"  dof={dof:4d}  used {used:.5f}  exact {exact:.5f}  "
              f"diff {err:+.5f}{flag}")
        if err >= 5e-3:
            PROBLEMS.append((f"tcrit(dof={dof})", used, exact))


# ------------------------------------------------------------------ regression
def zmean(mat, genes):
    g = [x for x in genes if x in mat.index]
    z = mat.loc[g].sub(mat.loc[g].mean(axis=1), axis=0) \
                  .div(mat.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0)


def paired_frame(cohort):
    e = pd.read_csv(f"results/{cohort}/expr_log.tsv.gz", sep="\t", index_col=0)
    ph = pd.read_csv(f"results/{cohort}/phenotype.tsv", sep="\t", index_col=0)
    rsi = pd.read_csv(f"results/{cohort}/rsi.tsv", sep="\t", index_col=0)
    keep = [c for c in e.columns if c in ph.index and c in rsi.index]
    e, ph = e[keep], ph.loc[keep]
    src = SOURCES[cohort]
    d1src = e if src["d1"] is None else pd.read_csv(src["d1"], sep="\t", index_col=0)[keep]
    c1src = pd.read_csv(src["c1"], sep="\t", index_col=0)[keep]
    d1 = zmean(d1src, D1)
    d2 = zmean(d1src, D1 + D2_EXTRA)
    c1 = zmean(c1src, C1)
    r = rsi.loc[keep, "RSI"].astype(float)
    red = zmean(e, cfg.MODULE_REDUCTION)
    drn = zmean(e, cfg.MODULE_DRAIN)
    is_t = ~ph["tissue"].astype(str).str.lower().str.contains(
        "adjacent|non-tumor|nontumor|normal")
    rows = []
    for p, g in pd.DataFrame({"p": ph["patient_id"], "t": is_t}).groupby("p"):
        ti, ni = g.index[g["t"]], g.index[~g["t"]]
        if len(ti) == 1 and len(ni) == 1:
            rows.append((ti[0], ni[0]))
    T = [a for a, _ in rows]
    N = [b for _, b in rows]
    return pd.DataFrame({
        "dRSI": r[T].values - r[N].values,
        "dD1": d1[T].values - d1[N].values,
        "dD2": d2[T].values - d2[N].values,
        "dC1": c1[T].values - c1[N].values,
        "dRED": red[T].values - red[N].values,
        "dDRN": drn[T].values - drn[N].values})


def check_regressions():
    print("\n=== 3. every reported intercept, CI and P against statsmodels ===")
    D = {c: json.load(open(f"results/{c}_differentiation.json"))
         for c in ("GSE76427", "GSE14520")}
    C = {c: json.load(open(f"results/{c}_composition.json"))
         for c in ("GSE76427", "GSE14520")}
    for cohort in ("GSE76427", "GSE14520"):
        df = paired_frame(cohort)
        print(f"\n  --- {cohort} ({len(df)} pairs) ---")
        models = [
            ("D1",    ["dD1"],        D[cohort]["paired_adjusted_D1"]["intercept"]),
            ("D2",    ["dD2"],        D[cohort]["paired_adjusted_D2"]["intercept"]),
            ("C1",    ["dC1"],        C[cohort]["composite_C1"]["fit"]["intercept"]),
            ("joint", ["dD1", "dC1"], C[cohort]["composite_joint"]["fit"]["intercept"]),
        ]
        for name, xs, rep in models:
            X = sm.add_constant(df[xs].to_numpy())
            fit = sm.OLS(df["dRSI"].to_numpy(), X).fit()
            ci = fit.conf_int(0.05)[0]
            print(f"    {name:6s} intercept  reported {rep['beta']:+.4f} "
                  f"[{rep['ci95'][0]:+.4f},{rep['ci95'][1]:+.4f}] P={rep['p']:.3g}"
                  f"   |  statsmodels {fit.params[0]:+.4f} "
                  f"[{ci[0]:+.4f},{ci[1]:+.4f}] P={fit.pvalues[0]:.3g}")
            cmp(f"{cohort}/{name}/beta", rep["beta"], fit.params[0], rtol=1e-3, atol=1e-4)
            cmp(f"{cohort}/{name}/p", rep["p"], fit.pvalues[0], rtol=2e-2, atol=1e-12)
            cmp(f"{cohort}/{name}/ci_lo", rep["ci95"][0], ci[0], rtol=5e-3, atol=5e-3)
            cmp(f"{cohort}/{name}/ci_hi", rep["ci95"][1], ci[1], rtol=5e-3, atol=5e-3)
        # module joint models
        for mod, col in (("reduction", "dRED"), ("drain", "dDRN")):
            rep = C[cohort][f"module_{mod}_joint"]
            X = sm.add_constant(df[["dD1", "dC1"]].to_numpy())
            fit = sm.OLS(df[col].to_numpy(), X).fit()
            cmp(f"{cohort}/module_{mod}/beta", rep["estimate"], fit.params[0],
                rtol=1e-3, atol=1e-4)
            cmp(f"{cohort}/module_{mod}/p", rep["p"], fit.pvalues[0],
                rtol=2e-2, atol=1e-12)
        # Wilcoxon signed rank on the paired differences
        w_ref = sps.wilcoxon(df["dRSI"].to_numpy(), alternative="two-sided",
                             method="approx")
        print(f"    Wilcoxon (SciPy, normal approx): P = {w_ref.pvalue:.4g}")


def main():
    check_primitives()
    check_tcrit()
    check_regressions()
    print(f"\n{'='*66}\n{CHECKS[0]} independent comparisons run.")
    if PROBLEMS:
        print(f"{len(PROBLEMS)} MISMATCHES:")
        for lab, a, b in PROBLEMS[:40]:
            print(f"  {lab}: hand-rolled {a:.10g} vs reference {b:.10g}")
    else:
        print("No mismatch: the hand-rolled pipeline agrees with SciPy/statsmodels.")
    json.dump({"n_checks": CHECKS[0],
               "n_mismatches": len(PROBLEMS),
               "mismatches": [{"what": a, "hand_rolled": b, "reference": c}
                              for a, b, c in PROBLEMS]},
              open("results/INDEPENDENT_STATS_CHECK.json", "w"), indent=1)


if __name__ == "__main__":
    main()
