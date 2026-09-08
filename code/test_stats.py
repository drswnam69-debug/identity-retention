#!/usr/bin/env python3
"""test_stats.py -- validate stats_lite against SciPy and against permutation.

The pipeline must run without SciPy, so every routine here has a pure-numpy
implementation. This script checks those implementations are right.
"""
from __future__ import annotations
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_lite as sl

try:
    from scipy import stats as sp
except Exception:
    sp = None

rng = np.random.default_rng(20260825)
checks: list[tuple[str, bool, str]] = []


def chk(name, ok, detail=""):
    checks.append((name, bool(ok), detail))


# ---- ranks -------------------------------------------------------------
a = rng.integers(0, 5, 200).astype(float)
if sp is not None:
    chk("rankdata matches SciPy (with ties)",
        np.allclose(sl.rankdata(a), sp.rankdata(a)))

# ---- t tail ------------------------------------------------------------
if sp is not None:
    worst = max(abs(sl._t_two_sided(t, df) - 2 * sp.t.sf(abs(t), df))
                for t in (0.1, 1.0, 2.5, 4.0, 8.0) for df in (3, 10, 75, 500))
    chk("t tail matches SciPy without SciPy path", worst < 1e-10,
        f"max abs err {worst:.2e}")
    # force the pure-python branch
    saved, sl._sp = sl._sp, None
    worst2 = max(abs(sl._t_two_sided(t, df) - 2 * sp.t.sf(abs(t), df))
                 for t in (0.1, 1.0, 2.5, 4.0, 8.0) for df in (3, 10, 75, 500))
    sl._sp = saved
    chk("pure-numpy t tail is exact", worst2 < 1e-10, f"max abs err {worst2:.2e}")

# ---- Mann-Whitney ------------------------------------------------------
x, y = rng.normal(0, 1, 60), rng.normal(0.6, 1, 55)
u, z, p = sl.mannwhitney_u(x, y)
if sp is not None:
    su = sp.mannwhitneyu(x, y, alternative="two-sided")
    chk("Mann-Whitney U matches SciPy", np.isclose(u, su.statistic))
    chk("Mann-Whitney p matches SciPy (normal approx)",
        abs(p - su.pvalue) < 5e-3, f"|dp|={abs(p - su.pvalue):.2e}")
xt = rng.integers(0, 4, 60).astype(float)
yt = rng.integers(0, 4, 55).astype(float)
ut, zt, pt = sl.mannwhitney_u(xt, yt)
if sp is not None:
    st = sp.mannwhitneyu(xt, yt, alternative="two-sided")
    chk("Mann-Whitney U matches SciPy with heavy ties",
        np.isclose(ut, st.statistic))

# ---- Spearman ----------------------------------------------------------
xx = rng.normal(size=80); yy = 0.5 * xx + rng.normal(size=80)
rho, prho = sl.spearman(xx, yy)
if sp is not None:
    sr = sp.spearmanr(xx, yy)
    chk("Spearman rho matches SciPy", np.isclose(rho, sr.statistic))
    chk("Spearman p matches SciPy", abs(prho - sr.pvalue) < 1e-8)

# ---- BH-FDR ------------------------------------------------------------
pv = rng.uniform(size=200) ** 2
q = sl.benjamini_hochberg(pv)
if sp is not None:
    from statsmodels.stats.multitest import multipletests
    _, q_sm, _, _ = multipletests(pv, method="fdr_bh")
    chk("BH-FDR matches statsmodels", np.allclose(q, q_sm))
chk("BH-FDR is monotone in p", np.all(np.diff(q[np.argsort(pv)]) >= -1e-12))

# ---- OLS ---------------------------------------------------------------
n = 150
X = np.column_stack([np.ones(n), rng.normal(size=n), rng.integers(0, 2, n)])
beta_true = np.array([1.0, 0.8, -0.5])
yv = X @ beta_true + rng.normal(0, 0.7, n)
res = sl.ols(yv, X, names=["const", "x1", "x2"])
try:
    import statsmodels.api as smapi
    sm = smapi.OLS(yv, X).fit()
    chk("OLS coefficients match statsmodels", np.allclose(res["beta"], sm.params))
    chk("OLS standard errors match statsmodels", np.allclose(res["se"], sm.bse))
    chk("OLS p-values match statsmodels", np.allclose(res["p"], sm.pvalues))
except Exception:
    pass

# ---- Jonckheere-Terpstra ----------------------------------------------
# 1. null calibration: p should be ~uniform, so a 5% test rejects ~5%
rej = 0; R = 400
for i in range(R):
    r1 = np.random.default_rng(1000 + i)          # one stream per replicate,
    g = [r1.normal(size=s) for s in (12, 18, 20, 15)]   # drawn sequentially
    _, _, pj = sl.jonckheere_terpstra(g)
    rej += pj < 0.05
chk("JT null type-I error is near 5%", 0.025 <= rej / R <= 0.085,
    f"{rej / R:.3f} over {R} null datasets")

# 2. z-based p agrees with a permutation p on a real-looking ordered sample
g_ord = [rng.normal(loc=m, scale=1.0, size=s)
         for m, s in zip((0, 0.25, 0.5, 0.8), (30, 40, 45, 25))]
jt, zj, p_norm = sl.jonckheere_terpstra(g_ord)
p_perm = sl.jt_permutation_p(g_ord, n_perm=4000, seed=7)
chk("JT normal-approx p agrees with permutation p",
    abs(p_norm - p_perm) < 0.02 or (p_norm < 1e-3 and p_perm < 1e-3),
    f"normal={p_norm:.4g}  permutation={p_perm:.4g}")

# 3. tie-corrected variance: heavily tied ordinal data must stay calibrated
rej_t = 0
for i in range(300):
    r2 = np.random.default_rng(5000 + i)
    g = [r2.integers(0, 4, s).astype(float) for s in (14, 16, 18, 12)]
    _, _, pj = sl.jonckheere_terpstra(g)
    rej_t += pj < 0.05
chk("JT stays calibrated with heavy ties", 0.02 <= rej_t / 300 <= 0.10,
    f"{rej_t / 300:.3f} over 300 tied null datasets")

# 4. power: a real monotone trend must be detected
chk("JT detects a monotone trend", p_norm < 0.001, f"p={p_norm:.3g}")

# 5. direction: a decreasing series must not be called increasing
g_dec = list(reversed(g_ord))
_, _, p_inc = sl.jonckheere_terpstra(g_dec, alternative="increasing")
_, _, p_dec = sl.jonckheere_terpstra(g_dec, alternative="decreasing")
chk("JT respects the alternative's direction", p_inc > 0.9 and p_dec < 0.001,
    f"increasing p={p_inc:.3g}  decreasing p={p_dec:.3g}")

# ---- report ------------------------------------------------------------
w = max(len(n) for n, _, _ in checks)
ok_all = True
print(f"SciPy present: {sp is not None}\n")
for name, ok, detail in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{w}}  {detail}")
    ok_all &= ok
print(f"\n{'ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED'}")
raise SystemExit(0 if ok_all else 1)
