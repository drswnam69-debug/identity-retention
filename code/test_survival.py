#!/usr/bin/env python3
"""test_survival.py -- validate survival.py against lifelines and scipy.

survival.py must run with only numpy and pandas; this checks the pure-numpy
implementations are correct where a reference is available, and by simulation
where one is not.
"""
from __future__ import annotations
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import survival as sv
import stats_lite as sl

try:
    import pandas as pd
    from statsmodels.duration.hazard_regression import PHReg
    from statsmodels.duration.survfunc import SurvfuncRight, survdiff
    HAVE_REF = True
except Exception:
    HAVE_REF = False
try:
    from scipy import stats as sp
except Exception:
    sp = None

rng = np.random.default_rng(76427)
checks = []
def chk(name, ok, detail=""): checks.append((name, bool(ok), detail))


def simulate(n=250, betas=(0.6, -0.4), ties=False, seed=0):
    r = np.random.default_rng(seed)
    X = np.column_stack([r.normal(size=n), r.integers(0, 2, n).astype(float)])
    lp = X @ np.array(betas)
    t = r.exponential(1 / np.exp(lp))
    c = r.exponential(1.2, n)
    time = np.minimum(t, c)
    event = (t <= c).astype(float)
    if ties:
        time = np.round(time, 1) + 0.1     # force many shared event times
    return time, event, X


# ---- Cox: coefficients, SEs, p-values ---------------------------------
for label, ties in (("no ties", False), ("heavy ties", True)):
    time, event, X = simulate(ties=ties, seed=1 if not ties else 2)
    got = sv.cox_ph(time, event, X, names=["x1", "x2"])
    if HAVE_REF:
        ref = PHReg(time, X, status=event, ties="breslow").fit()
        chk(f"Cox coefficients match statsmodels ({label})",
            np.allclose(got["coef"], ref.params, atol=1e-6),
            f"max diff {np.max(np.abs(got['coef'] - ref.params)):.2e}")
        chk(f"Cox standard errors match statsmodels ({label})",
            np.allclose(got["se"], ref.bse, rtol=1e-4),
            f"max rel diff {np.max(np.abs(got['se'] / ref.bse - 1)):.2e}")
        chk(f"Cox hazard ratios match statsmodels ({label})",
            np.allclose(got["hr"], np.exp(ref.params), atol=1e-6))
        chk(f"Cox p-values match statsmodels ({label})",
            np.allclose(got["p"], ref.pvalues, atol=1e-6),
            f"max diff {np.max(np.abs(got['p'] - ref.pvalues)):.2e}")

# recovery of the true coefficients
time, event, X = simulate(n=4000, betas=(0.6, -0.4), seed=7)
got = sv.cox_ph(time, event, X)
chk("Cox recovers the simulated coefficients",
    abs(got["coef"][0] - 0.6) < 0.08 and abs(got["coef"][1] + 0.4) < 0.12,
    f"got {got['coef'][0]:.3f}, {got['coef'][1]:.3f} for 0.6, -0.4")

# null calibration
rej = 0
for i in range(300):
    r = np.random.default_rng(500 + i)
    x = r.normal(size=150)
    t = r.exponential(1.0, 150); c = r.exponential(1.2, 150)
    res = sv.cox_ph(np.minimum(t, c), (t <= c).astype(float), x[:, None])
    rej += res["p"][0] < 0.05
chk("Cox null type-I error is near 5%", 0.02 <= rej / 300 <= 0.09,
    f"{rej / 300:.3f} over 300 null datasets")

# ---- Kaplan-Meier ------------------------------------------------------
time, event, X = simulate(n=200, seed=3)
t_km, s_km, _ = sv.kaplan_meier(time, event)
if HAVE_REF:
    sf = SurvfuncRight(time, event)
    ref = np.array([sf.surv_prob[np.searchsorted(sf.surv_times, tt, "right") - 1]
                    for tt in t_km])
    chk("Kaplan-Meier curve matches statsmodels",
        np.allclose(s_km, ref, atol=1e-9),
        f"max diff {np.max(np.abs(s_km - ref)):.2e}")

# ---- log-rank ----------------------------------------------------------
g = (X[:, 1] > 0.5).astype(int)
lr = sv.logrank(time, event, g)
if HAVE_REF:
    chi2_ref, p_ref = survdiff(time, event, g)
    chk("log-rank p matches statsmodels (2 groups)",
        abs(lr["p"] - p_ref) < 1e-6,
        f"ours {lr['p']:.6g} vs {p_ref:.6g}")

g3 = rng.integers(0, 3, len(time))
lr3 = sv.logrank(time, event, g3)
chk("log-rank runs for 3 groups and returns a valid p",
    0 <= lr3["p"] <= 1 and lr3["df"] == 2, f"p={lr3['p']:.3g}")

# chi-square tail against scipy
if sp is not None:
    worst = max(abs(sv._chi2_sf(x, df) - sp.chi2.sf(x, df))
                for x in (0.5, 2.0, 6.0, 15.0, 40.0) for df in (1, 2, 5, 12))
    chk("chi-square upper tail matches scipy", worst < 1e-12,
        f"max abs err {worst:.2e}")

# ---- proportional hazards check ---------------------------------------
time, event, X = simulate(n=600, betas=(0.7, 0.0), seed=11)
res = sv.cox_ph(time, event, X[:, [0]])
ph = sv.schoenfeld_ph_test(time, event, X[:, 0], float(res["coef"][0]))
chk("PH test does not flag a proportional-hazards dataset", ph["p"] > 0.01,
    f"p={ph['p']:.3g}")

# a deliberately non-proportional dataset must be flagged: crossing hazards,
# built from Weibull shapes that differ between the groups
r = np.random.default_rng(21)
n = 1200
x = r.integers(0, 2, n).astype(float)
u = r.uniform(size=n)
# shape 0.6 gives a falling hazard, shape 2.5 a rising one; same median
t = np.where(x > 0.5, (-np.log(u)) ** (1 / 2.5), (-np.log(u)) ** (1 / 0.6))
c = r.exponential(3.0, n)
tt, ee = np.minimum(t, c), (t <= c).astype(float)
res2 = sv.cox_ph(tt, ee, x[:, None])
ph2 = sv.schoenfeld_ph_test(tt, ee, x, float(res2["coef"][0]))
chk("PH test flags a non-proportional dataset", ph2["p"] < 0.05,
    f"p={ph2['p']:.3g}")

w = max(len(n) for n, _, _ in checks)
ok_all = True
print(f"statsmodels reference: {HAVE_REF} | scipy: {sp is not None}\n")
for name, ok, detail in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{w}}  {detail}")
    ok_all &= ok
print(f"\n{'ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED'}")
raise SystemExit(0 if ok_all else 1)
