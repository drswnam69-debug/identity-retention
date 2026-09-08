#!/usr/bin/env python3
"""19_verify_composition.py -- verification of the 6j joint-model machinery.

Five checks on synthetic data whose truth is known by construction. None can
pass by accident, and two of them are designed to FAIL loudly if the engine
ever starts reporting a collinear joint fit as if it were precise.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util

spec = importlib.util.spec_from_file_location(
    "comp", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "18_composition_adjust.py"))
comp = importlib.util.module_from_spec(spec)
sys.argv = [sys.argv[0], "--cohort", "__none__"]
try:
    spec.loader.exec_module(comp)
except SystemExit:
    pass
ols_ci, vif, C1, C2_EXTRA = comp.ols_ci, comp.vif, comp.C1, comp.C2_EXTRA

rng = np.random.default_rng(20260826)
ok = True


def check(label, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {label:<56} {detail}")


def joint(y, x1, x2):
    X = np.column_stack([np.ones(len(y)), x1, x2])
    return ols_ci(y, X, ["intercept", "D1", "C1"])


print("=== verification of the composition adjustment (PREREG 6j) ===")

N = 400

# 0 -- the covariate is disjoint from everything it must be disjoint from
import rsi_config as cfg
D1 = comp.D1
panel = set(cfg.MODULE_SUPPLY) | set(cfg.MODULE_REDUCTION) | set(cfg.MODULE_DRAIN)
check("C1/C2 share no gene with the RSI panel",
      not (set(C1) | set(C2_EXTRA)) & panel)
check("C1/C2 share no gene with D1/D2",
      not (set(C1) | set(C2_EXTRA)) & (set(D1) | set(comp.D2_EXTRA)))
check("C1 is the 21 genes registered in 6j", len(C1) == 21, f"n={len(C1)}")

# 1 -- shift entirely driven by composition
x1 = rng.normal(0, 1, N)               # dD1
x2 = rng.normal(0, 1, N)               # dC1
y = 1.4 * x2
f = joint(y, x1, x2)
check("composition-driven shift gives joint intercept 0",
      abs(f["intercept"]["beta"]) < 1e-6 and abs(f["C1"]["beta"] - 1.4) < 1e-6,
      f"intercept={f['intercept']['beta']:+.2e}, C1 slope={f['C1']['beta']:+.3f}")

# 2 -- shift driven by identity AND composition, plus a real residual rise
x1 = rng.normal(0, 1, N)
x2 = rng.normal(0, 1, N)
y = 0.60 - 0.9 * x1 + 0.7 * x2 + rng.normal(0, 0.25, N)
f = joint(y, x1, x2)
check("residual rise is recovered under two true confounders",
      abs(f["intercept"]["beta"] - 0.60) < 0.05 and f["intercept"]["ci95"][0] > 0,
      f"intercept={f['intercept']['beta']:+.3f} (true +0.600)")

# 3 -- single-covariate model is BIASED when the omitted covariate matters;
#      the joint model is the one that recovers the truth. This is the whole
#      argument for step 3 being decisive, so it is asserted, not assumed.
x1 = rng.normal(0, 1, N)
x2 = 0.8 * x1 + rng.normal(0, 0.6, N)          # correlated confounders
y = 0.50 - 1.0 * x1 + 0.8 * x2 + rng.normal(0, 0.2, N)
only_d1 = ols_ci(y, np.column_stack([np.ones(N), x1]), ["intercept", "D1"])
fj = joint(y, x1, x2)
check("joint model beats the single-covariate model here",
      abs(fj["intercept"]["beta"] - 0.50) < abs(only_d1["intercept"]["beta"] - 0.50),
      f"joint {fj['intercept']['beta']:+.3f} vs D1-only "
      f"{only_d1['intercept']['beta']:+.3f} (true +0.500)")

# 4 -- collinearity must be DETECTED, not silently absorbed
x1 = rng.normal(0, 1, N)
x2 = x1 + rng.normal(0, 0.05, N)               # near-duplicate covariates
X = np.column_stack([np.ones(N), x1, x2])
v = vif(X)
check("near-duplicate covariates are flagged by VIF", max(v) >= 5,
      f"VIF = {v[0]:.1f}, {v[1]:.1f}  (rule fires at 5)")

# 5 -- and independent covariates are NOT flagged
x1 = rng.normal(0, 1, N)
x2 = rng.normal(0, 1, N)
v = vif(np.column_stack([np.ones(N), x1, x2]))
check("independent covariates are not falsely flagged", max(v) < 1.5,
      f"VIF = {v[0]:.2f}, {v[1]:.2f}")

# 6 -- absent-confound branch: if C1 does not differ, adjusting changes nothing
x1 = rng.normal(0, 1, N)
x2 = rng.normal(0, 1, N)
y = 0.9 + rng.normal(0, 0.3, N)                # rise unrelated to either
f = joint(y, x1, x2)
check("with no true confounding the estimate is unmoved",
      abs(f["intercept"]["beta"] - 0.9) < 0.05,
      f"intercept={f['intercept']['beta']:+.3f} (true +0.900)")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
raise SystemExit(0 if ok else 1)
