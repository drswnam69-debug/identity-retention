#!/usr/bin/env python3
"""17_verify_adjustment.py -- verification of the 6g/6h adjustment machinery.

Three checks, none of which can pass by accident:

  1. If DeltaRSI is a pure function of DeltaD1, the registered intercept must be 0.
  2. If DeltaRSI is a constant shift independent of DeltaD1, the intercept must
     recover that constant.
  3. On the real GSE14520 pairs, flipping tumor/adjacent at random within
     each pair must destroy the effect. The observed intercept has to sit far
     outside that null.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv = [sys.argv[0]]
import importlib.util
spec = importlib.util.spec_from_file_location(
    "diffadj", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "12_differentiation_adjust.py"))
diffadj = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diffadj)
ols_ci = diffadj.ols_ci
D1 = diffadj.D1

rng = np.random.default_rng(20260826)
ok = True


def check(label, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {label:<58} {detail}")


print("=== verification of the differentiation adjustment (PREREG 6g/6h) ===")

# 1 -- fully explained by identity loss
d_d1 = rng.normal(0, 1, 400)
d_rsi = -1.7 * d_d1
fit = ols_ci(d_rsi, np.column_stack([np.ones(400), d_d1]),
             ["intercept", "slope"])
check("fully identity-driven shift gives intercept 0",
      abs(fit["intercept"]["beta"]) < 1e-6 and fit["slope"]["beta"] < -1.6,
      f"intercept={fit['intercept']['beta']:+.2e} slope={fit['slope']['beta']:+.2f}")

# 2 -- constant shift, independent of identity
d_d1 = rng.normal(0, 1, 400)
d_rsi = 0.85 + rng.normal(0, 0.3, 400)
fit = ols_ci(d_rsi, np.column_stack([np.ones(400), d_d1]),
             ["intercept", "slope"])
check("identity-independent shift is recovered intact",
      abs(fit["intercept"]["beta"] - 0.85) < 0.05
      and fit["intercept"]["ci95"][0] > 0,
      f"intercept={fit['intercept']['beta']:+.3f} (true +0.850)")

# 3 -- label permutation on the real cohort
base = "results/GSE14520"
if os.path.exists(f"{base}/rsi.tsv"):
    rsi = pd.read_csv(f"{base}/rsi.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(f"{base}/phenotype.tsv", sep="\t", index_col=0).loc[rsi.index]
    expr = pd.read_csv(f"{base}/expr_log.tsv.gz", sep="\t", index_col=0)[rsi.index]
    d1s, _ = diffadj.zmean(expr, D1)
    is_t = (ph["tissue"] == "HCC tumor").to_numpy()
    df = pd.DataFrame({"pid": ph["patient_id"].astype(str).values, "t": is_t,
                       "rsi": rsi["RSI"].to_numpy(float),
                       "d1": d1s.to_numpy()}, index=rsi.index)

    def paired(frame, col, flip=None):
        tt = frame[frame.t].set_index("pid")[col]; tt = tt[~tt.index.duplicated()]
        nn = frame[~frame.t].set_index("pid")[col]; nn = nn[~nn.index.duplicated()]
        common = sorted(set(tt.index) & set(nn.index))
        d = (tt.loc[common] - nn.loc[common]).to_numpy(float)
        return d * flip if flip is not None else d

    d_rsi = paired(df, "rsi")
    d_d1 = paired(df, "d1")
    obs = ols_ci(d_rsi, np.column_stack([np.ones_like(d_d1), d_d1]),
                 ["intercept", "slope"])["intercept"]["beta"]

    null = []
    for _ in range(2000):
        f = rng.choice([-1.0, 1.0], size=len(d_rsi))
        null.append(ols_ci(d_rsi * f,
                           np.column_stack([np.ones_like(d_d1), d_d1 * f]),
                           ["intercept", "slope"])["intercept"]["beta"])
    null = np.array(null)
    p_perm = (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)
    check("within-pair label flips destroy the adjusted effect",
          abs(null.mean()) < 0.05 and p_perm < 0.001,
          f"null mean {null.mean():+.4f} sd {null.std():.3f}; "
          f"observed {obs:+.3f}, permutation P < {p_perm:.4f}")

    prev = json.load(open("results/GSE14520_differentiation.json"))
    check("reproduces the reported intercept exactly",
          abs(obs - prev["paired_adjusted_D1"]["intercept"]["beta"]) < 5e-4,
          f"{obs:+.4f} vs reported "
          f"{prev['paired_adjusted_D1']['intercept']['beta']:+.4f}")
else:
    print("  (GSE14520 cohort not built; check 3 skipped)")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
raise SystemExit(0 if ok else 1)
