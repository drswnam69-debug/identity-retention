#!/usr/bin/env python3
"""54_pool_shift_summary.py -- archive each 6ae pool's own paired shift.

53_negative_control_and_curvature.py wrote every draw's covariate shift to
NEGATIVE_CONTROL_6ae_draws.csv but summarized only the retention each pool
leaves. The manuscript states the unrestricted pool's median shift, so that
median needs a record of its own rather than living only in the prose. This
reads the deposited per-draw file and writes the summary beside it; it
recomputes nothing and cannot disagree with the run that produced the draws.
"""
import csv
import json
import os
import statistics
from collections import defaultdict

HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRAWS = f"{HOME}/results/NEGATIVE_CONTROL_6ae_draws.csv"
OUT = f"{HOME}/results/POOL_SHIFT_SUMMARY_6ae.json"

by = defaultdict(list)
with open(DRAWS, encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        by[row["pool"]].append(float(row["covariate_mean_paired_shift"]))

out = {"source": os.path.basename(DRAWS), "pools": {}}
for pool, vals in sorted(by.items()):
    vals.sort()
    out["pools"][pool] = {
        "n_draws": len(vals),
        "median_own_paired_shift": round(statistics.median(vals), 4),
        "min_own_paired_shift": round(vals[0], 4),
        "max_own_paired_shift": round(vals[-1], 4),
    }

with open(OUT, "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=2)
for pool, d in out["pools"].items():
    print(f"  {pool:14s} n={d['n_draws']:4d}  median shift "
          f"{d['median_own_paired_shift']:+.4f}")
print("written", OUT)
