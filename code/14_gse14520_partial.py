#!/usr/bin/env python3
"""14_gse14520_partial.py -- GSE14520 under PREREGISTRATION 6h.

Acquisition was interrupted: the RSI panel genes are complete, the 22-gene D1
covariate is not (5 of 22). This script therefore reports ONLY what the
available genes support -- the UNADJUSTED module contrasts and the composite
RSI -- and refuses to compute the differentiation adjustment on a partial
covariate.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_lite as sl

D1_FULL = ["HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4", "ALB", "TTR", "TF",
           "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2", "CPS1",
           "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1"]
SUPPLY = ["HMGCR", "HMGCS1", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "PDSS1",
          "PDSS2", "COQ2", "COQ3", "COQ5", "COQ6", "COQ7", "COQ9"]
REDUCTION = ["CYB5R3", "CYB5R1", "AIFM2", "NQO1"]
DRAIN = ["MTARC1", "MTARC2", "POR"]

rows = []
i = 0
while os.path.exists(f"data/tmp14520/g{i}.txt"):
    for line in open(f"data/tmp14520/g{i}.txt"):
        f = line.split()
        rows.append((f[0], [float(x) for x in f[1:]]))
    i += 1
ph = [l.split() for l in open("data/tmp14520/pheno.txt")]
gsms = [p[0] for p in ph]
expr = pd.DataFrame({g: v for g, v in rows}).T
expr.columns = gsms
assert expr.shape[1] == 445

print("=== GSE14520 — PREREGISTRATION 6h, PARTIAL ACQUISITION ===")
print(f"  genes retrieved: {expr.shape[0]}   samples: {expr.shape[1]}")

d1_have = [g for g in D1_FULL if g in expr.index]
print(f"  D1 covariate: {len(d1_have)}/{len(D1_FULL)} genes present "
      f"({', '.join(d1_have)})")
if len(d1_have) < len(D1_FULL):
    print("  !! The differentiation adjustment is NOT run: the covariate "
          "registered in 6g is 22 genes and substituting a shorter one\n"
          "     would be exactly the post-hoc change 6g/6h forbid.")

is_t = np.array([p[1] == "T" for p in ph])
pid = [p[2] for p in ph]


def zmean(genes):
    g = [x for x in genes if x in expr.index]
    missing = sorted(set(genes) - set(g))
    z = expr.loc[g].sub(expr.loc[g].mean(axis=1), axis=0) \
                   .div(expr.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0).to_numpy(float), g, missing


res = {"plan": "PREREGISTRATION 6h", "cohort": "GSE14520",
       "status": "PARTIAL — differentiation covariate incomplete",
       "n": len(gsms), "n_tumor": int(is_t.sum()),
       "n_adjacent": int((~is_t).sum()),
       "d1_genes_available": d1_have,
       "d1_genes_missing": [g for g in D1_FULL if g not in expr.index],
       "adjustment_run": False}

print(f"\n  tumor {is_t.sum()}   adjacent {(~is_t).sum()}")
print("\n  UNADJUSTED module contrasts (what the acquired genes do support):")
mods = {}
for name, genes in (("SUPPLY", SUPPLY), ("REDUCTION", REDUCTION),
                    ("DRAIN", DRAIN)):
    s, g, miss = zmean(genes)
    d = s[is_t].mean() - s[~is_t].mean()
    _, _, p = sl.mannwhitney_u(s[is_t], s[~is_t])
    note = f"  (missing on HG-U133A: {miss})" if miss else ""
    print(f"    {name:<10} {len(g)}/{len(genes)} genes  "
          f"tumor {s[is_t].mean():+.3f} adjacent {s[~is_t].mean():+.3f}   "
          f"delta {d:+.3f}  P = {p:.3g}{note}")
    mods[name] = {"n_genes": len(g), "missing": miss,
                  "delta": round(float(d), 4), "p": p}
    globals()["s_" + name] = s
res["modules_unadjusted"] = mods

# composite RSI, exactly the locked definition
rsi = 0.5 * (s_SUPPLY + s_REDUCTION) - s_DRAIN
d = rsi[is_t].mean() - rsi[~is_t].mean()
_, _, p = sl.mannwhitney_u(rsi[is_t], rsi[~is_t])
dd = np.subtract.outer(rsi[is_t], rsi[~is_t]).ravel()
gt = sum((x > rsi[~is_t]).sum() for x in rsi[is_t])
lt = sum((x < rsi[~is_t]).sum() for x in rsi[is_t])
delta = (gt - lt) / (is_t.sum() * (~is_t).sum())
print(f"\n  H2 (composite RSI, tumor vs adjacent):")
print(f"    Hodges-Lehmann {np.median(dd):+.3f}   Cliff's delta {delta:+.3f}"
      f"   P = {p:.3g}")
res["H2_composite"] = {"hodges_lehmann": round(float(np.median(dd)), 4),
                       "cliffs_delta": round(float(delta), 4), "p": p}

# paired version
df = pd.DataFrame({"pid": pid, "t": is_t, "rsi": rsi,
                   "red": s_REDUCTION, "dra": s_DRAIN})
tt = df[df.t].set_index("pid"); tt = tt[~tt.index.duplicated()]
nn = df[~df.t].set_index("pid"); nn = nn[~nn.index.duplicated()]
common = sorted(set(tt.index) & set(nn.index))
print(f"\n  paired, {len(common)} patient-matched pairs:")
res["n_pairs"] = len(common)
for key, lab in (("rsi", "RSI"), ("red", "REDUCTION"), ("dra", "DRAIN")):
    dv = (tt.loc[common, key] - nn.loc[common, key]).to_numpy(float)
    n = len(dv)
    t = dv.mean() / (dv.std(ddof=1) / np.sqrt(n))
    pv = sl._t_two_sided(float(t), n - 1)
    print(f"    {lab:<10} mean paired difference {dv.mean():+.3f}   "
          f"P = {pv:.3g}")
    res[f"paired_{key}"] = {"mean_diff": round(float(dv.mean()), 4), "p": pv}

with open("results/gse14520_partial.json", "w") as fh:
    json.dump(res, fh, indent=2, default=float)
print("\n  wrote -> results/gse14520_partial.json")
