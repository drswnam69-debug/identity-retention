#!/usr/bin/env python3
"""38_tcga_premise.py -- PREREGISTRATION 6w premise check.

Reported BEFORE any contrast, under the standing rule that a public matrix must
be shown to contain what the analysis assumes. If the hepatocyte identity score
is not lower in tumor, the cohort is reported as not testable and nothing else
from it is reported.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


_da = _load("d", "12_differentiation_adjust.py")
_ca = _load("c", "18_composition_adjust.py")
_bm = _load("b", "30_signature_benchmark.py")
D1, D2_EXTRA, C1, zmean = _da.D1, _da.D2_EXTRA, _ca.C1, _bm.zmean
import rsi_config as cfg  # noqa: E402

expr = pd.read_csv("data/TCGA_LIHC_symbols.tsv.gz", sep="\t", index_col=0)
ph = pd.read_csv("results/TCGA_LIHC/phenotype.tsv", sep="\t", index_col=0)
expr = expr[ph.index]
prep = json.load(open("results/TCGA_LIHC/PREP_6w.json"))

print("=== TCGA-LIHC premise check (PREREG 6w) ===")
print("\n  [1] what this matrix is, before any contrast")
print(f"      source            {prep['source_file']}")
print(f"      declared unit     {prep['declared_unit']}, range "
      f"{prep['value_range'][0]} to {prep['value_range'][1]}")
print( "      processing already applied by the source: STAR quantification, "
       "aliquot averaging,")
print( "        then log2(x + 1). No further normalization is added here.")
print(f"      symbols           {prep['n_symbols']}  "
      f"({prep['n_duplicate_rows_collapsed']} duplicate rows collapsed by max)")
print(f"      samples kept      {prep['n_tumor_01']} tumor (01), "
      f"{prep['n_adjacent_11']} adjacent (11)")
print(f"      discarded         {prep['sample_types_discarded']} "
      f"(02 = recurrent tumor, excluded by the 6w rule)")
print(f"      PAIRS             {prep['n_pairs']}")

is_t = (ph["tissue"] == "tumor").to_numpy()
fr = pd.DataFrame({"pid": ph["patient_id"].astype(str).values, "t": is_t},
                  index=ph.index)


def paired(s):
    f = fr.assign(v=s.to_numpy())
    tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
    nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
    k = sorted(set(tt.index) & set(nn.index))
    return (tt.loc[k] - nn.loc[k]).to_numpy(float), k


print("\n  [2] coverage of the locked gene lists")
cov = {}
for name, genes in [("D1", D1), ("D2", D1 + D2_EXTRA), ("C1", C1),
                    ("SUPPLY", cfg.MODULE_SUPPLY),
                    ("REDUCTION", cfg.MODULE_REDUCTION),
                    ("DRAIN", cfg.MODULE_DRAIN)]:
    found = [g for g in genes if g in expr.index]
    cov[name] = {"n_found": len(found), "n_defined": len(genes),
                 "missing": [g for g in genes if g not in found]}
    print(f"      {name:<10} {len(found)}/{len(genes)}"
          + (f"   missing {cov[name]['missing']}" if cov[name]["missing"] else ""))

d1s, _ = zmean(expr, D1)
c1s, _ = zmean(expr, C1)
d_d1, pairs = paired(d1s)
d_c1, _ = paired(c1s)
p_d1 = float(wilcoxon(d_d1, alternative="two-sided").pvalue)
p_c1 = float(wilcoxon(d_c1, alternative="two-sided").pvalue)

print("\n  [3] THE PREMISE: does hepatocyte identity fall in tumor?")
print(f"      D1 paired difference, tumor minus adjacent: {d_d1.mean():+.4f}"
      f"   P = {p_d1:.3g}   ({len(pairs)} pairs)")
premise_ok = bool(d_d1.mean() < 0 and p_d1 < 0.05)
print(f"      direction: {'LOWER in tumor' if d_d1.mean() < 0 else 'HIGHER in tumor'}")
print(f"      >>> PREMISE {'HOLDS' if premise_ok else 'FAILS'} <<<")
if not premise_ok:
    print("      The confound this measure exists to remove is absent. Under 6w "
          "this cohort is\n      reported as NOT TESTABLE and no retention value "
          "from it is reported.")

print("\n  [4] the individual identity genes, for inspection")
rows = []
for g in D1:
    if g not in expr.index:
        continue
    z = expr.loc[g]
    z = (z - z.mean()) / z.std()
    dv, _ = paired(z)
    pv = float(wilcoxon(dv, alternative="two-sided").pvalue)
    rows.append((g, float(dv.mean()), pv))
rows.sort(key=lambda r: r[1])
for g, m, pv in rows:
    print(f"      {g:<10} {m:+.3f}   P = {pv:.2e}")
n_down = sum(1 for _, m, _ in rows if m < 0)
print(f"      {n_down} of {len(rows)} identity genes fall in tumor")

print(f"\n  [5] composition covariate C1: {d_c1.mean():+.4f}  P = {p_c1:.3g}")

res = {"plan": "PREREGISTRATION 6w premise check",
       "amendment_sha256_of_text_as_written": prep["amendment_sha256_of_text_as_written"],
       "cohort": "TCGA_LIHC", "n_pairs": len(pairs),
       "preparation": prep, "coverage": cov,
       "D1_paired_delta": round(float(d_d1.mean()), 4), "D1_p": p_d1,
       "D1_direction": "LOWER in tumor" if d_d1.mean() < 0 else "HIGHER in tumor",
       "identity_confound_present": premise_ok,
       "n_identity_genes_falling": n_down, "n_identity_genes": len(rows),
       "per_gene_D1": [{"gene": g, "paired_delta": round(m, 4), "p": pv}
                       for g, m, pv in rows],
       "C1_paired_delta": round(float(d_c1.mean()), 4), "C1_p": p_c1}
with open("results/TCGA_LIHC/PREMISE_6w.json", "w") as fh:
    json.dump(res, fh, indent=1)
print("\n  wrote results/TCGA_LIHC/PREMISE_6w.json")
