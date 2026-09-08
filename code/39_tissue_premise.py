#!/usr/bin/env python3
"""39_tissue_premise.py -- premise check for a new tissue (PREREG 6x, 6y).

Reported BEFORE any contrast. The tissue identity covariate was written out in full
in its amendment before the matrix was opened, is imported from its config module,
and is not modified here. If K1 is not lower in tumor, the kidney cohort is reported
as not testable and nothing else from it is reported.
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
D1, C1, zmean = _da.D1, _ca.C1, _bm.zmean
import argparse  # noqa: E402
_ap = argparse.ArgumentParser()
_ap.add_argument("--cohort", default="TCGA_KIRC")
_ap.add_argument("--identity", choices=["K1", "L1"], default="K1")
_ap.add_argument("--prep", default=None)
_a = _ap.parse_args()
if _a.identity == "K1":
    from k1_config import K1 as IDENT, K1_EFFECTORS as EFF, K1_REGULATORS as REG
    IDENT_NAME, PLAN = "K1", "6x"
else:
    from l1_config import L1 as IDENT, L1_EFFECTORS as EFF, L1_REGULATORS as REG
    IDENT_NAME, PLAN = "L1", "6y"
COH = _a.cohort
PREP = _a.prep or f"results/{COH}/PREP_{PLAN}.json"

expr = pd.read_csv(f"data/{COH}_symbols.tsv.gz", sep="\t", index_col=0)
ph = pd.read_csv(f"results/{COH}/phenotype.tsv", sep="\t", index_col=0)
expr = expr[ph.index]
prep = json.load(open(PREP))

print(f"=== {COH} premise check (PREREG {PLAN}), identity covariate {IDENT_NAME} ===")
print("\n  [1] what this matrix is, before any contrast")
print(f"      source          {prep['source_file']}")
print(f"      declared unit   {prep['declared_unit']}, range "
      f"{prep['value_range'][0]} to {prep['value_range'][1]}, missing "
      f"{prep['n_missing_values']}")
print( "      processing already applied: STAR quantification, aliquot averaging,")
print( "        then log2(x + 1). Nothing further is added here.")
print(f"      samples kept    {prep['n_tumor_01']} tumor (01), "
      f"{prep['n_adjacent_11']} adjacent (11); discarded "
      f"{prep['sample_types_discarded']}")
print(f"      PAIRS           {prep['n_pairs']}")

is_t = (ph["tissue"] == "tumor").to_numpy()
fr = pd.DataFrame({"pid": ph["patient_id"].astype(str).values, "t": is_t},
                  index=ph.index)


def paired(s):
    f = fr.assign(v=s.to_numpy())
    tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
    nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
    k = sorted(set(tt.index) & set(nn.index))
    return (tt.loc[k] - nn.loc[k]).to_numpy(float), k


print("\n  [2] coverage of the locked covariates")
cov = {}
for name, genes in [(IDENT_NAME, IDENT), ("C1", C1), ("D1 (liver, not used)", D1)]:
    found = [g for g in genes if g in expr.index]
    cov[name] = {"n_found": len(found), "n_defined": len(genes),
                 "missing": [g for g in genes if g not in found]}
    print(f"      {name:<22} {len(found)}/{len(genes)}"
          + (f"   missing {cov[name]['missing']}" if cov[name]["missing"] else ""))

k1s, _ = zmean(expr, IDENT)
c1s, _ = zmean(expr, C1)
d_k1, pairs = paired(k1s)
d_c1, _ = paired(c1s)
p_k1 = float(wilcoxon(d_k1, alternative="two-sided").pvalue)
p_c1 = float(wilcoxon(d_c1, alternative="two-sided").pvalue)

print(f"\n  [3] THE PREMISE: does {IDENT_NAME} identity fall in tumor?")
print(f"      {IDENT_NAME} paired difference, tumor minus adjacent: {d_k1.mean():+.4f}"
      f"   P = {p_k1:.3g}   ({len(pairs)} pairs)")
premise_ok = bool(d_k1.mean() < 0 and p_k1 < 0.05)
print(f"      direction: {'LOWER in tumor' if d_k1.mean() < 0 else 'HIGHER in tumor'}")
print(f"      >>> PREMISE {'HOLDS' if premise_ok else 'FAILS'} <<<")
if not premise_ok:
    print("      Under 6x the kidney cohort is reported as NOT TESTABLE and the "
          "manuscript states\n      that the covariate construction did not "
          "transfer to a second tissue.")

print(f"\n  [4] the individual {IDENT_NAME} genes")
rows = []
for g in IDENT:
    if g not in expr.index:
        continue
    z = expr.loc[g]
    z = (z - z.mean()) / z.std()
    dv, _ = paired(z)
    rows.append((g, "regulator" if g in REG else "effector",
                 float(dv.mean()), float(wilcoxon(dv).pvalue)))
rows.sort(key=lambda r: r[2])
for g, cls, m, pv in rows:
    print(f"      {g:<10} {cls:<10} {m:+.3f}   P = {pv:.2e}")
eff = [r[2] for r in rows if r[1] == "effector"]
reg = [r[2] for r in rows if r[1] == "regulator"]
print(f"      effectors  mean {np.mean(eff):+.3f}  ({sum(1 for x in eff if x<0)}/{len(eff)} fall)")
print(f"      regulators mean {np.mean(reg):+.3f}  ({sum(1 for x in reg if x<0)}/{len(reg)} fall)")

print(f"\n  [5] composition covariate C1: {d_c1.mean():+.4f}  P = {p_c1:.3g}")
if PLAN == "6x":
    print("      6x predicted in advance that C1 would RISE in tumor here, the "
          "opposite of liver:")
    print(f"      prediction {'CONFIRMED' if d_c1.mean() > 0 and p_c1 < 0.05 else 'NOT confirmed'}")
else:
    print("      6y registered no directional prediction for C1; reported as observed.")

res = {"plan": f"PREREGISTRATION {PLAN} premise check",
       "amendment_sha256_of_text_as_written":
           ("aba7b918a4c9842e5ea027fae14d61c6e279f5c4df2da235e9fbbfe3e732259c"
            if PLAN == "6x" else
            "5fd0c9b26945dc0d134cf5aceda83096b406314dc49672824a80be2cbcc64340"),
       "cohort": COH, "identity_covariate": IDENT_NAME, "n_pairs": len(pairs), "preparation": prep,
       "identity_genes": IDENT, "coverage": cov,
       "identity_paired_delta": round(float(d_k1.mean()), 4), "identity_p": p_k1,
       "identity_direction": "LOWER in tumor" if d_k1.mean() < 0 else "HIGHER in tumor",
       "identity_confound_present": premise_ok,
       "per_gene_identity": [{"gene": g, "class": c, "paired_delta": round(m, 4), "p": pv}
                       for g, c, m, pv in rows],
       "effector_mean": round(float(np.mean(eff)), 4),
       "regulator_mean": round(float(np.mean(reg)), 4),
       "n_effectors_falling": int(sum(1 for x in eff if x < 0)),
       "n_regulators_falling": int(sum(1 for x in reg if x < 0)),
       "C1_paired_delta": round(float(d_c1.mean()), 4), "C1_p": p_c1,
       "C1_rise_predicted": PLAN == "6x",
       "C1_prediction_confirmed": bool(d_c1.mean() > 0 and p_c1 < 0.05)}
with open(f"results/{COH}/PREMISE_{PLAN}.json", "w") as fh:
    json.dump(res, fh, indent=1)
print(f"\n  wrote results/{COH}/PREMISE_{PLAN}.json")
