#!/usr/bin/env python3
"""16_gse14520_cohort.py -- assemble GSE14520 as a standard cohort directory.

PREREGISTRATION 6h (cohort 2) with the acquisition correction of 6i.

Input is the per-GSM SOFT transfer held in data/tmp14520/: one file per group
of five genes, each line "SYMBOL v1 v2 ... v445" in the sample order of
pheno.txt, values as served by GEO (RMA log2) and collapsed probeset->symbol
by max, the same rule io_utils uses.

This script writes results/GSE14520/{expr_log.tsv.gz, phenotype.tsv, rsi.tsv}
so that 12_differentiation_adjust.py can run on it unmodified, and it ASSERTS
that the panel values reproduce the unadjusted numbers already reported under
6h -- adding HNF1A must not move a single index value.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi as rsimod
import rsi_config as cfg

TMP = "data/tmp14520"
OUT = "results/GSE14520"
D1 = ["HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4",
      "ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG",
      "F2", "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1"]
D2_EXTRA = ["CYP2E1", "CYP3A4", "CYP1A2", "CYP2C9"]

ph = [l.split() for l in open(f"{TMP}/pheno.txt") if l.strip()]
gsms = [p[0] for p in ph]
assert len(gsms) == 445 and len(set(gsms)) == 445

rows = {}
for path in sorted(glob.glob(f"{TMP}/g*.txt")) + \
        ([f"{TMP}/hnf1a.txt"] if os.path.exists(f"{TMP}/hnf1a.txt") else []):
    for line in open(path):
        f = line.split()
        if not f:
            continue
        assert f[0] not in rows, f"duplicate gene {f[0]}"
        assert len(f) - 1 == 445, f"{f[0]}: {len(f)-1} values, expected 445"
        rows[f[0]] = [float(x) for x in f[1:]]

expr = pd.DataFrame(rows).T
expr.columns = gsms
expr.index.name = "gene"
print(f"=== GSE14520 cohort assembly (PREREG 6h/6i) ===")
print(f"  genes {expr.shape[0]}   samples {expr.shape[1]}")

# --- covariate completeness, reported not assumed --------------------------
for name, genes in (("D1", D1), ("D2", D1 + D2_EXTRA),
                    ("SUPPLY", cfg.MODULE_SUPPLY),
                    ("REDUCTION", cfg.MODULE_REDUCTION),
                    ("DRAIN", cfg.MODULE_DRAIN)):
    have = [g for g in genes if g in expr.index]
    miss = [g for g in genes if g not in expr.index]
    print(f"  {name:<10} {len(have)}/{len(genes)}" +
          (f"   absent from HG-U133A: {miss}" if miss else ""))

tissue = ["HCC tumor" if p[1] == "T" else "Adjacent non-tumor liver"
          for p in ph]
pheno = pd.DataFrame({"tissue": tissue,
                      "patient_id": [p[2] for p in ph],
                      "lcs_id": [p[3] for p in ph]}, index=gsms)
pheno.index.name = "sample"

rsi_tab, report = rsimod.compute_rsi(expr)
report["cohort"] = "GSE14520"
report["source"] = "GEO per-GSM SOFT, GPL3921, RMA log2 as served"
report["note_6i"] = ("HNF1A recovered from the retired symbol TCF1 "
                     "(210515_at, 216930_at; Entrez 6927)")

# --- assertion: the panel is unchanged by the 6i correction ----------------
prev = "results/gse14520_partial.json"
if os.path.exists(prev):
    old = json.load(open(prev))
    is_t = np.array([t == "HCC tumor" for t in tissue])
    r = rsi_tab["RSI"].to_numpy(float)
    dd = np.subtract.outer(r[is_t], r[~is_t]).ravel()
    hl = float(np.median(dd))
    assert abs(hl - old["H2_composite"]["hodges_lehmann"]) < 5e-4, \
        (f"the 6i correction moved the composite H2 "
         f"({hl:.4f} vs {old['H2_composite']['hodges_lehmann']:.4f}); "
         "adding HNF1A must not touch the index")
    for mod, key in (("REDUCTION", "REDUCTION"), ("DRAIN", "DRAIN")):
        col = {"REDUCTION": "z_reduction", "DRAIN": "z_drain"}[mod]
        s = rsi_tab[col].to_numpy(float)
        d = float(s[is_t].mean() - s[~is_t].mean())
        assert abs(d - old["modules_unadjusted"][key]["delta"]) < 5e-4, \
            f"{mod} delta moved: {d:.4f} vs {old['modules_unadjusted'][key]['delta']:.4f}"
    print("  [assert] panel values identical to the 6h unadjusted report "
          "(H2 and both module deltas match to 4 dp)")

os.makedirs(OUT, exist_ok=True)
expr.to_csv(f"{OUT}/expr_log.tsv.gz", sep="\t")
pheno.to_csv(f"{OUT}/phenotype.tsv", sep="\t")
rsi_tab.to_csv(f"{OUT}/rsi.tsv", sep="\t")
json.dump(report, open(f"{OUT}/rsi_report.json", "w"), indent=2, default=str)
n_pair = len(set(pheno.loc[[t == "HCC tumor" for t in tissue], "patient_id"]) &
             set(pheno.loc[[t != "HCC tumor" for t in tissue], "patient_id"]))
print(f"  tumor {sum(1 for t in tissue if t == 'HCC tumor')}   "
      f"adjacent {sum(1 for t in tissue if t != 'HCC tumor')}   "
      f"matched pairs {n_pair}")
print(f"  wrote -> {OUT}/")
