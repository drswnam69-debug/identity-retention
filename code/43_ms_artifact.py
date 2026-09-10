#!/usr/bin/env python3
"""43_ms_artifact.py -- PREREGISTRATION 6ab.

Is the transcript-to-protein transfer failure an artifact of protein measurement?
Two predictors fixed in 6ab: the mean abundance of a signature's proteins in the
deposited matrix, and the proportion of its genes that appear in that matrix at
all. Neither is recomputed from anything; both are read off the same matrix and
the same signature list already used.
"""
from __future__ import annotations
import json, os, sys, importlib.util
import numpy as np, pandas as pd
from scipy.stats import spearmanr

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_cands = [_here, _os.path.join(_here, "rsi", "code"), _os.path.join(_here, "code"),
          _os.path.join(_os.path.dirname(_here), "code")]
if _os.environ.get("IR_ROOT"):
    _cands.insert(0, _os.path.join(_os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        if _c not in _sys.path:
            _sys.path.insert(0, _c)
        break
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, \
    DATA as IR_DATA, CODE as IR_CODE, FIGURES as IR_FIGURES, DOCS as IR_DOCS
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
def _l(n,f):
    s=importlib.util.spec_from_file_location(n,os.path.join(HERE,f))
    m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
import rsi_config as cfg

ALIAS={"MTARC1":"MARC1","MTARC2":"MARC2","G6PC1":"G6PC"}
RHO=-0.30

df=pd.read_csv(f"{IR_DOCS}/gao_proteins.tsv",sep="\t",index_col=0,low_memory=False)
df.index=[str(i).strip() for i in df.index]
df=df[[c for c in df.columns if str(c).strip()]].apply(pd.to_numeric,errors="coerce")
abund=df.mean(axis=1)
idx=set(df.index)
print("=== is the transfer failure a mass-spectrometry artifact? (PREREG 6ab) ===")
print(f"  matrix {df.shape[0]} proteins x {df.shape[1]} samples, "
      f"missing cells {int(df.isna().sum().sum())}")
print(f"  abundance: median {abund.median():.3f}, "
      f"range {abund.min():.3f} to {abund.max():.3f}")

G=json.load(open("results/GAO_MRNA_PROTEIN_6aa.json"))
sets=json.load(open("genesets/eligible.json"))
rows=[]
for r in G["signatures"]:
    genes=sets[r["name"]]
    found=[g if g in idx else ALIAS[g] for g in genes if g in idx or ALIAS.get(g) in idx]
    if not found: continue
    rows.append({"name":r["name"],"discrepancy":r["discrepancy"],
                 "n_genes":len(genes),"n_quantified":len(found),
                 "coverage":round(len(found)/len(genes),4),
                 "mean_abundance":round(float(abund.loc[found].mean()),4)})
R=pd.DataFrame(rows)
print(f"\n  {len(R)} signatures, the set already evaluable at both levels")
print(f"  coverage: median {R['coverage'].median():.3f} "
      f"(range {R['coverage'].min():.3f} to {R['coverage'].max():.3f})")

res={}
for pred,lab in (("mean_abundance","abundance"),("coverage","proteome coverage")):
    rho,p=spearmanr(R[pred],R["discrepancy"])
    ok = bool(rho<=RHO and p<0.05)
    res[pred]={"spearman_rho":round(float(rho),4),"spearman_p":float(p),
               "meets_prediction":ok}
    print(f"  [1] discrepancy vs {lab:<18} rho = {rho:+.3f}  (P = {p:.3g})"
          f"   prediction (rho <= -0.30) {'MET' if ok else 'not met'}")

print("\n  [2] where the study's own modules sit in the dynamic range")
own={}
n=len(abund)
for mod in ("REDUCTION","DRAIN"):
    print(f"    {mod}")
    for g in getattr(cfg,f"MODULE_{mod}"):
        k=g if g in idx else ALIAS.get(g)
        if k is None or k not in idx:
            print(f"      {g:<10} not quantified in this matrix"); own[g]=None; continue
        a=float(abund.loc[k]); pct=float((abund<a).mean()*100)
        own[g]={"abundance":round(a,4),"percentile":round(pct,1)}
        print(f"      {g:<10} abundance {a:+.3f}   {pct:.0f}th percentile")

verdict=("ARTIFACT_SUPPORTED" if any(v["meets_prediction"] for v in res.values())
         else "ARTIFACT_NOT_SUPPORTED")
print(f"\n  [3] 6ab verdict: {verdict}")

out={"plan":"PREREGISTRATION 6ab",
     "amendment_sha256_of_text_as_written":
        "c8f82f75e46b30d9746b800545c7eb3ac1bb44883b62eaeeab19add1287b6814",
     "n_proteins":int(df.shape[0]),"n_samples":int(df.shape[1]),
     "n_missing_cells":int(df.isna().sum().sum()),
     "n_signatures":int(len(R)),
     "coverage_median":round(float(R["coverage"].median()),4),
     "coverage_range":[round(float(R["coverage"].min()),4),
                       round(float(R["coverage"].max()),4)],
     "predictors":res,"own_modules":own,"threshold":RHO,"verdict":verdict,
     "signatures":R.to_dict("records")}
json.dump(out,open("results/MS_ARTIFACT_6ab.json","w"),indent=1)
print("\n  wrote results/MS_ARTIFACT_6ab.json")
