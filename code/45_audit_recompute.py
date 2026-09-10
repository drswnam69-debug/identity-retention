#!/usr/bin/env python3
"""Recompute two manuscript claims that no archived result file carries:
   (A) per-protein identity retention for the three N-reductive proteins in Gao;
   (B) the restricted-D1 sensitivity in GSE14520 for BOTH modules.
Uses the same scoring and adjustment as 32_protein_validation.py and
12_differentiation_adjust.py, so the values are comparable to the archived ones.
"""
import importlib.util, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

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
HERE = f"{IR_ROOT}/code"
def imp(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
sys.path.insert(0, HERE)
import rsi_config as cfg

sl = imp("stats_lite", "stats_lite.py")

D1 = ["ALB","TTR","TF","SERPINA1","AHSG","APOH","FGA","FGB","FGG","F2","CPS1",
      "OTC","ARG1","TAT","G6PC1","PCK1","ASGR1","HNF4A","HNF1A","FOXA1","FOXA2","NR1H4"]
C1 = ["PTPRC","CD53","LAPTM5","CD3E","CD2","CD68","AIF1","ITGAM","LCP1",
      "PECAM1","VWF","CDH5","ENG","CLEC4G","COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRB","ACTA2"]
ALIAS = {"MTARC1":"MARC1","MTARC2":"MARC2","G6PC1":"G6PC"}

def resolve(genes, index):
    out=[]
    for g in genes:
        if g in index: out.append(g)
        elif ALIAS.get(g) in index: out.append(ALIAS[g])
    return out

def zmean(mat, genes):
    g = resolve(genes, set(mat.index))
    if not g: return None,0
    sub = mat.loc[g]; sd = sub.std(axis=1); sub = sub[(sd>0).values]
    if sub.empty: return None,0
    z = sub.sub(sub.mean(axis=1),axis=0).div(sub.std(axis=1),axis=0)
    return z.mean(axis=0), sub.shape[0]

def ols(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n,p = X.shape
    s2 = resid @ resid/(n-p)
    se = np.sqrt(np.diag(s2*np.linalg.inv(X.T@X)))
    t = beta/se
    pv = [float(sl._t_two_sided(float(tt), n-p)) for tt in t]
    return beta, se, t, pv

# ---------------------------------------------------------------- A. Gao
print("=== A. per-protein retention, Gao 2019 (D1, 19 proteins) ===")
df = pd.read_csv(f"{IR_DOCS}/gao_proteins.tsv", sep="\t", index_col=0, low_memory=False)
df.index=[str(i).strip() for i in df.index]
df = df[[c for c in df.columns if str(c).strip()]].apply(pd.to_numeric, errors="coerce")
tum=[c for c in df.columns if str(c).startswith("T")]; non=[c for c in df.columns if str(c).startswith("N")]
pt={c[1:]:c for c in tum}; pn={c[1:]:c for c in non}
pids=sorted(set(pt)&set(pn))
if float(np.nanmax(df.values))>100: df=np.log2(df.clip(lower=0)+1.0)
sD1,_=zmean(df,D1)
paired=lambda s: np.array([s[pt[p]]-s[pn[p]] for p in pids],float)
dD1=paired(sD1); n=len(pids)
gao={}
for g in ["MTARC1","MTARC2","POR"]:
    gg = ALIAS.get(g,g) if ALIAS.get(g,g) in df.index else g
    row=df.loc[gg]; z=(row-row.mean())/row.std()
    dv=paired(z); un=float(dv.mean())
    b,se,t,pv=ols(dv, np.column_stack([np.ones(n),dD1]))
    ret=abs(b[0])/abs(un)
    gao[g]={"unadjusted":round(un,4),"wilcoxon_p":float(wilcoxon(dv).pvalue),
            "intercept":round(float(b[0]),4),"intercept_p":pv[0],"retention":round(float(ret),4)}
    print(f"  {g:8s} unadjusted {un:+.3f}  intercept {b[0]:+.3f}  retention {ret*100:.1f}%")

# ---------------------------------------------------------------- B. GSE14520
print("\n=== B. GSE14520, D1 restricted to the 19 genes measured as protein ===")
expr = pd.read_csv(f"{IR_DATA}/GSE14520_symbols.tsv.gz", sep="\t", index_col=0)
ph = pd.read_csv(f"{IR_RESULTS}/GSE14520/phenotype.tsv", sep="\t")
print("  phenotype columns:", list(ph.columns)[:12])
json.dump({"gao_per_protein_retention_D1":gao}, open(f"{IR_RESULTS}/AUDIT_RECOMPUTE_2026-09-09.json","w"), indent=1)
