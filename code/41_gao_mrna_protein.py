#!/usr/bin/env python3
"""41_gao_mrna_protein.py -- PREREGISTRATION 6aa, Part A.

Does measured mRNA to protein concordance explain the transfer failure?

The correlations are the source study's own, published per gene across the same
159 paired cases. Nothing is recomputed and nothing is filtered.
"""
from __future__ import annotations
import importlib.util, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
_da = _load("d", "12_differentiation_adjust.py"); _ca = _load("c", "18_composition_adjust.py")
_bm = _load("b", "30_signature_benchmark.py")
ols_ci, D1, C1, zmean = _da.ols_ci, _da.D1, _ca.C1, _bm.zmean
import rsi_config as cfg
MIN_SET, MAX_SET, NULL_P = _bm.MIN_SET, _bm.MAX_SET, _bm.NULL_P
MIN_ON = 10
ALIAS = {"MTARC1": "MARC1", "MTARC2": "MARC2", "G6PC1": "G6PC"}
RHO_THRESHOLD = -0.30

corr = pd.read_csv("data/gao_mrna_protein_corr.tsv", sep="\t")
corr["spearman"] = pd.to_numeric(corr["spearman"], errors="coerce")
CO = dict(zip(corr["gene"], corr["spearman"]))
allc = corr["spearman"].dropna().to_numpy()
print("=== Gao mRNA to protein concordance (PREREG 6aa part A) ===")
print(f"  {len(allc)} genes; median rho {np.median(allc):.3f} "
      f"(IQR {np.percentile(allc,25):.3f} to {np.percentile(allc,75):.3f}), "
      f"range {allc.min():.3f} to {allc.max():.3f}")

def look(g):
    return CO.get(g, CO.get(ALIAS.get(g), None))

print("\n  [A1] the study's own genes")
for name, genes in (("REDUCTION", cfg.MODULE_REDUCTION), ("DRAIN", cfg.MODULE_DRAIN)):
    print(f"    {name}")
    for g in genes:
        v = look(g)
        pct = (float((allc < v).mean()) * 100) if v is not None else None
        print(f"      {g:<10} " + (f"rho = {v:+.3f}   {pct:.0f}th percentile"
                                   if v is not None else "not in the table"))
d1v = [look(g) for g in D1]; d1v = [v for v in d1v if v is not None]
print(f"    D1: {len(d1v)}/{len(D1)} genes, median rho {np.median(d1v):+.3f}")

# --- A2: protein-level benchmark, then the association ----------------------
df = pd.read_csv("/home/claude/gao_proteins.tsv", sep="\t", index_col=0, low_memory=False)
df.index = [str(i).strip() for i in df.index]
df = df[[c for c in df.columns if str(c).strip()]].apply(pd.to_numeric, errors="coerce")
tum = [c for c in df.columns if str(c).startswith("T")]
non = [c for c in df.columns if str(c).startswith("N")]
pt = {c[1:]: c for c in tum}; pn = {c[1:]: c for c in non}
pids = sorted(set(pt) & set(pn)); n = len(pids)
print(f"\n  [A2] Gao protein matrix: {df.shape[0]} proteins, {n} pairs")

def resolve(genes):
    idx = set(df.index)
    return [g if g in idx else ALIAS[g] for g in genes
            if g in idx or ALIAS.get(g) in idx]
def zm(genes):
    g = resolve(genes)
    if not g: return None, 0
    sub = df.loc[g]; sd = sub.std(axis=1); sub = sub[(sd > 0).values]
    if sub.empty: return None, 0
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), sub.shape[0]
def paired(s):
    return np.array([s[pt[p]] - s[pn[p]] for p in pids], float)

d_d1 = paired(zm(D1)[0]); d_c1 = paired(zm(C1)[0])
X = np.column_stack([np.ones(n), d_d1, d_c1])
sets = {k: v for k, v in json.load(open("genesets/eligible.json")).items()
        if MIN_SET <= len(v) <= MAX_SET}
arch = {r["name"]: r for r in json.load(
    open("results/SIGNATURE_BENCHMARK_GSE14520.json"))["signatures"]
    if r["status"] == "ok"}
rows = []
for name, genes in sorted(sets.items()):
    s, n_on = zm(genes)
    if s is None or n_on < MIN_ON: continue
    dv = paired(s); un = float(dv.mean())
    try: wp = float(wilcoxon(dv).pvalue)
    except Exception: continue
    if not np.isfinite(wp) or wp >= NULL_P: continue
    f = ols_ci(dv, X, ["intercept", "dD1", "dC1"])
    prot = abs(f["intercept"]["beta"]) / abs(un)
    cs = [look(g) for g in genes]; cs = [c for c in cs if c is not None]
    if len(cs) < 5 or name not in arch: continue
    rows.append({"name": name, "n_proteins": n_on,
                 "retention_protein": round(prot, 4),
                 "retention_transcript": arch[name]["retention_joint"],
                 "discrepancy": round(abs(prot - arch[name]["retention_joint"]), 4),
                 "n_genes_with_corr": len(cs),
                 "mean_mrna_protein_rho": round(float(np.mean(cs)), 4)})
R = pd.DataFrame(rows)
print(f"      {len(R)} signatures evaluable at BOTH levels with correlations")
rho, p = spearmanr(R["mean_mrna_protein_rho"], R["discrepancy"])
print(f"\n  [A2] Spearman rho between mean mRNA-protein concordance and "
      f"|transcript − protein| retention discrepancy")
print(f"      rho = {rho:+.3f}  (P = {p:.3g})   threshold fixed in 6aa: rho <= -0.30")
verdict = ("SUPPORTED" if (rho <= RHO_THRESHOLD and p < 0.05)
           else ("POSITIVE" if (rho >= 0.30 and p < 0.05) else "NOT_SUPPORTED"))
print(f"  [A3] 6aa part A verdict: {verdict}")
print(f"      median discrepancy {R['discrepancy'].median():.3f}; "
      f"median protein retention {R['retention_protein'].median():.3f} vs "
      f"transcript {R['retention_transcript'].median():.3f}")

res = {"plan": "PREREGISTRATION 6aa part A",
       "amendment_sha256_of_text_as_written":
           "2cbe2a9ab7928f33899312e5d52f72f7142859b4c0e453e5bbb0a53b2eebe42e",
       "n_genes_with_correlation": int(len(allc)),
       "all_genes_median_rho": round(float(np.median(allc)), 4),
       "all_genes_iqr": [round(float(np.percentile(allc, 25)), 4),
                         round(float(np.percentile(allc, 75)), 4)],
       "own_genes": {g: (round(look(g), 4) if look(g) is not None else None)
                     for g in list(cfg.MODULE_REDUCTION) + list(cfg.MODULE_DRAIN)},
       "D1_median_rho": round(float(np.median(d1v)), 4), "D1_n_with_corr": len(d1v),
       "n_pairs": n, "n_signatures_both_levels": int(len(R)),
       "spearman_rho": round(float(rho), 4), "spearman_p": float(p),
       "threshold": RHO_THRESHOLD, "verdict": verdict,
       "median_discrepancy": round(float(R["discrepancy"].median()), 4),
       "median_retention_protein": round(float(R["retention_protein"].median()), 4),
       "median_retention_transcript": round(float(R["retention_transcript"].median()), 4),
       "signatures": R.to_dict("records")}
json.dump(res, open("results/GAO_MRNA_PROTEIN_6aa.json", "w"), indent=1)
print("\n  wrote results/GAO_MRNA_PROTEIN_6aa.json")
