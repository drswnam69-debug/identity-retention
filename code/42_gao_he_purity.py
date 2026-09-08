#!/usr/bin/env python3
"""42_gao_he_purity.py -- PREREGISTRATION 6aa, Part B.

Does a pathologist's tumor purity, read from hematoxylin and eosin staining,
change the identity-retention fractions? If it does not, then the strongest form
of the claim in Section 4.1 holds: the identity adjustment is not doing the job a
purity measure does, and this is shown against purity as the field measures it
rather than against another gene score.
"""
from __future__ import annotations
import importlib.util, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
_da = _load("d", "12_differentiation_adjust.py"); _ca = _load("c", "18_composition_adjust.py")
ols_ci, D1, C1 = _da.ols_ci, _da.D1, _ca.C1
import rsi_config as cfg
ALIAS = {"MTARC1": "MARC1", "MTARC2": "MARC2", "G6PC1": "G6PC"}
MOVE_THRESHOLD = 0.10                      # 6aa decision rule

df = pd.read_csv("/home/claude/gao_proteins.tsv", sep="\t", index_col=0, low_memory=False)
df.index = [str(i).strip() for i in df.index]
df = df[[c for c in df.columns if str(c).strip()]].apply(pd.to_numeric, errors="coerce")
cl = pd.read_csv("data/gao_clinical.tsv", sep="\t")
cl["he_purity"] = pd.to_numeric(cl["he_purity"], errors="coerce")

tum = [c for c in df.columns if str(c).startswith("T")]
non = [c for c in df.columns if str(c).startswith("N")]
pt = {c[1:]: c for c in tum}; pn = {c[1:]: c for c in non}
pids_all = sorted(set(pt) & set(pn))
pur = {}
for _, r in cl.iterrows():
    t = str(r["tumor_id"]).strip()
    if t.startswith("T") and pd.notna(r["he_purity"]):
        pur[t[1:]] = float(r["he_purity"])
pids = [p for p in pids_all if p in pur]
print("=== Gao histological tumor purity (PREREG 6aa part B) ===")
print(f"  pairs in the protein matrix: {len(pids_all)}")
print(f"  of those with a recorded HE purity: {len(pids)}")
pv = np.array([pur[p] for p in pids])
print(f"  purity: median {np.median(pv):.2f}, range {pv.min():.2f} to {pv.max():.2f}")

def resolve(genes):
    idx = set(df.index)
    return [g if g in idx else ALIAS[g] for g in genes
            if g in idx or ALIAS.get(g) in idx]
def zm(genes):
    g = resolve(genes)
    sub = df.loc[g]; sd = sub.std(axis=1); sub = sub[(sd > 0).values]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0)
    return z.mean(axis=0), sub.shape[0]
d1s, _ = zm(D1); c1s, _ = zm(C1)

print("\n  [B1/B2] does purity track the two covariates, in tumors?")
t_d1 = np.array([d1s[pt[p]] for p in pids]); t_c1 = np.array([c1s[pt[p]] for p in pids])
r_c1, p_c1 = spearmanr(pv, t_c1); r_d1, p_d1s = spearmanr(pv, t_d1)
print(f"      HE purity vs C1 (composition): rho = {r_c1:+.3f}  (P = {p_c1:.3g})")
print(f"      HE purity vs D1 (identity)   : rho = {r_d1:+.3f}  (P = {p_d1s:.3g})")

# paired differences on the patients with purity
def paired(s):
    return np.array([s[pt[p]] - s[pn[p]] for p in pids], float)
d_d1 = paired(d1s); d_c1 = paired(c1s)
# the adjacent liver is taken as fully non-tumor, so the paired change in tumor
# content is the tumor's purity itself; this is stated in the manuscript
d_pur = pv - 0.0
n = len(pids)
X_base = np.column_stack([np.ones(n), d_d1, d_c1])
X_pur = np.column_stack([np.ones(n), d_d1, d_c1, (d_pur - d_pur.mean())])

print(f"\n  [B3] retention with and without histological purity, {n} pairs")
out = {}
for mod in ("REDUCTION", "DRAIN"):
    dv = paired(zm(getattr(cfg, f"MODULE_{mod}"))[0])
    un = float(dv.mean())
    f0 = ols_ci(dv, X_base, ["intercept", "dD1", "dC1"])
    f1 = ols_ci(dv, X_pur, ["intercept", "dD1", "dC1", "purity"])
    r0 = abs(f0["intercept"]["beta"]) / abs(un)
    r1 = abs(f1["intercept"]["beta"]) / abs(un)
    out[mod] = {"unadjusted_shift": round(un, 4),
                "retention_without_purity": round(r0, 4),
                "retention_with_purity": round(r1, 4),
                "move": round(abs(r1 - r0), 4),
                "purity_coefficient": f1["purity"]}
    print(f"      {mod:<10} {r0:.4f} -> {r1:.4f}   moves {abs(r1-r0):.4f}"
          f"   (purity term P = {f1['purity']['p']:.3g})")
moves = [out[m]["move"] for m in out]
verdict = ("PURITY_ADDS_NOTHING" if max(moves) <= MOVE_THRESHOLD
           else "PARTLY_CONFOUNDED_WITH_PURITY")
print(f"\n  [B4] 6aa part B verdict: {verdict}  "
      f"(largest move {max(moves):.4f}, threshold {MOVE_THRESHOLD})")

res = {"plan": "PREREGISTRATION 6aa part B",
       "amendment_sha256_of_text_as_written":
           "2cbe2a9ab7928f33899312e5d52f72f7142859b4c0e453e5bbb0a53b2eebe42e",
       "n_pairs_total": len(pids_all), "n_pairs_with_purity": n,
       "purity_median": round(float(np.median(pv)), 4),
       "purity_range": [round(float(pv.min()), 4), round(float(pv.max()), 4)],
       "purity_vs_C1": {"rho": round(float(r_c1), 4), "p": float(p_c1)},
       "purity_vs_D1": {"rho": round(float(r_d1), 4), "p": float(p_d1s)},
       "modules": out, "move_threshold": MOVE_THRESHOLD, "verdict": verdict}
json.dump(res, open("results/GAO_HE_PURITY_6aa.json", "w"), indent=1)
print("\n  wrote results/GAO_HE_PURITY_6aa.json")
