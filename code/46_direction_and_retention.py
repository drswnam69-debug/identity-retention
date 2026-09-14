#!/usr/bin/env python3
"""6ad. Does anything about a signature predict its retention?
Reads only archived benchmark files; computes no new module score."""
import json
import statistics as st
from scipy.stats import spearmanr, mannwhitneyu


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
R = f"{IR_RESULTS}"
COH = [("GSE14520", "liver, array, 213 pairs"), ("TCGA_LIHC", "liver, RNA-seq, 50 pairs"),
       ("GSE76427", "liver, array, 52 pairs"), ("TCGA_LUAD", "lung, 58 pairs"),
       ("TCGA_KIRC", "kidney, 72 pairs")]
out = {}
for key, lab in COH:
    d = json.load(open(f"{R}/SIGNATURE_BENCHMARK_{key}.json"))
    sig = [s for s in d["signatures"] if s.get("status") == "ok"]
    if len(sig) < 5:
        continue
    ret = [s["retention_joint"] for s in sig]
    sh = [s["unadjusted_mean_delta"] for s in sig]
    sz = [s["n_genes"] for s in sig]
    r_dir, p_dir = spearmanr(sh, ret)
    r_sz, p_sz = spearmanr(sz, ret)
    rise = [s["retention_joint"] for s in sig if s["unadjusted_mean_delta"] > 0]
    fall = [s["retention_joint"] for s in sig if s["unadjusted_mean_delta"] < 0]
    rec = {"label": lab, "n": len(sig),
           "rho_retention_vs_signed_shift": round(float(r_dir), 4),
           "p_retention_vs_signed_shift": float(p_dir),
           "rho_retention_vs_set_size": round(float(r_sz), 4),
           "p_retention_vs_set_size": float(p_sz),
           "n_rise": len(rise), "n_fall": len(fall),
           "median_rise": round(st.median(rise), 4) if rise else None,
           "median_fall": round(st.median(fall), 4) if fall else None}
    if len(rise) > 2 and len(fall) > 2:
        rec["fold_range_within_rise"] = round(max(rise) / min(rise), 1)
        rec["fold_range_within_fall"] = round(max(fall) / min(fall), 1)
        u, pu = mannwhitneyu(rise, fall)
        rec["mannwhitney_p_rise_vs_fall"] = float(pu)
    out[key] = rec
    print(f"{key:10s} n={rec['n']:3d}  direction rho={rec['rho_retention_vs_signed_shift']:+.3f} "
          f"(P={rec['p_retention_vs_signed_shift']:.2e})  size rho={rec['rho_retention_vs_set_size']:+.3f} "
          f"(P={rec['p_retention_vs_set_size']:.2f})")

# name suffix, in the largest cohort only
d = json.load(open(f"{R}/SIGNATURE_BENCHMARK_GSE14520.json"))
sig = [s for s in d["signatures"] if s.get("status") == "ok"]
up = [s["retention_joint"] for s in sig if s["name"].endswith("_UP")]
dn = [s["retention_joint"] for s in sig if s["name"].endswith("_DN")]
u, p = mannwhitneyu(up, dn)
out["name_suffix_GSE14520"] = {"n_UP": len(up), "median_UP": round(st.median(up), 4),
                               "n_DN": len(dn), "median_DN": round(st.median(dn), 4),
                               "mannwhitney_p": float(p)}
print(f"name suffix: _UP n={len(up)} med={st.median(up):.3f} vs _DN n={len(dn)} "
      f"med={st.median(dn):.3f}  P={p:.4f}")
json.dump(out, open(f"{R}/DIRECTION_6ad.json", "w"), indent=1)
print("\nwrote", f"{R}/DIRECTION_6ad.json")
