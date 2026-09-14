#!/usr/bin/env python3
"""28_diagnostics_6o.py -- PREREGISTRATION 6o.

6o-A: is DRAIN's collapse under D2 attributable to POR, as the manuscript says?
      Refit D2 on DRAIN' = mean z(MTARC1, MTARC2), which contains no POR.
6o-B: does removing AIFM2 raise or lower REDUCTION's retained fraction?
      Test it in GSE76427, where AIFM2 is measured.

Decision rules were fixed in 6o before this ran.
"""
import importlib.util, json, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rsi_config as cfg

def _load(n, f):
    sp = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(sp); sp.loader.exec_module(m); return m

_da = _load("d", "12_differentiation_adjust.py")
_ca = _load("c", "18_composition_adjust.py")
ols_ci, D1, D2_EXTRA, C1 = _da.ols_ci, _da.D1, _da.D2_EXTRA, _ca.C1

SRC = {"GSE76427": ("data/GSE76427_d1_panel.tsv", "data/GSE76427_composition_panel.tsv"),
       "GSE14520": (None, "data/GSE14520_composition_panel.tsv")}

def zmean(mat, genes):
    g = [x for x in genes if x in mat.index]
    z = mat.loc[g].sub(mat.loc[g].mean(axis=1), axis=0).div(mat.loc[g].std(axis=1), axis=0)
    return z.mean(axis=0), g

def frame(cohort):
    e = pd.read_csv(f"results/{cohort}/expr_log.tsv.gz", sep="\t", index_col=0)
    ph = pd.read_csv(f"results/{cohort}/phenotype.tsv", sep="\t", index_col=0)
    rsi = pd.read_csv(f"results/{cohort}/rsi.tsv", sep="\t", index_col=0)
    keep = [c for c in e.columns if c in ph.index and c in rsi.index]
    e, ph = e[keep], ph.loc[keep]
    d1s, c1s = SRC[cohort]
    d1src = e if d1s is None else pd.read_csv(d1s, sep="\t", index_col=0)[keep]
    is_t = ~ph["tissue"].astype(str).str.lower().str.contains("adjacent|non-tumor|nontumor|normal")
    rows = []
    for p, g in pd.DataFrame({"p": ph["patient_id"], "t": is_t}).groupby("p"):
        ti, ni = g.index[g["t"]], g.index[~g["t"]]
        if len(ti) == 1 and len(ni) == 1: rows.append((ti[0], ni[0]))
    T = [a for a, _ in rows]; N = [b for _, b in rows]
    d2, _ = zmean(d1src, D1 + D2_EXTRA)
    d1, _ = zmean(d1src, D1)
    c1, _ = zmean(pd.read_csv(c1s, sep="\t", index_col=0)[keep], C1)
    return e, T, N, (d1[T].values - d1[N].values), (d2[T].values - d2[N].values), \
           (c1[T].values - c1[N].values)

def fit(dy, covs, names):
    X = np.column_stack([np.ones(len(dy))] + list(covs))
    return ols_ci(np.asarray(dy, float), X, ["intercept"] + names)["intercept"]

out = {"plan": "PREREGISTRATION 6o"}

print("=== 6o-A: is DRAIN's D2 collapse attributable to POR? ===")
print("    (if DRAIN' without POR also collapses, the POR explanation is falsified)\n")
a = {}
for c in ("GSE76427", "GSE14520"):
    e, T, N, dD1, dD2, dC1 = frame(c)
    row = {}
    for tag, genes in (("DRAIN_locked", cfg.MODULE_DRAIN), ("DRAIN_prime", ["MTARC1", "MTARC2"])):
        s, used = zmean(e, genes)
        d = s[T].values - s[N].values
        un = float(d.mean())
        f1 = fit(d, [dD1], ["D1"]); f2 = fit(d, [dD2], ["D2"])
        row[tag] = {"genes": used, "unadjusted_mean": round(un, 4),
                    "D1": {"beta": f1["beta"], "p": f1["p"], "ret": round(f1["beta"]/un, 4)},
                    "D2": {"beta": f2["beta"], "p": f2["p"], "ret": round(f2["beta"]/un, 4)}}
        print(f"  {c} {tag:13s} unadj {un:+.3f} | D1 {f1['beta']:+.3f} "
              f"({100*f1['beta']/un:5.1f}%) P={f1['p']:.3g} | "
              f"D2 {f2['beta']:+.3f} ({100*f2['beta']/un:5.1f}%) P={f2['p']:.3g}"
              f"  {'** D2 NOT significant' if f2['p'] >= 0.05 else ''}")
    a[c] = row
    print()
lost_locked = [c for c in a if a[c]["DRAIN_locked"]["D2"]["p"] >= 0.05]
lost_prime  = [c for c in a if a[c]["DRAIN_prime"]["D2"]["p"] >= 0.05]
if set(lost_prime) == set(lost_locked) and lost_locked:
    verdict_a = ("FALSIFIED - DRAIN' loses significance under D2 in the same cohorts "
                 "as the locked DRAIN, so the collapse is not attributable to POR")
elif not lost_prime and lost_locked:
    verdict_a = "SUPPORTED - only the POR-containing module collapses under D2"
else:
    verdict_a = "INDETERMINATE"
print(f"  VERDICT 6o-A: {verdict_a}")
out["A_por_explains_D2_collapse"] = {"per_cohort": a, "verdict": verdict_a}

print("\n=== 6o-B: does dropping AIFM2 raise or lower REDUCTION retention? ===")
print("    (tested in GSE76427, the paired cohort where AIFM2 is measured)\n")
e, T, N, dD1, dD2, dC1 = frame("GSE76427")
b = {}
for tag, genes in (("REDUCTION_4gene", cfg.MODULE_REDUCTION),
                   ("REDUCTION_3gene_noAIFM2", [g for g in cfg.MODULE_REDUCTION if g != "AIFM2"])):
    s, used = zmean(e, genes)
    d = s[T].values - s[N].values
    un = float(d.mean())
    fj = fit(d, [dD1, dC1], ["D1", "C1"]); f1 = fit(d, [dD1], ["D1"])
    b[tag] = {"genes": used, "unadjusted_mean": round(un, 4),
              "D1_ret": round(f1["beta"]/un, 4), "joint_ret": round(fj["beta"]/un, 4)}
    print(f"  {tag:26s} ({len(used)} genes) unadj {un:+.3f} | "
          f"D1 retention {100*f1['beta']/un:5.1f}% | joint retention {100*fj['beta']/un:5.1f}%")
d_d1 = 100*(b["REDUCTION_3gene_noAIFM2"]["D1_ret"] - b["REDUCTION_4gene"]["D1_ret"])
d_jt = 100*(b["REDUCTION_3gene_noAIFM2"]["joint_ret"] - b["REDUCTION_4gene"]["joint_ret"])
print(f"\n  dropping AIFM2 changes retention by {d_d1:+.1f} pp (D1) and {d_jt:+.1f} pp (joint)")
if abs(d_d1) < 2 and abs(d_jt) < 2:
    verdict_b = ("INDETERMINABLE - both shifts under 2 pp; the manuscript's directional "
                 "claim is not supported and must be replaced by 'direction unknown'")
elif d_d1 < 0 and d_jt < 0:
    verdict_b = "CONSERVATIVE SUPPORTED - dropping AIFM2 lowers retention"
elif d_d1 > 0 and d_jt > 0:
    verdict_b = ("CLAIM WRONG - dropping AIFM2 RAISES retention, so the three-gene score "
                 "is generous, not conservative")
else:
    verdict_b = "INDETERMINATE - the two models disagree in sign"
print(f"  VERDICT 6o-B: {verdict_b}")
out["B_aifm2_direction"] = {"GSE76427": b, "delta_pp_D1": round(d_d1, 2),
                            "delta_pp_joint": round(d_jt, 2), "verdict": verdict_b}

json.dump(out, open("results/DIAGNOSTICS_6o.json", "w"), indent=1)
print("\nwrote results/DIAGNOSTICS_6o.json")
