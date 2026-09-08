#!/usr/bin/env python3
"""40_prognostic_value.py -- PREREGISTRATION 6z.

Does the identity-retention fraction say anything about whether a signature's
prognostic value is independent of the tumor's differentiation state?

Retention comes from the paired samples (already archived under 6t/6w). Survival
comes from every primary tumor with an outcome. The two estimation sets are
different by construction, and neither informs the other.

Everything below is fixed in 6z: the cohort, the endpoint, the signature set, the
adjustment variable, the models, the quantity of interest, the power rule and the
outcomes accepted.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(HERE, f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


_da = _load("d", "12_differentiation_adjust.py")
_bm = _load("b", "30_signature_benchmark.py")
D1, zmean = _da.D1, _bm.zmean
import survival as sv  # noqa: E402

MIN_EVENTS, MIN_PROGNOSTIC = 30, 10       # 6z power rule
RHO_THRESHOLD = 0.30                      # 6z decision threshold


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in enumerate(o[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def main() -> None:
    surv = pd.read_csv("data/TCGA-LIHC.survival.tsv.gz", sep="\t")
    surv = surv.dropna(subset=["OS", "OS.time"])
    surv = surv[surv["sample"].astype(str).str.split("-").str[3].str[:2] == "01"]
    surv = surv.drop_duplicates(subset=["_PATIENT"])
    expr = pd.read_csv("data/TCGA_LIHC_TUMORS_symbols.tsv.gz", sep="\t", index_col=0)
    keep = [s for s in surv["sample"] if s in expr.columns]
    surv = surv.set_index("sample").loc[keep]
    expr = expr[keep]
    time = surv["OS.time"].to_numpy(float)
    event = surv["OS"].to_numpy(int)
    n, n_ev = len(keep), int(event.sum())
    print("=== TCGA-LIHC prognostic analysis (PREREG 6z) ===")
    print(f"  {n} primary tumors with overall survival, {n_ev} deaths")
    if n_ev < MIN_EVENTS:
        raise SystemExit("fewer than 30 deaths; reported as underpowered under 6z")

    d1, n_d1 = zmean(expr, D1)
    d1z = ((d1 - d1.mean()) / d1.std()).to_numpy(float)
    fit_d1 = sv.cox_ph(time, event, d1z.reshape(-1, 1), ["D1"])
    print(f"\n  [1] does identity itself predict outcome? D1 on {n_d1} genes")
    hr = float(fit_d1["hr"][0]); pv = float(fit_d1["p"][0])
    lo, hi = float(fit_d1["hr_lo"][0]), float(fit_d1["hr_hi"][0])
    ph = sv.schoenfeld_ph_test(time, event, d1z, float(fit_d1["coef"][0]))
    print(f"      HR per SD = {hr:.3f} ({lo:.3f} to {hi:.3f}), "
          f"P = {pv:.3g};  Schoenfeld rho = {ph['rho']:.3f}, P = {ph['p']:.3g}")

    arch = json.load(open("results/SIGNATURE_BENCHMARK_TCGA_LIHC.json"))
    sets = json.load(open("genesets/eligible.json"))
    ok = [r for r in arch["signatures"] if r["status"] == "ok"]
    print(f"\n  [2] {len(ok)} signatures, the evaluable set already reported under 6t")

    rows = []
    for r in ok:
        genes = sets[r["name"]]
        s, n_on = zmean(expr, genes)
        if s is None:
            continue
        z = ((s - s.mean()) / s.std()).to_numpy(float)
        f_un = sv.cox_ph(time, event, z.reshape(-1, 1), ["sig"])
        f_ad = sv.cox_ph(time, event, np.column_stack([z, d1z]), ["sig", "D1"])
        def grab(f, i=0):
            return float(f["coef"][i]), float(f["p"][i])
        b_un, p_un = grab(f_un)
        b_ad, p_ad = grab(f_ad)
        rows.append({"name": r["name"], "retention_joint": r["retention_joint"],
                     "n_on_platform": n_on,
                     "logHR_unadj": round(b_un, 4), "p_unadj": p_un,
                     "HR_unadj": round(float(np.exp(b_un)), 4),
                     "logHR_adj": round(b_ad, 4), "p_adj": p_ad,
                     "HR_adj": round(float(np.exp(b_ad)), 4),
                     "surviving_fraction": (round(abs(b_ad) / abs(b_un), 4)
                                            if abs(b_un) > 1e-8 else None)})
    df = pd.DataFrame(rows)
    df["q_unadj"] = bh(df["p_unadj"].values)
    df["q_adj"] = bh(df["p_adj"].values)
    prog = df[df["p_unadj"] < 0.05]
    print(f"      prognostic before adjustment at P < 0.05: {len(prog)}")
    print(f"      surviving Benjamini-Hochberg at 5%:       {int((df['q_unadj']<0.05).sum())}")
    if len(prog) < MIN_PROGNOSTIC:
        raise SystemExit("fewer than 10 prognostic signatures; underpowered under 6z")
    lost = prog[prog["p_adj"] >= 0.05]
    print(f"      of the {len(prog)}, losing significance after adjusting for D1: {len(lost)}")

    sub = df[df["surviving_fraction"].notna()]
    rho, p_rho = spearmanr(sub["retention_joint"], sub["surviving_fraction"])
    print(f"\n  [3] PRIMARY: Spearman rho between retention and the surviving "
          f"fraction of log HR")
    print(f"      rho = {rho:+.3f}  (P = {p_rho:.3g})  over {len(sub)} signatures")
    med = df["retention_joint"].median()
    hi_r = prog[prog["retention_joint"] >= med]; lo_r = prog[prog["retention_joint"] < med]
    print(f"      split at the median retention {med:.3f}:")
    for lab, g in (("high retention", hi_r), ("low retention", lo_r)):
        if len(g):
            print(f"        {lab:<15} {len(g)} prognostic, "
                  f"{int((g['p_adj']>=0.05).sum())} lose significance "
                  f"({(g['p_adj']>=0.05).mean()*100:.0f}%), "
                  f"median surviving fraction {g['surviving_fraction'].median():.3f}")

    if rho >= RHO_THRESHOLD and p_rho < 0.05:
        verdict = "ASSOCIATED"
    elif rho <= -RHO_THRESHOLD and p_rho < 0.05:
        verdict = "INVERSE"
    else:
        verdict = "NOT_ASSOCIATED"
    print(f"\n  [4] 6z verdict: {verdict}")

    res = {"plan": "PREREGISTRATION 6z",
           "amendment_sha256_of_text_as_written":
               "8252ff3c8bb3d7fbf18e054b3c3bcee0c68febf6defa97174e7f1773ffca9fc7",
           "cohort": "TCGA_LIHC tumors", "endpoint": "overall survival",
           "n_tumors": n, "n_deaths": n_ev,
           "D1": {"n_genes": n_d1, "HR_per_SD": round(hr, 4),
                  "ci95": [round(lo, 4), round(hi, 4)],
                  "p": pv, "schoenfeld_rho": ph["rho"], "schoenfeld_p": ph["p"]},
           "n_signatures": len(df),
           "n_prognostic_p05": int(len(prog)),
           "n_prognostic_bh05": int((df["q_unadj"] < 0.05).sum()),
           "n_losing_significance_after_adjustment": int(len(lost)),
           "median_retention": round(float(med), 4),
           "spearman_rho": round(float(rho), 4), "spearman_p": float(p_rho),
           "rho_threshold": RHO_THRESHOLD, "verdict": verdict,
           "signatures": df.to_dict("records")}
    with open("results/PROGNOSTIC_6z_TCGA_LIHC.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print("\n  wrote results/PROGNOSTIC_6z_TCGA_LIHC.json")


if __name__ == "__main__":
    main()
