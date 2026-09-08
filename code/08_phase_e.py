#!/usr/bin/env python3
"""08_phase_e.py -- H4 (RSI and survival) and H6 (ferroptosis signature).

This is the first step that reads an outcome variable. The model was fixed in
PREREGISTRATION 6d before any survival value was fetched.

Example
-------
python3 code/08_phase_e.py --cohort GSE76427 --outdir results \
        --ferroptosis data/ferroptosis_sets.tsv \
        --h6-cohorts GSE135251 GSE130970 GSE167523 GSE164760
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg      # noqa: E402
import stats_lite as sl       # noqa: E402
import survival as sv         # noqa: E402


def zscore(v: pd.Series) -> pd.Series:
    return (v - v.mean()) / v.std(ddof=1)


def load_sets(path: str) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            name, genes = line.rstrip("\n").split("\t", 1)
            sets[name] = [g.strip().upper() for g in genes.split(",")
                          if g.strip()]
    return sets


def module_score(expr: pd.DataFrame, genes: list[str]):
    present = [g for g in genes if g in expr.index]
    if len(present) < 5:
        return None, present
    sub = expr.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0).dropna(how="all")
    return z.mean(axis=0), present


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="GSE76427")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--ferroptosis", default=None)
    ap.add_argument("--h6-cohorts", nargs="*", default=[])
    ap.add_argument("--skip-h4", action="store_true",
                    help="run H6 only; H4 needs the prognostic cohort")
    args = ap.parse_args()

    res: dict = {"lock_hash": cfg.lock_hash(),
                 "plan": "PREREGISTRATION 6d, fixed before any outcome was read"}

    # ================================================================ H4
    if args.skip_h4:
        h4 = None
    else:
        h4 = _run_h4(args, res)

    # ================================================================ H6
    _run_h6(args, res)

    out_path = os.path.join(args.outdir,
                            "phase_e_h6.json" if args.skip_h4
                            else "phase_e.json")
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    if "H6" in res:
        print("\n  H6  ferroptosis signatures (RSI-panel genes removed)")
        for name, meta in res["H6"]["sets"].items():
            print(f"    {name}: {meta['n_after_overlap_removal']} genes "
                  f"after dropping {meta['dropped_overlapping_with_RSI_panel']}")
        for c, rec in res["H6"]["cohorts"].items():
            bits = [f"{k.replace('GOBP_','').replace('_REGULATION_OF_FERROPTOSIS','')[:12]} "
                    f"rho={v['rho']:+.3f} P={v['p']:.2g}"
                    for k, v in rec.items() if "rho" in v]
            print(f"    {c:<12} " + "   ".join(bits))
    print(f"\n  wrote -> {out_path}")
    return 0


def _run_h4(args, res):
    d = os.path.join(args.outdir, args.cohort)
    scores = pd.read_csv(os.path.join(d, "rsi.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(d, "phenotype.tsv"), sep="\t",
                        index_col=0).loc[scores.index]
    expr = pd.read_csv(os.path.join(d, "expr_log.tsv.gz"), sep="\t",
                       index_col=0)[scores.index]

    tissue = pheno["tissue"].astype(str).str.lower()
    tum = tissue.str.contains("tumor") & ~tissue.str.contains("non-tumor")
    df = pd.DataFrame({
        "rsi": scores.loc[tum, "RSI"],
        "time": pd.to_numeric(pheno.loc[tum, "duryears_os"], errors="coerce"),
        "event": pd.to_numeric(pheno.loc[tum, "event_os"], errors="coerce"),
        "age": pd.to_numeric(pheno.loc[tum, "age"], errors="coerce"),
        "sex": pd.to_numeric(pheno.loc[tum, "sex"], errors="coerce"),
        "bclc": pheno.loc[tum, "bclc"].astype(str),
        "mean_expr": expr.loc[:, tum.index[tum]].mean(axis=0),
        "patient": pheno.loc[tum, "patient_id"].astype(str)
                   if "patient_id" in pheno else pheno.index[tum],
    }).dropna(subset=["rsi", "time", "event"])
    df = df[df["time"] > 0]
    df = df.drop_duplicates(subset="patient", keep="first")

    rsi_sd = zscore(df["rsi"])
    h4: dict = {"cohort": args.cohort, "n": int(len(df)),
                "n_events": int(df["event"].sum()),
                "endpoint": "overall survival",
                "note": "one sample per patient; tumor specimens only"}

    m = sv.cox_ph(df["time"], df["event"], rsi_sd.to_numpy()[:, None],
                  names=["RSI_per_SD"])
    h4["primary_continuous"] = {
        "HR_per_SD": round(float(m["hr"][0]), 4),
        "ci95": [round(float(m["hr_lo"][0]), 4), round(float(m["hr_hi"][0]), 4)],
        "p": float(m["p"][0]), "n": m["n"], "n_events": m["n_events"]}
    h4["ph_check"] = sv.schoenfeld_ph_test(df["time"], df["event"],
                                           rsi_sd.to_numpy(),
                                           float(m["coef"][0]))

    # tertiles, cut at the cohort's own 33rd and 67th percentiles
    q1, q2 = df["rsi"].quantile([1 / 3, 2 / 3])
    tert = pd.cut(df["rsi"], [-np.inf, q1, q2, np.inf],
                  labels=["T1", "T2", "T3"])
    h4["tertiles"] = {
        "cutpoints": [round(float(q1), 4), round(float(q2), 4)],
        "n": tert.value_counts().reindex(["T1", "T2", "T3"]).to_dict(),
        "events": {k: int(df.loc[tert == k, "event"].sum())
                   for k in ("T1", "T2", "T3")},
        "median_survival": {k: round(sv.median_survival(
            df.loc[tert == k, "time"], df.loc[tert == k, "event"]), 3)
            for k in ("T1", "T2", "T3")},
        "logrank": sv.logrank(df["time"], df["event"], tert.astype(str)),
    }

    # adjusted model, reported with an explicit events-per-variable warning
    cov = [rsi_sd.to_numpy()]
    names = ["RSI_per_SD"]
    for col, label in (("age", "age"), ("sex", "sex_female")):
        v = df[col]
        if v.notna().mean() > 0.8:
            cov.append((v == 2).astype(float).to_numpy() if col == "sex"
                       else zscore(v).fillna(0).to_numpy())
            names.append(label)
    bclc = df["bclc"].replace({".": np.nan})
    if bclc.notna().mean() > 0.8:
        cov.append((bclc.isin(["C", "D"])).astype(float).to_numpy())
        names.append("BCLC_C_or_D")
    try:
        ma = sv.cox_ph(df["time"], df["event"], np.column_stack(cov),
                       names=names)
        h4["adjusted"] = {
            "terms": ma["names"],
            "HR": [round(float(x), 4) for x in ma["hr"]],
            "ci95": [[round(float(a), 4), round(float(b), 4)]
                     for a, b in zip(ma["hr_lo"], ma["hr_hi"])],
            "p": [float(x) for x in ma["p"]],
            "events_per_variable": round(ma["n_events"] / len(names), 2),
            "warning": ("fewer than 10 events per variable; this model is "
                        "underpowered and is secondary"
                        if ma["n_events"] / len(names) < 10 else None)}
    except Exception as exc:
        h4["adjusted"] = {"error": str(exc)}

    # global-axis covariate, per PREREGISTRATION 6d
    rho_ax, p_ax = sl.spearman(df["rsi"].to_numpy(float),
                               df["mean_expr"].to_numpy(float))
    h4["global_axis"] = {"rho": round(rho_ax, 4), "p": p_ax}
    if abs(rho_ax) > 0.2:
        mg = sv.cox_ph(df["time"], df["event"],
                       np.column_stack([rsi_sd.to_numpy(),
                                        zscore(df["mean_expr"]).to_numpy()]),
                       names=["RSI_per_SD", "mean_expr_z"])
        h4["adjusted_for_global_axis"] = {
            "HR_per_SD": round(float(mg["hr"][0]), 4),
            "ci95": [round(float(mg["hr_lo"][0]), 4),
                     round(float(mg["hr_hi"][0]), 4)],
            "p": float(mg["p"][0])}
    res["H4"] = h4
    p = h4["primary_continuous"]
    print(f"\n=== {args.cohort} — PHASE E ===")
    print(f"  H4  n={p['n']} tumors, {p['n_events']} deaths")
    print(f"    Cox, RSI per SD:  HR = {p['HR_per_SD']:.3f}  "
          f"(95% CI {p['ci95'][0]:.3f}-{p['ci95'][1]:.3f})   P = {p['p']:.3g}")
    ph = h4["ph_check"]
    print(f"    proportional hazards: rho={ph['rho']}, P={ph['p']:.3g}")
    t = h4["tertiles"]
    print(f"    tertiles n={t['n']}  events={t['events']}")
    print(f"    log-rank across tertiles: P = {t['logrank']['p']:.3g}")
    if "adjusted" in h4 and "terms" in h4["adjusted"]:
        a = h4["adjusted"]
        print(f"    adjusted ({', '.join(a['terms'])}), "
              f"{a['events_per_variable']} events/variable:")
        for n_, hr, ci, pv in zip(a["terms"], a["HR"], a["ci95"], a["p"]):
            print(f"      {n_:<14} HR {hr:.3f} ({ci[0]:.3f}-{ci[1]:.3f})  P={pv:.3g}")
    return h4


def _run_h6(args, res):
    if args.ferroptosis and os.path.exists(args.ferroptosis):
        sets = load_sets(args.ferroptosis)
        panel = {g for gs in cfg.ALL_MODULES.values() for g in gs}
        h6: dict = {"sets": {}, "cohorts": {}}
        cleaned = {}
        for name, genes in sets.items():
            dropped = sorted(set(genes) & panel)
            cleaned[name] = [g for g in genes if g not in panel]
            h6["sets"][name] = {"n_original": len(genes),
                                "n_after_overlap_removal": len(cleaned[name]),
                                "dropped_overlapping_with_RSI_panel": dropped}
        for c in (args.h6_cohorts or [args.cohort]):
            path = os.path.join(args.outdir, c, "expr_log.tsv.gz")
            if not os.path.exists(path):
                continue
            e = pd.read_csv(path, sep="\t", index_col=0)
            sc = pd.read_csv(os.path.join(args.outdir, c, "rsi.tsv"),
                             sep="\t", index_col=0)
            e = e[sc.index]
            rec = {}
            for name, genes in cleaned.items():
                s, present = module_score(e, genes)
                if s is None:
                    rec[name] = {"skipped": f"only {len(present)} genes present"}
                    continue
                rho, p = sl.spearman(s.to_numpy(float),
                                     sc["RSI"].to_numpy(float))
                rec[name] = {"n_genes": len(present), "rho": round(rho, 4),
                             "p": p, "n": int(len(sc))}
            h6["cohorts"][c] = rec
        res["H6"] = h6


if __name__ == "__main__":
    raise SystemExit(main())
