#!/usr/bin/env python3
"""03_phase_b.py -- H1: does RSI rise monotonically across the MASLD spectrum?

Primary test: Jonckheere-Terpstra for an ordered alternative, which is the
test that matches the hypothesis ("monotone increase"), rather than a general
Kruskal-Wallis. Sub-scores and individual genes are reported alongside, with
BH-FDR across the gene panel.

No outcome variable (survival, tumor status) is touched here.

Example
-------
python3 code/03_phase_b.py --cohort GSE135251 --outdir results
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg      # noqa: E402
import phenotype as ph        # noqa: E402
import stats_lite as sl       # noqa: E402


def cohort_axis(outdir: str, cohort: str):
    """The ordered severity ladder this cohort was prepared with."""
    with open(os.path.join(outdir, cohort, "qc.json")) as fh:
        qc = json.load(fh)
    return qc.get("stage_order", ph.STAGE_ORDER), qc.get("axis", "group")


def fisher_ci(rho: float, n: int, conf: float = 0.95):
    """95% CI for a Spearman rho via the Fisher z transform."""
    if not np.isfinite(rho) or n < 5 or abs(rho) >= 1:
        return (float("nan"), float("nan"))
    z = math.atanh(rho)
    se = 1.06 / math.sqrt(n - 3)          # Bonett-Wright SE for Spearman
    crit = 1.959963984540054 if conf == 0.95 else 1.959963984540054
    return (math.tanh(z - crit * se), math.tanh(z + crit * se))


def ordered_groups(values: pd.Series, stage: pd.Series, order):
    """Split values by stage, returned in the pre-specified order."""
    present = [s for s in order if (stage == s).any()]
    return present, [values[stage == s].to_numpy(float) for s in present]


def trend(values: pd.Series, stage: pd.Series, order) -> dict:
    labels, groups = ordered_groups(values, stage, order)
    jt, z, p = sl.jonckheere_terpstra(groups, alternative="increasing")
    rank = stage.map({s: i for i, s in enumerate(order)}).astype(float)
    rho, p_rho = sl.spearman(values.to_numpy(float), rank.to_numpy(float))
    lo, hi = fisher_ci(rho, int(np.isfinite(values).sum()))
    return {
        "stages": labels,
        "n_per_stage": [int(len(g)) for g in groups],
        "median_per_stage": [round(float(np.median(g)), 4) if len(g) else None
                             for g in groups],
        "JT": round(jt, 1), "z": round(z, 4), "p_increasing": p,
        "spearman_rho": round(rho, 4), "spearman_p": p_rho,
        "rho_ci95": [round(lo, 4), round(hi, 4)],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    out = os.path.join(args.outdir, args.cohort)
    scores = pd.read_csv(os.path.join(out, "rsi.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(out, "phenotype.tsv"), sep="\t",
                        index_col=0)
    pheno = pheno.loc[scores.index]
    order, axis = cohort_axis(args.outdir, args.cohort)
    stage = pheno["stage"].astype(str)
    keep = stage.isin(order)
    if keep.sum() < len(stage):
        print(f"  ! dropping {int((~keep).sum())} samples with no stage")
    scores, pheno, stage = scores[keep], pheno[keep], stage[keep]

    res: dict = {"cohort": args.cohort, "n": int(len(scores)),
                 "axis": axis, "stage_order": order,
                 "lock_hash": cfg.lock_hash(), "hypothesis": cfg.HYPOTHESES["H1"]}

    # ---- H1, the pre-specified primary test ----------------------------
    res["primary_RSI"] = trend(scores["RSI"], stage, order)

    # ---- which component moves the index -------------------------------
    res["components"] = {c: trend(scores[c], stage, order)
                         for c in ("z_supply", "z_reduction", "z_drain")}

    # ---- context modules, reported but not part of the index -----------
    res["modules"] = {c: trend(scores[c], stage, order)
                      for c in ("z_shared", "z_effector", "z_upstream")
                      if c in scores and scores[c].notna().any()}

    # ---- individual genes, BH-FDR across the panel ---------------------
    gene_cols = [c for c in scores.columns if c.startswith("gene_")]
    gene_res = {c[5:]: trend(scores[c], stage, order) for c in gene_cols}
    if gene_res:
        names = list(gene_res)
        q = sl.benjamini_hochberg([gene_res[g]["p_increasing"] for g in names])
        for g, qv in zip(names, q):
            gene_res[g]["q_bh"] = float(qv)
    res["genes"] = gene_res

    # ---- covariate-adjusted model where covariates exist ---------------
    rank = stage.map({s: i for i, s in enumerate(order)}).astype(float)
    cols, names = [np.ones(len(scores)), rank.to_numpy(float)], ["const", "stage"]
    for cand, label in (("age", "age"), ("bmi", "bmi")):
        if cand in pheno and pd.to_numeric(pheno[cand], errors="coerce").notna().sum() > 0.8 * len(pheno):
            cols.append(pd.to_numeric(pheno[cand], errors="coerce").to_numpy(float))
            names.append(label)
    sexcol = next((c for c in ("sex", "gender", "Sex") if c in pheno), None)
    if sexcol is not None:
        s = pheno[sexcol].astype(str).str.lower().str[0]
        if s.isin(list("mf")).mean() > 0.8:
            cols.append((s == "m").astype(float).to_numpy())
            names.append("sex_male")
    try:
        model = sl.ols(scores["RSI"].to_numpy(float), np.column_stack(cols),
                       names=names)
        res["adjusted_model"] = {
            "terms": model["names"],
            "beta": [round(float(b), 4) for b in model["beta"]],
            "se": [round(float(s), 4) for s in model["se"]],
            "p": [float(p) for p in model["p"]],
            "n": model["n"], "r2": round(model["r2"], 4),
            "note": f"stage is coded 0-{len(order) - 1} along the "
                    f"{axis} axis",
        }
    except Exception as exc:                      # pragma: no cover
        res["adjusted_model"] = {"error": str(exc)}

    # ---- NAS components, where the cohort provides them ----------------
    nas = {}
    for cand in ("nas", "steatosis_grade", "lobular_inflammation_grade",
                 "cytological_ballooning_grade", "fibrosis"):
        if cand in pheno:
            v = pd.to_numeric(pheno[cand], errors="coerce")
            if v.notna().sum() > 10:
                rho, p = sl.spearman(scores["RSI"].to_numpy(float),
                                     v.to_numpy(float))
                nas[cand] = {"spearman_rho": round(rho, 4), "p": p,
                             "n": int(v.notna().sum())}
    if nas:
        res["histology_components"] = nas

    with open(os.path.join(out, "phase_b.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=float)

    # ---- console summary ------------------------------------------------
    pr = res["primary_RSI"]
    print(f"\n=== {args.cohort} — PHASE B (H1), axis: {axis} ===")
    print(f"  stages: {', '.join(f'{s}={n}' for s, n in zip(pr['stages'], pr['n_per_stage']))}")
    print(f"  RSI medians: {pr['median_per_stage']}")
    print(f"  Jonckheere-Terpstra (increasing): z={pr['z']:.3f}  "
          f"p={pr['p_increasing']:.3g}")
    print(f"  Spearman rho vs stage: {pr['spearman_rho']:.3f} "
          f"(95% CI {pr['rho_ci95'][0]:.3f} to {pr['rho_ci95'][1]:.3f})")
    print("  component trends (z, p):")
    for c, r in res["components"].items():
        print(f"    {c:<12} z={r['z']:>7.3f}  p={r['p_increasing']:.3g}")
    if gene_res:
        sig = [(g, r) for g, r in gene_res.items() if r.get("q_bh", 1) < 0.05]
        print(f"  genes with BH-FDR q<0.05 for an increasing trend: "
              f"{len(sig)}/{len(gene_res)}")
        for g, r in sorted(sig, key=lambda kv: kv[1]["p_increasing"])[:8]:
            print(f"    {g:<10} rho={r['spearman_rho']:>6.3f}  q={r['q_bh']:.3g}")
    print(f"  wrote -> {out}/phase_b.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
