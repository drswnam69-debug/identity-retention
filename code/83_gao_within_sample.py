#!/usr/bin/env python3
"""83_gao_within_sample.py

The cross-level comparison in this project is between cohorts: transcript level
comes from GSE14520, GSE76427 and TCGA-LIHC, protein level from Gao, Jiang and
Yi. Measurement level is therefore confounded with patients, etiology, platform
and era, and "the same quantity means different things at the two levels" cannot
be separated from "these are different people measured differently."

Gao 2019 is the one cohort measured at both levels on the same tissue. Its
deposit does NOT include the messenger RNA matrix, so the paired shift cannot be
recomputed at transcript level there. What it does include is the study's own
per-gene messenger RNA to protein Spearman correlation across those same 159
paired cases, and that is enough for one sharp test.

THE TEST, stated before it is run.
A gene whose two measurements move the same way across these samples has a high
positive correlation. A gene whose two measurements move in OPPOSITE directions
across these samples must have a correlation near zero or below it. So if the
panel's plasma proteins really reverse between levels while its metabolic
enzymes do not, then within this one cohort:

  P1  the plasma class should have a markedly lower correlation than the
      metabolic class;
  P2  the plasma class should sit near or below zero, not merely lower;
  P3  the metabolic class should sit at or above the genome-wide median.

If P1 fails, the class-structured reversal is not visible within a single
cohort and the cross-cohort result is not separable from a cohort difference.
If P1 holds but P2 fails, the classes differ in how tightly the two levels
track, which is a weaker claim than reversal.

Writes results/GAO_WITHIN_SAMPLE.json.
"""
import importlib.util, json, os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("p70", f"{HERE}/70_panel_coherence.py")
p70 = importlib.util.module_from_spec(spec); spec.loader.exec_module(p70)

ALIAS = {"G6PC1": "G6PC", "MTARC1": "MARC1"}


def main() -> int:
    corr = pd.read_csv(f"{ROOT}/data/gao_mrna_protein_corr.tsv", sep="\t")
    corr["spearman"] = pd.to_numeric(corr["spearman"], errors="coerce")
    rho = dict(zip(corr["gene"], corr["spearman"]))
    allv = corr["spearman"].dropna().to_numpy()

    def get(g):
        return rho.get(g, rho.get(ALIAS.get(g)))

    out = {"what": "within-cohort cross-level check in Gao 2019; the test was "
                   "stated in this script's docstring before it was run",
           "source": "data/gao_mrna_protein_corr.tsv, the source study's own "
                     "per-gene Spearman correlation across the same 159 pairs",
           "background": {
               "n_genes": int(allv.size),
               "median": float(np.median(allv)),
               "iqr": [float(np.percentile(allv, 25)), float(np.percentile(allv, 75))],
               "min": float(allv.min()), "max": float(allv.max()),
               "share_below_zero": float((allv < 0).mean())}}

    print(f"background: {allv.size} genes, median rho {np.median(allv):+.3f}, "
          f"range {allv.min():+.3f} to {allv.max():+.3f}, "
          f"{100*(allv<0).mean():.2f}% below zero\n")

    by_class = {}
    for cls, genes in p70.CLASS.items():
        vals, members = [], {}
        for g in genes:
            v = get(g)
            if v is not None and np.isfinite(v):
                vals.append(float(v))
                members[g] = {"rho": float(v),
                              "percentile": float(100 * (allv < v).mean())}
        if not vals:
            continue
        by_class[cls] = {
            "n": len(vals), "members": members,
            "median_rho": float(np.median(vals)),
            "mean_rho": float(np.mean(vals)),
            "median_percentile": float(np.median(
                [m["percentile"] for m in members.values()])),
            "n_below_zero": int(sum(v < 0 for v in vals))}
        print(f"{cls:26s} n={len(vals):2d}  median rho {np.median(vals):+.3f}  "
              f"(percentile {by_class[cls]['median_percentile']:5.1f})  "
              f"below zero: {by_class[cls]['n_below_zero']}")
        for g, m in sorted(members.items(), key=lambda kv: kv[1]["rho"]):
            print(f"     {g:9s} {m['rho']:+.3f}   {m['percentile']:5.1f}th")

    out["by_class"] = by_class

    pl = [m["rho"] for m in by_class["secreted plasma protein"]["members"].values()]
    me = [m["rho"] for m in by_class["metabolic enzyme"]["members"].values()]
    u = stats.mannwhitneyu(pl, me, alternative="two-sided")
    out["P1_plasma_vs_metabolic"] = {
        "n_plasma": len(pl), "n_metabolic": len(me),
        "median_plasma": float(np.median(pl)), "median_metabolic": float(np.median(me)),
        "difference_of_medians": float(np.median(pl) - np.median(me)),
        "mannwhitney_p": float(u.pvalue),
        "min_attainable_p": float(2.0 / __import__("math").comb(len(pl) + len(me), min(len(pl), len(me)))),
        "holds": bool(np.median(pl) < np.median(me))}
    # P2: is the plasma class near or below zero?
    w = stats.wilcoxon(pl) if len(pl) > 5 else None
    out["P2_plasma_near_zero"] = {
        "median": float(np.median(pl)), "min": float(min(pl)), "max": float(max(pl)),
        "n_below_zero": int(sum(v < 0 for v in pl)),
        "wilcoxon_p_against_zero": float(w.pvalue) if w is not None else None,
        "median_percentile_in_background": by_class["secreted plasma protein"]["median_percentile"],
        "holds": bool(np.median(pl) <= 0.1)}
    out["P3_metabolic_at_or_above_median"] = {
        "median": float(np.median(me)),
        "background_median": float(np.median(allv)),
        "holds": bool(np.median(me) >= np.median(allv))}

    print(f"\nP1 plasma below metabolic: {out['P1_plasma_vs_metabolic']['holds']}  "
          f"({np.median(pl):+.3f} vs {np.median(me):+.3f}, "
          f"Mann-Whitney P={u.pvalue:.3g}, floor {out['P1_plasma_vs_metabolic']['min_attainable_p']:.3g})")
    print(f"P2 plasma near or below zero: {out['P2_plasma_near_zero']['holds']}  "
          f"(median {np.median(pl):+.3f}, {out['P2_plasma_near_zero']['n_below_zero']} of "
          f"{len(pl)} below zero)")
    print(f"P3 metabolic at or above the background median: "
          f"{out['P3_metabolic_at_or_above_median']['holds']}  "
          f"({np.median(me):+.3f} vs {np.median(allv):+.3f})")

    with open(f"{ROOT}/results/GAO_WITHIN_SAMPLE.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"\nwrote {ROOT}/results/GAO_WITHIN_SAMPLE.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
