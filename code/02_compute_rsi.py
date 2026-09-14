#!/usr/bin/env python3
"""02_compute_rsi.py -- compute RSI and run PHASE A internal-validity checks.

This script deliberately does NOT touch any outcome variable (survival,
tumor status). It reports only the internal properties of the index, so that
the definition can be locked before outcomes are opened.

Example
-------
python3 code/02_compute_rsi.py --cohort GSE135251 --outdir results
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg   # noqa: E402
import rsi                 # noqa: E402
import stats_lite as sl    # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--sensitivity-cholesterol", action="store_true",
                    help="sensitivity S2: add cholesterol-arm genes to SUPPLY")
    args = ap.parse_args()

    out = os.path.join(args.outdir, args.cohort)
    expr = pd.read_csv(os.path.join(out, "expr_log.tsv.gz"), sep="\t",
                       index_col=0)

    scores, report = rsi.compute_rsi(
        expr, include_cholesterol=args.sensitivity_cholesterol)

    suffix = "_S2chol" if args.sensitivity_cholesterol else ""
    scores.to_csv(os.path.join(out, f"rsi{suffix}.tsv"), sep="\t")

    # --- internal validity (no outcomes) -----------------------------------
    modules_used = dict(cfg.ALL_MODULES)
    if args.sensitivity_cholesterol:
        modules_used["SUPPLY"] = (list(cfg.MODULE_SUPPLY)
                                  + list(cfg.SENSITIVITY_S2_CHOLESTEROL))

    validity = {}
    for name in ("SUPPLY", "REDUCTION", "DRAIN"):
        genes = modules_used[name]
        validity[name] = {
            "n_genes_specified": len(genes),
            "n_genes_present": len([g for g in genes if g in expr.index]),
            "cronbach_alpha": round(rsi.cronbach_alpha(expr, genes), 4),
            "pc1_variance_explained": round(
                rsi.pc1_variance_explained(expr, genes), 4),
        }

    comp = scores[["z_supply", "z_reduction", "z_drain", "RSI"]]
    validity["component_correlation_spearman"] = (
        comp.rank().corr().round(3).to_dict())  # rank-Pearson == Spearman
    validity["rsi_summary"] = {
        "mean": round(float(scores["RSI"].mean()), 4),
        "sd": round(float(scores["RSI"].std(ddof=1)), 4),
        "min": round(float(scores["RSI"].min()), 4),
        "max": round(float(scores["RSI"].max()), 4),
    }
    # Is the index just re-reading each sample's overall expression level?
    # Every module z-score is vulnerable to a global axis; the point of
    # subtracting the drain is that it should cancel. This checks that it does.
    mean_expr = expr[scores.index].mean(axis=0)
    validity["global_axis_check"] = {
        "note": "Spearman rho against each sample's mean expression. Module "
                "scores are expected to correlate; RSI should not, because "
                "the drain subtraction cancels the shared component.",
        **{k: round(sl.spearman(scores[k].to_numpy(float),
                                mean_expr.to_numpy(float))[0], 4)
           for k in ("z_supply", "z_reduction", "z_drain", "RSI")},
    }

    validity["variant"] = ("S2_cholesterol_added"
                           if args.sensitivity_cholesterol else "locked")
    validity["lock_date"] = cfg.LOCK_DATE
    validity["lock_hash"] = cfg.lock_hash()
    validity["missing_genes"] = {
        "supply": report["supply_missing"],
        "reduction": report["reduction_missing"],
        "drain": report["drain_missing"],
    }

    with open(os.path.join(out, f"internal_validity{suffix}.json"), "w") as fh:
        json.dump(validity, fh, indent=2)

    print(f"\n=== {args.cohort} — PHASE A internal validity ===")
    for name in ("SUPPLY", "REDUCTION", "DRAIN"):
        v = validity[name]
        print(f"  {name:<10} {v['n_genes_present']}/{v['n_genes_specified']} "
              f"genes  alpha={v['cronbach_alpha']:.3f}  "
              f"PC1={v['pc1_variance_explained']:.1%}")
    miss = {k: v for k, v in validity["missing_genes"].items() if v}
    if miss:
        print(f"  ! genes absent from the matrix: {miss}")
    print(f"  RSI  mean={validity['rsi_summary']['mean']:.3f} "
          f"sd={validity['rsi_summary']['sd']:.3f} "
          f"range=[{validity['rsi_summary']['min']:.2f}, "
          f"{validity['rsi_summary']['max']:.2f}]")
    ga = validity["global_axis_check"]
    print(f"  global-axis check (rho vs mean expression): "
          f"supply {ga['z_supply']:+.2f}  reduction {ga['z_reduction']:+.2f}  "
          f"drain {ga['z_drain']:+.2f}  ->  RSI {ga['RSI']:+.2f}")
    print(f"  lock hash {cfg.lock_hash()[:16]}...")
    print(f"  wrote -> {out}/rsi{suffix}.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
