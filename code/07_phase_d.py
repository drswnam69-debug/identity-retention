#!/usr/bin/env python3
"""07_phase_d.py -- H5: NMD-target burden versus RSI.

EXPLORATORY. The NMD Consensus method is transcript-level; collapsing to gene
level cancels the target-vs-control-isoform contrast it measures. See
PREREGISTRATION.md 6c. This does not settle H5 either way.

Higher score = more NMD substrate accumulation = LOWER inferred NMD activity.
The framework predicts burden and RSI move together (rho > 0); the registered
wording said rho < 0 and was in error, so the test is two-sided.

Example
-------
python3 code/07_phase_d.py --cohorts GSE135251 GSE130970 GSE167523 \
        --outdir results --nmd-set data/S4_NMD_Consensus_target_gene_set.txt
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stats_lite as sl       # noqa: E402
import rsi_config as cfg      # noqa: E402

MITO = ["ATP5A1", "COX11", "DAP3", "ETHE1", "HIBADH", "MRPL48", "MRPL49",
        "MRPS27", "MRPS35", "NDUFS4"]
STRESS = ["BAX", "CIRBP", "CRTC2", "GADD45B"]


def load_nmd_set(path: str) -> list[str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    tg = df[df["final_consensus"] == "NMD_target"]
    genes = sorted({g.strip().upper() for g in tg["gene_symbol"].dropna()})
    return genes


def score(expr: pd.DataFrame, genes: list[str]):
    present = [g for g in genes if g in expr.index]
    if len(present) < 10:
        return None, present
    sub = expr.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0).dropna(how="all")
    return z.mean(axis=0), present


def fisher_ci(rho: float, n: int):
    if not np.isfinite(rho) or n < 5 or abs(rho) >= 1:
        return (float("nan"), float("nan"))
    z, se = math.atanh(rho), 1.06 / math.sqrt(n - 3)
    return math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)


def corr(a: pd.Series, b: pd.Series) -> dict:
    rho, p = sl.spearman(a.to_numpy(float), b.to_numpy(float))
    lo, hi = fisher_ci(rho, int(min(a.notna().sum(), b.notna().sum())))
    return {"rho": round(rho, 4), "p_two_sided": p,
            "ci95": [round(lo, 4), round(hi, 4)], "n": int(len(a))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="+", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--nmd-set", required=True)
    args = ap.parse_args()

    genes = load_nmd_set(args.nmd_set)
    panel = {g for gs in cfg.ALL_MODULES.values() for g in gs}
    assert not (set(genes) & panel), "NMD set overlaps the RSI panel"

    out: dict = {
        "gene_set": os.path.basename(args.nmd_set),
        "n_genes_in_set": len(genes),
        "citation": "Palou-Marquez G, Supek F. Genome Biol 2025;26:316",
        "status": "EXPLORATORY - gene-level use degrades a transcript-level "
                  "method; see PREREGISTRATION 6c",
        "direction": "higher score = more substrate accumulation = lower "
                     "inferred NMD activity; framework expects rho > 0 with "
                     "RSI; tested two-sided",
        "cohorts": {},
    }

    for c in args.cohorts:
        d = os.path.join(args.outdir, c)
        expr = pd.read_csv(os.path.join(d, "expr_log.tsv.gz"), sep="\t",
                           index_col=0)
        scores = pd.read_csv(os.path.join(d, "rsi.tsv"), sep="\t", index_col=0)
        expr = expr[scores.index]

        rec: dict = {}
        s_all, present = score(expr, genes)
        rec["n_genes_matched"] = len(present)
        if s_all is None:
            rec["skipped"] = ("fewer than 10 set genes present in this "
                              f"matrix ({len(present)})")
            out["cohorts"][c] = rec
            print(f"\n=== {c} — skipped: only {len(present)} of {len(genes)} "
                  f"NMD-set genes present ===")
            continue

        rec["H5_primary"] = corr(s_all, scores["RSI"])
        rec["vs_components"] = {
            k: corr(s_all, scores[k])
            for k in ("z_supply", "z_reduction", "z_drain")}

        # pre-specified sensitivities
        s_n1, p1 = score(expr, [g for g in genes if g not in MITO])
        s_n2, p2 = score(expr, [g for g in genes if g not in STRESS])
        rec["N1_no_mito"] = {**corr(s_n1, scores["RSI"]), "n_genes": len(p1)}
        rec["N2_no_stress"] = {**corr(s_n2, scores["RSI"]), "n_genes": len(p2)}
        rec["N3_vs_machinery"] = {
            g: corr(s_all, expr.loc[g]) for g in ("SMG1", "SMG8", "SMG9",
                                                  "UPF1", "UPF2")
            if g in expr.index}

        # is the score just tracking overall expression level?
        rec["score_vs_library_mean"] = corr(s_all, expr.mean(axis=0))

        out["cohorts"][c] = rec
        pr = rec["H5_primary"]
        print(f"\n=== {c} — PHASE D (H5, exploratory) ===")
        print(f"  {len(present)}/{len(genes)} NMD-set genes matched")
        print(f"  NMD burden vs RSI:  rho={pr['rho']:+.3f}  "
              f"(95% CI {pr['ci95'][0]:+.3f} to {pr['ci95'][1]:+.3f})  "
              f"P={pr['p_two_sided']:.3g}")
        print(f"    N1 no mito/OXPHOS  rho={rec['N1_no_mito']['rho']:+.3f}  "
              f"P={rec['N1_no_mito']['p_two_sided']:.3g}")
        print(f"    N2 no ISR/stress   rho={rec['N2_no_stress']['rho']:+.3f}  "
              f"P={rec['N2_no_stress']['p_two_sided']:.3g}")
        for k, v in rec["vs_components"].items():
            print(f"    vs {k:<12} rho={v['rho']:+.3f}  P={v['p_two_sided']:.3g}")
        for g, v in rec["N3_vs_machinery"].items():
            print(f"    vs {g:<12} rho={v['rho']:+.3f}  P={v['p_two_sided']:.3g}")

    # random-effects pooling of the primary correlation
    rows = [(c, r["H5_primary"]["rho"], r["H5_primary"]["n"])
            for c, r in out["cohorts"].items() if "H5_primary" in r]
    if len(rows) > 1:
        z = np.array([math.atanh(r) for _, r, _ in rows])
        v = np.array([1.06 ** 2 / (n - 3) for _, _, n in rows])
        w = 1 / v
        zf = (w * z).sum() / w.sum()
        q = float((w * (z - zf) ** 2).sum())
        df = len(z) - 1
        cst = w.sum() - (w ** 2).sum() / w.sum()
        tau2 = max(0.0, (q - df) / cst) if cst > 0 else 0.0
        wr = 1 / (v + tau2)
        zr = (wr * z).sum() / wr.sum()
        se = math.sqrt(1 / wr.sum())
        out["pooled"] = {
            "rho": round(math.tanh(zr), 4),
            "ci95": [round(math.tanh(zr - 1.96 * se), 4),
                     round(math.tanh(zr + 1.96 * se), 4)],
            "i2_percent": round(max(0.0, (q - df) / q) * 100 if q > 0 else 0, 1),
            "n_total": int(sum(n for _, _, n in rows)),
            "cohorts": [c for c, _, _ in rows],
        }
        p = out["pooled"]
        print(f"\n  pooled over {len(rows)} cohorts (n={p['n_total']}): "
              f"rho={p['rho']:+.3f} (95% CI {p['ci95'][0]:+.3f} to "
              f"{p['ci95'][1]:+.3f}), I2={p['i2_percent']:.0f}%")

    path = os.path.join(args.outdir, "phase_d.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=float)
    print(f"\n  wrote -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
