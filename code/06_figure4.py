#!/usr/bin/env python3
"""06_figure4.py -- manuscript Figure 4: the HCC transition and the field effect.

Panel A  RSI across the five GSE164760 tissue groups, with the two
         pre-specified contrasts marked
Panel B  gene-level effect sizes for the field-effect contrast, with
         bootstrap confidence intervals
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (MM, WIDTH_MM, C, STAGE_RAMP, blank_axes,  # noqa: E402
                      panel_label, save_all)

ORDER = ["Healthy liver", "NASH liver", "Cirrhotic liver",
         "Non-tumoral NASH liver adjacent to HCC", "NASH-HCC tumor"]
SHORT = {"Healthy liver": "Healthy", "NASH liver": "NASH\n(no HCC)",
         "Cirrhotic liver": "Cirrhotic",
         "Non-tumoral NASH liver adjacent to HCC": "Adjacent\nnon-tumor",
         "NASH-HCC tumor": "NASH-HCC\ntumor"}


def cliffs(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(((x[:, None] > y[None, :]).sum()
                  - (x[:, None] < y[None, :]).sum()) / (len(x) * len(y)))


def cliffs_ci(x, y, n_boot=4000, seed=0):
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x, float), np.asarray(y, float)
    d = [cliffs(rng.choice(x, len(x), replace=True),
                rng.choice(y, len(y), replace=True)) for _ in range(n_boot)]
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def fmt_p(p):
    return f"$P$ = {p:.3g}" if p >= 1e-4 else f"$P$ < 10$^{{-4}}$"


def bracket(ax, x0, x1, y, label, color):
    ax.plot([x0, x0, x1, x1], [y - 0.06, y, y, y - 0.06], color=color, lw=0.9,
            clip_on=False)
    ax.text((x0 + x1) / 2, y + 0.03, label, ha="center", va="bottom",
            fontsize=6.2, color=color, clip_on=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="GSE164760")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--figdir", default="results/figures")
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)

    d = os.path.join(args.outdir, args.cohort)
    scores = pd.read_csv(os.path.join(d, "rsi.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(d, "phenotype.tsv"), sep="\t",
                        index_col=0).loc[scores.index]
    tissue = pheno["tissue"].astype(str).str.strip()
    pc = json.load(open(os.path.join(d, "phase_c.json")))

    fig = plt.figure(figsize=(WIDTH_MM * MM, 145 * MM))

    # ============================================================= panel A
    A = fig.add_axes([0.085, 0.665, 0.885, 0.250])
    panel_label(blank_axes(fig, [0, 0.925, 1, 0.075]), 0.006, 0.80, "A")
    fig.text(0.052, 0.982, f"RSI across the HCC transition  ({args.cohort}, "
             f"n = {len(scores)}, HG-U219 array)", fontsize=8.2,
             fontweight="bold", va="top", color=C["ink"])

    data = [scores.loc[tissue == t, "RSI"].to_numpy(float) for t in ORDER]
    ramp = [STAGE_RAMP[i] for i in (0, 1, 2, 4, 5)]
    rng = np.random.default_rng(4)
    for i, (v, col) in enumerate(zip(data, ramp), start=1):
        A.boxplot([v], positions=[i], widths=0.55, patch_artist=True,
                  showfliers=False, zorder=2,
                  medianprops=dict(color=C["ink"], lw=1.3),
                  boxprops=dict(facecolor=col, edgecolor=C["ink"], lw=0.7),
                  whiskerprops=dict(color=C["ink"], lw=0.7),
                  capprops=dict(color=C["ink"], lw=0.7))
        A.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=4.5,
                  color=C["ink"], alpha=0.42, lw=0, zorder=3)

    A.set_xticks(range(1, 6))
    A.set_xticklabels([f"{SHORT[t]}\n$n$ = {len(v)}"
                       for t, v in zip(ORDER, data)], fontsize=6.4)
    A.set_ylabel("RSI", fontsize=7.4)
    A.tick_params(axis="y", labelsize=6.6, length=2.5, width=0.6)
    A.tick_params(axis="x", length=0)
    A.axhline(0, color=C["rule"], lw=0.6, ls=(0, (3, 3)), zorder=1)
    for sp in ("top", "right"):
        A.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        A.spines[sp].set_color(C["rule"])
    A.set_xlim(0.45, 5.55)
    top = max(v.max() for v in data)
    A.set_ylim(min(v.min() for v in data) - 0.15, top + 1.15)

    bracket(A, 4, 5, top + 0.22,
            f"H2  $\\delta$ = {pc['H2']['cliffs_delta']:+.2f},  "
            f"{fmt_p(pc['H2']['p_two_sided'])}   not supported", C["muted"])
    bracket(A, 2, 4, top + 0.62,
            f"H3  field effect   $\\delta$ = {pc['H3']['cliffs_delta']:+.2f},  "
            f"{fmt_p(pc['H3']['p_two_sided'])}", C["patho"])

    # ============================================================= panel B
    B = fig.add_axes([0.300, 0.130, 0.480, 0.335])
    panel_label(blank_axes(fig, [0, 0.500, 1, 0.075]), 0.006, 0.80, "B")
    fig.text(0.052, 0.552, "Gene-level effect sizes for the field-effect "
             "contrast", fontsize=8.2, fontweight="bold", va="top",
             color=C["ink"])

    adj = scores.loc[tissue == ORDER[3]]
    nash = scores.loc[tissue == ORDER[1]]
    rows = []
    for g, v in pc["genes"].items():
        col = f"gene_{g}"
        if col not in scores.columns:
            continue
        lo, hi = cliffs_ci(adj[col].to_numpy(float), nash[col].to_numpy(float))
        rows.append((g, v["H3"]["cliffs_delta"], lo, hi,
                     v["H3"].get("q_bh", 1.0)))
    rows.sort(key=lambda r: r[1])

    ys = np.arange(len(rows))
    for y, (g, delta, lo, hi, q) in zip(ys, rows):
        sig = q < 0.05
        col = C["protect"] if delta > 0 else C["patho"]
        col = col if sig else C["muted"]
        B.plot([lo, hi], [y, y], color=col, lw=1.0, zorder=2)
        B.plot([delta], [y], "o", ms=5.0 if sig else 3.6, color=col,
               mec="white", mew=0.8, zorder=3)
        B.text(-0.02, y, g, fontsize=6.4, ha="right", va="center",
               color=C["ink"] if sig else C["muted"],
               fontweight="bold" if sig else "normal",
               transform=B.get_yaxis_transform())
        if sig:
            B.text(1.02, y, f"$q$ = {q:.3f}", fontsize=6.1, ha="left",
                   va="center", color=col, transform=B.get_yaxis_transform())

    B.axvline(0, color=C["rule"], lw=0.8, ls=(0, (3, 3)), zorder=1)
    B.set_yticks([])
    B.set_ylim(-0.8, len(rows) - 0.2)
    B.set_xlabel("Cliff's $\\delta$   (adjacent non-tumor  vs  NASH without "
                 "HCC)", fontsize=7.0, labelpad=2)
    B.tick_params(axis="x", labelsize=6.5, length=2.5, width=0.6)
    for sp in ("top", "right", "left"):
        B.spines[sp].set_visible(False)
    B.spines["bottom"].set_color(C["rule"])

    fig.text(0.052, 0.014,
             "Bars are bootstrap 95% CIs. Filled markers pass BH-FDR within "
             "this contrast. Both pre-specified hypotheses were tested on the "
             "composite\nindex; the gene-level panel is descriptive. The array "
             "cohort is analyzed as an independent arm and never merged with "
             "the RNA-seq cohorts.",
             fontsize=6.0, color=C["muted"], va="bottom", linespacing=1.5)

    paths = save_all(fig, os.path.join(args.figdir,
                                       "Figure4_field_effect"))
    print("wrote:")
    for p in paths:
        print("  " + p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
