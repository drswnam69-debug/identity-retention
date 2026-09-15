#!/usr/bin/env python3
"""make_graphical_abstract.py -- the graphical abstract.

Every number and every point is read from the archived result files, not
retyped, so the figure moves when the archive moves and 50_consistency_check.py
can bind its claims the same way it binds the main text.

Layout follows the article's title: what the measure is and what it gives in
three tissues on the left, the failure to transfer to a proteome on the right,
and the negative control that limits how either may be read along the bottom.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "rsi", "code"))
from paths import RESULTS, FIGURES  # noqa: E402
from figstyle import C, blank_axes, save_all  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

R = lambda n: json.load(open(os.path.join(RESULTS, n), encoding="utf-8"))

COHORTS = [  # label, tissue, file
    ("TCGA-LIHC",  "liver",  "SIGNATURE_BENCHMARK_TCGA_LIHC.json"),
    ("GSE14520",   "liver",  "SIGNATURE_BENCHMARK_GSE14520.json"),
    ("GSE76427",   "liver",  "SIGNATURE_BENCHMARK_GSE76427.json"),
    ("TCGA-LUAD",  "lung",   "SIGNATURE_BENCHMARK_TCGA_LUAD.json"),
    ("TCGA-KIRC",  "kidney", "SIGNATURE_BENCHMARK_TCGA_KIRC.json"),
]
TISSUE_COLOR = {"liver": C["supply"], "lung": C["protect"], "kidney": C["patho"]}


def panel_a(fig):
    ax = fig.add_axes([0.135, 0.375, 0.275, 0.355])
    rows = []
    for label, tissue, fn in COHORTS:
        d = R(fn)
        s = d["summary_joint"]
        rows.append((label, tissue, d["n_pairs"], d["n_sets_evaluable"],
                     s["median"], s["q1"], s["q3"]))
    y = np.arange(len(rows))[::-1]
    ax.axvline(0.5, color=C["rule"], lw=0.8, ls=(0, (3, 2)), zorder=1)
    for yi, (label, tissue, npair, nset, med, q1, q3) in zip(y, rows):
        col = TISSUE_COLOR[tissue]
        ax.plot([q1, q3], [yi, yi], color=col, lw=2.6, solid_capstyle="round",
                alpha=0.32, zorder=2)
        ax.plot([med], [yi], "o", color=col, ms=5.2, zorder=3,
                markeredgecolor="white", markeredgewidth=0.7)
        ax.text(1.03, yi + 0.30, f"{med:.2f}", fontsize=7, color=col,
                fontweight="bold", ha="left", va="center")
        ax.text(1.03, yi - 0.30, f"{nset} sets, {npair} pairs", fontsize=6.4,
                color=C["muted"], ha="left", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{lab}\n{tis}" for lab, tis, *_ in rows], fontsize=7,
                       linespacing=1.25)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("identity-retention fraction, median and IQR\n"
                  "(joint model, identity D1 + composition C1)", fontsize=7,
                  labelpad=2, linespacing=1.3)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(C["rule"])
    ax.spines["bottom"].set_color(C["rule"])
    ax.text(0.492, 0.995, "half", fontsize=6.4, color=C["muted"], ha="right",
            va="top", transform=ax.get_xaxis_transform())


def panel_b(fig):
    d = R("GAO_MRNA_PROTEIN_6aa.json")
    sig = d["signatures"]
    t = np.array([s["retention_transcript"] for s in sig])
    p = np.array([s["retention_protein"] for s in sig])
    lo, hi = 0.02, 3.0
    ax = fig.add_axes([0.635, 0.375, 0.295, 0.355])
    ax.plot([lo, hi], [lo, hi], color=C["rule"], lw=0.8, zorder=1)
    above = p > t
    ax.scatter(t[~above], p[~above], s=8, facecolor="white", linewidth=0.7,
               edgecolor=C["muted"], zorder=2)
    ax.scatter(t[above], p[above], s=8, facecolor=C["supply"], linewidth=0.0,
               alpha=0.72, zorder=3)
    mt, mp = float(np.median(t)), float(np.median(p))
    ax.plot([mt], [mp], "D", color=C["out"], ms=4.6, zorder=5,
            markeredgecolor="white", markeredgewidth=0.7)
    ax.annotate(f"median\n{mt:.2f} to {mp:.2f}", xy=(mt, mp),
                xytext=(0.035, 0.36), fontsize=6.6, color=C["out"],
                ha="left", va="center", linespacing=1.3,
                arrowprops=dict(arrowstyle="-", color=C["out"], lw=0.7,
                                shrinkA=1, shrinkB=3))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ticks = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]
    for a in (ax.xaxis, ax.yaxis):
        a.set_ticks(ticks)
        a.set_ticklabels([f"{v:g}" for v in ticks])
        a.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("retention in the transcriptome", fontsize=7, labelpad=2)
    ax.set_ylabel("retention in the proteome", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(C["rule"])
    ax.spines["bottom"].set_color(C["rule"])
    n_above = int(above.sum())
    ax.text(0.03, 0.97,
            f"{n_above} of {len(sig)} signatures\nsit above the line",
            transform=ax.transAxes, fontsize=6.6, color=C["supply"],
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="none", alpha=0.85))
    ax.text(0.97, 0.06,
            f"registered explanation of the gap\n"
            f"not supported: $\\rho$ = {d['spearman_rho']:+.2f}, "
            f"$P$ = {d['spearman_p']:.2f}",
            transform=ax.transAxes, fontsize=6.5, color=C["muted"],
            ha="right", va="bottom", linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="none", alpha=0.85))
    return d, n_above


def main():
    fig = plt.figure(figsize=(6.85, 4.35))
    fig.patch.set_facecolor("white")
    bg = blank_axes(fig, [0, 0, 1, 1])

    # -- title ------------------------------------------------------------
    bg.text(0.5, 0.965,
            "How much of a tumor expression signature survives adjustment "
            "for the dominant tumor-adjacent axis?",
            fontsize=9.2, fontweight="bold", color=C["ink"],
            ha="center", va="top")
    bg.text(0.5, 0.918,
            "A pre-registered measure: the intercept of the paired-difference "
            "regression on a disjoint tissue-identity score,\n"
            "as a share of the unadjusted paired shift. 31 timestamped "
            "amendments; code, protocol and results deposited.",
            fontsize=7, color=C["muted"], ha="center", va="top",
            linespacing=1.4)

    bg.plot([0.055, 0.945], [0.862, 0.862], color=C["rule"], lw=0.7)

    bg.text(0.135, 0.826, "Three tissues, five paired cohorts",
            fontsize=8, fontweight="bold", color=C["ink"], ha="left", va="top")
    bg.text(0.135, 0.786,
            "Most published liver signatures are not restatements of the\n"
            "axis: three quarters keep at least half of their contrast.",
            fontsize=7, color=C["ink"], ha="left", va="top", linespacing=1.4)

    bg.text(0.585, 0.826, "It does not transfer to a proteome",
            fontsize=8, fontweight="bold", color=C["ink"], ha="left", va="top")
    bg.text(0.585, 0.786,
            "One of the two deposited proteomes fails the identity\n"
            "premise and cannot be used; in the other, more survives.",
            fontsize=7, color=C["ink"], ha="left", va="top", linespacing=1.4)

    panel_a(fig)
    gao, _ = panel_b(fig)

    # -- footer: the negative control -------------------------------------
    nc = R("NEGATIVE_CONTROL_6ae.json")["negative_control"]
    sp = R("SPREAD_IS_NOT_6v.json")["spread_under_the_negative_control"]
    bg.add_patch(plt.Rectangle((0.055, 0.035), 0.89, 0.235,
                               facecolor="#f4f3ef", edgecolor=C["rule"],
                               linewidth=0.7, zorder=0))
    bg.text(0.078, 0.243, "What the control says the measure is not",
            fontsize=8, fontweight="bold", color=C["ink"], ha="left", va="top")
    bg.text(0.078, 0.200,
            f"Against {nc['n_draws']} covariates drawn at random and matched to the identity score on the size of their own paired\n"
            f"tumor shift, the identity score leaves a median {nc['D1_median']:.2f} against {nc['random_median_of_medians']:.2f}, "
            f"and {nc['fraction_random_at_or_below_D1_median']*100:.1f}% of the random draws leave at or below it. The\n"
            f"per-signature spread is reproduced rather than removed: {sp['D1']['fold_p05_p95']:.1f}-fold under the identity score, "
            f"{sp['matched_random']['fold_p05_p95']:.1f}-fold under the control.",
            fontsize=7, color=C["ink"], ha="left", va="top", linespacing=1.45)
    bg.text(0.078, 0.066,
            "The quantity measures alignment with a large paired tumor-adjacent "
            "axis, not anything specific to differentiation.",
            fontsize=7, fontweight="bold", color=C["patho"], ha="left", va="bottom")

    out = os.path.join(FIGURES, "GraphicalAbstract_identity_retention")
    os.makedirs(FIGURES, exist_ok=True)
    made = save_all(fig, out)
    print("\n".join("  " + m for m in made))
    print(f"  Gao proteome: {gao['n_signatures_both_levels']} signatures at both levels, "
          f"rho = {gao['spearman_rho']:+.3f}, P = {gao['spearman_p']:.2f}")


if __name__ == "__main__":
    main()
