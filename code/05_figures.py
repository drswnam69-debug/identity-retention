#!/usr/bin/env python3
"""05_figures.py -- manuscript Figure 2: RSI across the MASLD spectrum.

Panel A  RSI by stage (boxes + points), with the pre-specified JT trend test
Panel B  which component moves the index
Panel C  cross-cohort forest of Spearman rho, with a random-effects pooled estimate

Example
-------
python3 code/05_figures.py --cohorts GSE135251 GSE130970 GSE167523 \
        --outdir results --figdir results/figures
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phenotype as ph        # noqa: E402
import stats_lite as sl       # noqa: E402


def cohort_axis(outdir: str, cohort: str):
    """The ordered severity ladder this cohort was prepared with."""
    with open(os.path.join(outdir, cohort, "qc.json")) as fh:
        qc = json.load(fh)
    return qc.get("stage_order", ph.STAGE_ORDER), qc.get("axis", "group")
from figstyle import (MM, WIDTH_MM, C, STAGE_RAMP, blank_axes,  # noqa: E402
                      panel_label, save_all)

SHORT = {"control": "Control", "NAFL": "NAFL", "NASH_F0-F1": "F0–1",
         "NASH_F2": "F2", "NASH_F3": "F3", "NASH_F4": "F4"}


def short(label: str) -> str:
    return SHORT.get(label, label)

SERIES = [("z_supply", "Supply", C["supply"]),
          ("z_reduction", "Reduction", C["protect"]),
          ("z_drain", "Drain", C["patho"])]


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "n/a"
    if p < 1e-4:
        e = int(math.floor(math.log10(p)))
        return f"$P$ = {p / 10 ** e:.1f} × 10$^{{{e}}}$"
    return f"$P$ = {p:.3g}"


def load(cohort: str, outdir: str):
    d = os.path.join(outdir, cohort)
    scores = pd.read_csv(os.path.join(d, "rsi.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(d, "phenotype.tsv"), sep="\t",
                        index_col=0).loc[scores.index]
    stage = pheno["stage"].astype(str)
    keep = stage.isin(ph.STAGE_ORDER)
    with open(os.path.join(d, "phase_b.json")) as fh:
        stats = json.load(fh)
    return scores[keep], stage[keep], stats


def pooled_random_effects(rhos, ns):
    """DerSimonian-Laird pooling on Fisher z. Returns (rho, lo, hi, i2)."""
    z = np.array([math.atanh(r) for r in rhos])
    v = np.array([(1.06 ** 2) / (n - 3) for n in ns])
    w = 1 / v
    zf = float((w * z).sum() / w.sum())
    q = float((w * (z - zf) ** 2).sum())
    df = len(z) - 1
    c = float(w.sum() - (w ** 2).sum() / w.sum())
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1 / (v + tau2)
    zr = float((wr * z).sum() / wr.sum())
    se = math.sqrt(1 / wr.sum())
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    return math.tanh(zr), math.tanh(zr - 1.96 * se), math.tanh(zr + 1.96 * se), i2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="+", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--figdir", default="results/figures")
    ap.add_argument("--discovery", default=None,
                    help="cohort shown in panels A and B (default: first)")
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)
    disc = args.discovery or args.cohorts[0]

    scores, stage, stats = load(disc, args.outdir)
    order, axis = cohort_axis(args.outdir, disc)
    fig = plt.figure(figsize=(WIDTH_MM * MM, 150 * MM))

    # =========================================================== panel A
    A = fig.add_axes([0.085, 0.668, 0.860, 0.272])
    panel_label(fig.add_axes([0, 0.955, 1, 0.045], frameon=False,
                             xticks=[], yticks=[]), 0.006, 0.9, "A")
    fig.text(0.052, 0.978, f"RSI across the MASLD spectrum  ({disc}, "
             f"n = {len(scores)})", fontsize=8.2, fontweight="bold",
             va="top", color=C["ink"])

    present = [s for s in order if (stage == s).any()]
    ramp = [STAGE_RAMP[min(len(STAGE_RAMP) - 1,
                           round(i * (len(STAGE_RAMP) - 1)
                                 / max(1, len(present) - 1)))]
            for i in range(len(present))]
    data = [scores.loc[stage == s, "RSI"].to_numpy(float) for s in present]
    rng = np.random.default_rng(0)

    for i, (vals, col) in enumerate(zip(data, ramp), start=1):
        bp = A.boxplot([vals], positions=[i], widths=0.56, patch_artist=True,
                       showfliers=False, zorder=2,
                       medianprops=dict(color=C["ink"], lw=1.3),
                       boxprops=dict(facecolor=col, edgecolor=C["ink"], lw=0.7),
                       whiskerprops=dict(color=C["ink"], lw=0.7),
                       capprops=dict(color=C["ink"], lw=0.7))
        x = i + rng.uniform(-0.16, 0.16, len(vals))
        A.scatter(x, vals, s=4.0, color=C["ink"], alpha=0.40, lw=0, zorder=3)
        A.text(i, A.get_ylim()[0], "", fontsize=6)

    A.set_xticks(range(1, len(present) + 1))
    A.set_xticklabels([f"{short(s)}\n$n$ = {len(v)}"
                       for s, v in zip(present, data)], fontsize=6.6)
    A.set_ylabel("RSI", fontsize=7.4)
    A.tick_params(axis="y", labelsize=6.6, length=2.5, width=0.6)
    A.tick_params(axis="x", length=0)
    A.axhline(0, color=C["rule"], lw=0.6, ls=(0, (3, 3)), zorder=1)
    for sp in ("top", "right"):
        A.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        A.spines[sp].set_color(C["rule"])
    A.set_xlim(0.4, len(present) + 0.6)

    pr = stats["primary_RSI"]
    A.text(0.015, 0.965,
           "Jonckheere–Terpstra, increasing:  "
           f"$z$ = {pr['z']:.2f},  {fmt_p(pr['p_increasing'])}\n"
           f"Spearman $\\rho$ = {pr['spearman_rho']:.2f} "
           f"(95% CI {pr['rho_ci95'][0]:.2f} to {pr['rho_ci95'][1]:.2f})",
           transform=A.transAxes, fontsize=6.5, ha="left", va="top",
           color=C["ink"], linespacing=1.5)

    # =========================================================== panel B
    B = fig.add_axes([0.085, 0.383, 0.860, 0.200])
    panel_label(fig.add_axes([0, 0.585, 1, 0.045], frameon=False,
                             xticks=[], yticks=[]), 0.006, 0.9, "B")
    fig.text(0.052, 0.612, "Which component moves the index",
             fontsize=8.2, fontweight="bold", va="top", color=C["ink"])

    xs = np.arange(1, len(present) + 1)
    for key, label, col in SERIES:
        med = np.array([np.median(scores.loc[stage == s, key]) for s in present])
        lo = np.array([np.percentile(scores.loc[stage == s, key], 25)
                       for s in present])
        hi = np.array([np.percentile(scores.loc[stage == s, key], 75)
                       for s in present])
        B.fill_between(xs, lo, hi, color=col, alpha=0.13, lw=0, zorder=1)
        B.plot(xs, med, color=col, lw=2.0, zorder=3, solid_capstyle="round")
        B.plot(xs, med, "o", ms=3.6, color=col, mec="white", mew=0.8, zorder=4)
        B.plot([xs[-1] + 0.10, xs[-1] + 0.22], [med[-1], med[-1]], color=col,
               lw=2.0, zorder=3)
        B.text(xs[-1] + 0.27, med[-1], label, fontsize=6.6, va="center",
               ha="left", color=C["ink"])

    B.axhline(0, color=C["rule"], lw=0.6, ls=(0, (3, 3)), zorder=0)
    B.set_xticks(xs)
    B.set_xticklabels([short(s) for s in present], fontsize=6.6)
    B.set_ylabel("module $z$-score\n(median, IQR)", fontsize=7.2,
                 linespacing=1.4)
    B.tick_params(axis="y", labelsize=6.6, length=2.5, width=0.6)
    B.tick_params(axis="x", length=0)
    for sp in ("top", "right"):
        B.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        B.spines[sp].set_color(C["rule"])
    B.set_xlim(0.6, len(present) + 1.45)

    handles = [plt.Line2D([], [], color=col, lw=2.0, label=lab)
               for _, lab, col in SERIES]
    leg = B.legend(handles=handles, loc="lower left", fontsize=6.3,
                   frameon=False, ncol=3, handlelength=1.5,
                   columnspacing=1.2, borderpad=0.1,
                   bbox_to_anchor=(0.0, -0.02))
    for t in leg.get_texts():
        t.set_color(C["ink"])

    cs = stats["components"]
    B.text(0.985, 0.04, "  ".join(
        f"{lab}: $z$ = {cs[k]['z']:.1f}" for k, lab, _ in SERIES),
        transform=B.transAxes, fontsize=6.2, ha="right", va="bottom",
        color=C["muted"])

    # =========================================================== panel C
    rows = []
    for c in args.cohorts:
        try:
            with open(os.path.join(args.outdir, c, "phase_b.json")) as fh:
                st = json.load(fh)
        except FileNotFoundError:
            continue
        p = st["primary_RSI"]
        rows.append((c, p["spearman_rho"], p["rho_ci95"][0], p["rho_ci95"][1],
                     st["n"]))
    if not rows:
        raise SystemExit("no cohort results found")

    # panel C is only as tall as it has rows to show
    n_lines = len(rows) + (1 if len(rows) > 1 else 0)
    c_h = min(0.150, 0.028 * n_lines + 0.030)
    c_bottom = 0.128 + (0.150 - c_h)
    Cx = fig.add_axes([0.300, c_bottom, 0.500, c_h])
    panel_label(fig.add_axes([0, 0.290, 1, 0.045], frameon=False,
                             xticks=[], yticks=[]), 0.006, 0.9, "C")
    fig.text(0.052, 0.318, "Replication across cohorts",
             fontsize=8.2, fontweight="bold", va="top", color=C["ink"])

    ys = list(range(len(rows), 0, -1))
    for y, (c, rho, lo, hi, n) in zip(ys, rows):
        Cx.plot([lo, hi], [y, y], color=C["ink"], lw=1.0, zorder=2)
        Cx.plot([rho], [y], "s", ms=5.0, color=C["supply"], mec="white",
                mew=0.8, zorder=3)
        Cx.text(-0.02, y, f"{c}", fontsize=6.5, ha="right", va="center",
                color=C["ink"], transform=Cx.get_yaxis_transform())
        Cx.text(1.02, y, f"{rho:.2f} ({lo:.2f}, {hi:.2f})   $n$ = {n}",
                fontsize=6.3, ha="left", va="center", color=C["muted"],
                transform=Cx.get_yaxis_transform())

    if len(rows) > 1:
        prho, plo, phi, i2 = pooled_random_effects([r[1] for r in rows],
                                                   [r[4] for r in rows])
        y0 = 0.15
        Cx.plot([plo, phi], [y0, y0], color=C["ink"], lw=1.4, zorder=2)
        Cx.plot([prho], [y0], "D", ms=6.0, color=C["ink"], zorder=3)
        Cx.text(-0.02, y0, "Pooled", fontsize=6.5, fontweight="bold",
                ha="right", va="center", color=C["ink"],
                transform=Cx.get_yaxis_transform())
        Cx.text(1.02, y0, f"{prho:.2f} ({plo:.2f}, {phi:.2f})   "
                f"$I^2$ = {i2:.0f}%", fontsize=6.3, ha="left", va="center",
                color=C["ink"], transform=Cx.get_yaxis_transform())
        ys = ys + [y0]

    Cx.axvline(0, color=C["rule"], lw=0.8, ls=(0, (3, 3)), zorder=1)
    Cx.set_ylim(min(ys) - 0.6, max(ys) + 0.6)
    Cx.set_yticks([])
    Cx.set_xlabel("Spearman $\\rho$ (RSI vs. each cohort's severity axis)",
                  fontsize=7.0,
                  labelpad=2)
    Cx.tick_params(axis="x", labelsize=6.5, length=2.5, width=0.6)
    for sp in ("top", "right", "left"):
        Cx.spines[sp].set_visible(False)
    Cx.spines["bottom"].set_color(C["rule"])

    note = ("Higher RSI = greater reductive-supply dominance. The index was "
            "locked before any outcome variable was opened; lock hash "
            f"{stats['lock_hash'][:12]}.")
    if len(rows) == 1:
        note = ("Single cohort shown; pooling requires two or more. " + note)
    import textwrap
    fig.text(0.052, 0.018, "\n".join(textwrap.wrap(note, 118)), fontsize=6.0,
             color=C["muted"], va="bottom", linespacing=1.5)

    paths = save_all(fig, os.path.join(args.figdir, "Figure2_RSI_spectrum"))
    print("wrote:")
    for p in paths:
        print("  " + p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
