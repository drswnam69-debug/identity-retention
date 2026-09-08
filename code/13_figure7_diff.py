#!/usr/bin/env python3
"""13_figure7_diff.py -- Figure 7: the differentiation adjustment (PREREG 6g)."""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F

OUT = "results/figures"
d = json.load(open("results/GSE76427_differentiation.json"))


def pstr(p):
    if p >= 0.01:
        return f"P = {p:.3f}".rstrip("0").rstrip(".")
    e = int(np.floor(np.log10(p)))
    return f"$P = {p/10**e:.1f}\\times10^{{{e}}}$"


fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 86 * F.MM))

# ---- a: the composite, unadjusted -> D1 -> D2
axA = fig.add_axes([0.090, 0.245, 0.235, 0.650])
fig.text(0.010, 0.955, "a", fontsize=10, fontweight="bold")
levels = [("unadjusted", 1.129, None, None, F.C["muted"]),
          ("adjusted\nfor D1", d["paired_adjusted_D1"]["intercept"]["beta"],
           *d["paired_adjusted_D1"]["intercept"]["ci95"], F.C["ink"]),
          ("adjusted for D2\n(over-adjusted)",
           d["paired_adjusted_D2"]["intercept"]["beta"],
           *d["paired_adjusted_D2"]["intercept"]["ci95"], F.C["muted"])]
xs = np.arange(3)
for x, (lab, est, lo, hi, col) in zip(xs, levels):
    if lo is not None:
        axA.plot([x, x], [lo, hi], color=col, lw=1.2, solid_capstyle="round")
    axA.plot([x], [est], "o", ms=5.5, mfc=col, mec=col)
    axA.text(x + 0.20, est, f"{est:+.2f}", ha="left", va="center",
             fontsize=6.4, fontweight="bold", color=col)
axA.axhline(0, color=F.C["rule"], lw=0.7)
axA.axhline(1.129, color=F.C["muted"], lw=0.6, ls=(0, (3, 2)))
axA.set_xticks(xs)
axA.set_xticklabels([l[0] for l in levels], fontsize=6.2)
axA.set_ylabel("ΔRSI, tumor − adjacent\n(52 matched pairs)", fontsize=6.8)
axA.set_ylim(-0.15, 1.42)
axA.set_xlim(-0.55, 2.55)
axA.set_title("Half the rise survives adjustment", fontsize=7.0,
              fontweight="bold", pad=4)
for s_ in ("top", "right"):
    axA.spines[s_].set_visible(False)
axA.tick_params(labelsize=6.4, length=2.5, width=0.6)

# ---- b: module split
axB = fig.add_axes([0.430, 0.245, 0.235, 0.650])
fig.text(0.372, 0.955, "b", fontsize=10, fontweight="bold")
mods = [("REDUCTION", "reduction", F.C["protect"]),
        ("DRAIN", "drain", F.C["patho"])]
w = 0.32
for i, (lab, key, col) in enumerate(mods):
    m = d[f"module_{key}"]
    un = m["unadjusted_mean_delta"]
    ad = m["adjusted"]
    lo, hi = m["ci95"]
    axB.bar(i - w / 2, un, w, color=col, alpha=0.30, edgecolor=col, lw=0.8)
    axB.bar(i + w / 2, ad, w, color=col, edgecolor=col, lw=0.8)
    axB.plot([i + w / 2, i + w / 2], [lo, hi], color=F.C["ink"], lw=1.0)
    frac = abs(ad / un) * 100
    yy = max(un, ad, hi) + 0.09 if un > 0 else min(un, ad, lo) - 0.19
    axB.text(i, yy, f"{frac:.0f}% retained", ha="center", fontsize=6.3,
             fontweight="bold", color=col)
    axB.text(i, yy + (0.14 if un > 0 else -0.14), pstr(m["p"]),
             ha="center", fontsize=5.9, color=F.C["ink"])
axB.axhline(0, color=F.C["ink"], lw=0.8)
axB.set_xticks(range(len(mods)))
axB.set_xticklabels([m[0] for m in mods], fontsize=6.8)
axB.set_ylabel("Δ module score, tumor − adjacent", fontsize=6.8)
axB.set_ylim(-1.15, 1.30)
axB.set_title("The reduction arm is the durable one", fontsize=7.0,
              fontweight="bold", pad=4)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(facecolor=F.C["muted"], alpha=0.30,
                          edgecolor=F.C["muted"], label="unadjusted"),
                    Patch(facecolor=F.C["muted"], edgecolor=F.C["muted"],
                          label="D1-adjusted")],
           fontsize=5.9, frameon=False, loc="lower left", handlelength=1.1,
           handletextpad=0.4, borderpad=0.1)
for s_ in ("top", "right"):
    axB.spines[s_].set_visible(False)
axB.tick_params(labelsize=6.4, length=2.5, width=0.6)

# ---- c: what it means
axC = F.blank_axes(fig, [0.730, 0.215, 0.255, 0.680])
fig.text(0.700, 0.955, "c", fontsize=10, fontweight="bold")
axC.text(0.0, 1.00, "What the adjustment settles", fontsize=7.2,
         fontweight="bold", va="top")
lines = [
    (F.C["ink"], "The confound is real.", 
     "Hepatocyte identity (D1) is\nmarkedly lower in tumor,\n$P = 2.8\\times10^{-14}$."),
    (F.C["patho"], "It explains most of DRAIN.",
     "mARC1/mARC2/POR track\nhepatocyte identity — only\n41% of the fall survives."),
    (F.C["protect"], "It does not explain REDUCTION.",
     "71% of the rise survives and\nstays highly significant.\nCYB5R3-side signal is durable."),
]
y = 0.860
for col, head, body in lines:
    axC.plot([0.0, 0.028], [y + 0.012, y + 0.012], color=col, lw=2.2,
             solid_capstyle="round")
    axC.text(0.045, y + 0.012, head, fontsize=6.5, fontweight="bold",
             va="center", color=col)
    axC.text(0.045, y - 0.075, body, fontsize=5.9, va="top",
             color=F.C["muted"], linespacing=1.55)
    y -= 0.270
axC.plot([0.0, 1.0], [0.075, 0.075], color=F.C["rule"], lw=0.7)
axC.text(0.0, 0.030,
         "Composite retention is 50.3% — it clears\nthe pre-registered 50% threshold by only\n"
         "0.3 points, so that label is knife-edge.\nThe module split above is the real result.",
         fontsize=5.7, va="top", color=F.C["ink"], linespacing=1.6)

fig.text(0.090, 0.075,
         "GSE76427, 52 patient-matched pairs. D1 = 22 hepatocyte-identity genes "
         "(no redox enzyme, no P450);\nD2 = D1 + CYP2E1/3A4/1A2/2C9. Covariate and "
         "decision rule fixed in PREREGISTRATION §6g before extraction.",
         fontsize=5.5, color=F.C["muted"], va="top", linespacing=1.6)

paths = F.save_all(fig, os.path.join(OUT, "figure7_differentiation"))
print("wrote:", *paths, sep="\n  ")
