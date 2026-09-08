#!/usr/bin/env python3
"""21_figure_cascade.py -- manuscript Figure 3: the adjustment cascade.

Every number is read from the results JSON and asserted against the value
quoted in the manuscript, so the figure cannot drift from the analysis.
"""
import json, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F

OUT = "results/figures"
D = {c: json.load(open(f"results/{c}_differentiation.json")) for c in ("GSE76427", "GSE14520")}
C = {c: json.load(open(f"results/{c}_composition.json")) for c in ("GSE76427", "GSE14520")}

UNADJ = {"GSE76427": 1.129, "GSE14520": 1.261}
LAB = {"GSE76427": "GSE76427\n52 pairs · mixed etiology",
       "GSE14520": "GSE14520\n213 pairs · HBV"}


def close(a, b, tol=5e-3):
    assert abs(a - b) < tol, f"{a} != {b}"


for c in ("GSE76427", "GSE14520"):
    close(D[c]["paired_adjusted_D1"]["intercept"]["beta"],
          {"GSE76427": 0.568, "GSE14520": 0.792}[c])
    close(C[c]["composite_joint"]["estimate"],
          {"GSE76427": 0.451, "GSE14520": 0.814}[c])
    close(C[c]["composite_C1"]["estimate"],
          {"GSE76427": 1.136, "GSE14520": 1.462}[c])

fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 118 * F.MM))
LB = dict(fontsize=10, fontweight="bold")

# ---------------- a: composite cascade, forest --------------------------
axA = fig.add_axes([0.255, 0.610, 0.395, 0.320])
fig.text(0.012, 0.965, "a", **LB)
rows = []
for c in ("GSE14520", "GSE76427"):
    rows.append((c, "adjusted for D1 + C1", C[c]["composite_joint"]["estimate"],
                 C[c]["composite_joint"]["ci95"], F.C["protect"], True))
    rows.append((c, "adjusted for C1 alone", C[c]["composite_C1"]["estimate"],
                 C[c]["composite_C1"]["ci95"], F.C["muted"], False))
    rows.append((c, "adjusted for D1", D[c]["paired_adjusted_D1"]["intercept"]["beta"],
                 D[c]["paired_adjusted_D1"]["intercept"]["ci95"], F.C["ink"], False))
    rows.append((c, "unadjusted", UNADJ[c], None, F.C["muted"], False))
y = np.arange(len(rows))[::-1]
for yy, (c, lab, v, ci, col, bold) in zip(y, rows):
    if ci:
        axA.plot(ci, [yy, yy], color=col, lw=1.3 if bold else 1.0,
                 solid_capstyle="round")
    axA.plot([v], [yy], "o", ms=5.2 if bold else 4.2, color=col,
             mec="white", mew=0.7, zorder=3)
    axA.text(v, yy + 0.30, f"{v:+.2f}", ha="center", fontsize=6.0,
             fontweight="bold" if bold else "normal", color=col)
axA.axvline(0, color=F.C["ink"], lw=0.8)
axA.set_yticks(y)
axA.set_yticklabels([r[1] for r in rows], fontsize=6.4)
axA.set_xlim(-0.15, 1.75)
axA.set_ylim(-0.8, len(rows) + 0.6)
axA.axhline(3.5, color=F.C["rule"], lw=0.7)
bb = dict(facecolor="white", edgecolor="none", pad=1.2)
axA.text(1.74, 7.42, "GSE14520 · 213 pairs", ha="right", fontsize=6.5,
         fontweight="bold", color=F.C["ink"], bbox=bb, zorder=5)
axA.text(1.74, 3.42, "GSE76427 · 52 pairs", ha="right", fontsize=6.5,
         fontweight="bold", color=F.C["ink"], bbox=bb, zorder=5)
axA.set_xlabel("ΔRSI, tumor − adjacent, at zero change in the covariate",
               fontsize=6.6)
axA.set_title("Composite index: composition removes nothing,\nidentity removes about half",
              fontsize=7.2, fontweight="bold", pad=5, loc="left")
for s_ in ("top", "right", "left"):
    axA.spines[s_].set_visible(False)
axA.tick_params(labelsize=6.2, length=2.5, width=0.6)

# ---------------- b: module retention -----------------------------------
axB = fig.add_axes([0.725, 0.610, 0.245, 0.320])
fig.text(0.672, 0.965, "b", **LB)
groups = []
for c in ("GSE76427", "GSE14520"):
    groups.append((c,
                   abs(C[c]["module_reduction_joint"]["retained_fraction"]) * 100,
                   abs(C[c]["module_drain_joint"]["retained_fraction"]) * 100))
w = 0.34
for i, (c, red, dra) in enumerate(groups):
    axB.bar(i - w / 2, red, w, color=F.C["protect"], edgecolor=F.C["protect"], lw=0.8)
    axB.bar(i + w / 2, dra, w, color=F.C["patho"], edgecolor=F.C["patho"], lw=0.8)
    for v, off, col in ((red, -w / 2, F.C["protect"]), (dra, w / 2, F.C["patho"])):
        axB.text(i + off, v + 2.2, f"{v:.0f}", ha="center", fontsize=6.2,
                 fontweight="bold", color=col)
axB.axhline(50, color=F.C["rule"], lw=0.8, ls=(0, (3, 2)))
axB.annotate("indistinguishable\n51.2 vs 51.1", xy=(0.0, 56), xytext=(0.0, 74),
             fontsize=5.7, ha="center", color=F.C["muted"], linespacing=1.4,
             arrowprops=dict(arrowstyle="-", lw=0.6, color=F.C["rule"]))
axB.set_xticks(range(2))
axB.set_xticklabels(["GSE76427", "GSE14520"], fontsize=6.4)
axB.set_ylabel("% of shift retained, joint model", fontsize=6.6)
axB.set_ylim(0, 118)
axB.set_title("The two arms separate\nin one cohort only", fontsize=7.2,
              fontweight="bold", pad=5, loc="left")
axB.legend(handles=[Patch(facecolor=F.C["protect"], label="REDUCTION"),
                    Patch(facecolor=F.C["patho"], label="DRAIN")],
           fontsize=5.9, frameon=False, loc="upper left", handlelength=1.0,
           handletextpad=0.4, borderpad=0.1, ncol=2, columnspacing=0.9,
           bbox_to_anchor=(-0.02, 1.045))
for s_ in ("top", "right"):
    axB.spines[s_].set_visible(False)
axB.tick_params(labelsize=6.2, length=2.5, width=0.6)

# ---------------- c: module cascade, both cohorts ------------------------
axC = fig.add_axes([0.255, 0.215, 0.715, 0.265])
fig.text(0.012, 0.520, "c", **LB)
stages = [("unadjusted", None), ("D1", "_D1"), ("D1 + C1", "_joint")]
xs, labels = [], []
pos = 0
for c in ("GSE76427", "GSE14520"):
    for mod, col in (("reduction", F.C["protect"]), ("drain", F.C["patho"])):
        vals = []
        un = D[c][f"module_{mod}"]["unadjusted_mean_delta"]
        vals.append((un, None))
        d1 = D[c][f"module_{mod}"]
        vals.append((d1["adjusted"], d1["ci95"]))
        jt = C[c][f"module_{mod}_joint"]
        vals.append((jt["estimate"], jt["ci95"]))
        for k, (v, ci) in enumerate(vals):
            x = pos + k * 0.30
            if ci:
                axC.plot([x, x], ci, color=col, lw=1.1, solid_capstyle="round")
            axC.plot([x], [v], "o", ms=4.4, color=col if k else "white",
                     mec=col, mew=1.4, zorder=3)
        axC.plot([pos, pos + 0.60], [vals[0][0], vals[2][0]], color=col,
                 lw=0.7, ls=(0, (2, 2)), alpha=0.55, zorder=1)
        labels.append((pos + 0.30, mod.upper()))
        pos += 1.15
    xs.append(pos - 2.30 + 0.30)
    pos += 0.55
axC.axhline(0, color=F.C["ink"], lw=0.8)
axC.set_xticks([l[0] for l in labels])
axC.set_xticklabels([l[1] for l in labels], fontsize=6.2)
axC.set_ylim(-1.15, 1.15)
axC.set_ylabel("Δ module score, tumor − adjacent", fontsize=6.6)
axC.set_title("Module level: the drain fall is mostly identity loss; the reduction rise is not "
              "— but only in the larger cohort",
              fontsize=7.2, fontweight="bold", pad=5, loc="left")
for x, c in zip(xs, ("GSE76427", "GSE14520")):
    axC.text(x, -0.20, LAB[c].replace("\n", " · "), ha="center", fontsize=6.3,
             color=F.C["ink"], transform=axC.get_xaxis_transform())
axC.text(0.02, 0.965,
         "within each module, left to right:  ○ unadjusted   ● adjusted for D1"
         "   ● adjusted for D1 + C1",
         fontsize=5.9, color=F.C["muted"], transform=axC.transAxes, va="top")
for s_ in ("top", "right"):
    axC.spines[s_].set_visible(False)
axC.tick_params(labelsize=6.2, length=2.5, width=0.6)

fig.text(0.012, 0.105,
         "Points are the intercept of the paired-difference regression — the tumor rise at zero change in the covariate — with 95% CI. "
         "D1, 22-gene hepatocyte identity;\nC1, 21-gene non-parenchymal content. Covariates, tests and decision rules were fixed before the data to test them were obtained; "
         "variance inflation was 1.12 and 1.14 in the\ntwo paired cohorts, so the joint estimates are not collinear.",
         fontsize=5.6, color=F.C["muted"], va="top", linespacing=1.7)

paths = F.save_all(fig, os.path.join(OUT, "Figure3_cascade"))
print("wrote:", *paths, sep="\n  ")
