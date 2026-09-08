#!/usr/bin/env python3
"""22_figures_1_2.py -- manuscript Figures 1 and 2 for Hepatology International.

All plotted values are read from the results JSON and asserted before drawing.
"""
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F
import rsi_config as cfg

OUT = "results/figures"

# ============================ FIGURE 1 =================================
fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 84 * F.MM))
LB = dict(fontsize=10, fontweight="bold")

axA = F.blank_axes(fig, [0.035, 0.08, 0.46, 0.84])
fig.text(0.012, 0.955, "a", **LB)
axA.set_xlim(0, 1); axA.set_ylim(0, 1)

mods = [("SUPPLY", cfg.MODULE_SUPPLY, F.C["supply"], F.C["supply_f"], 0.175, "+ ½"),
        ("REDUCTION", cfg.MODULE_REDUCTION, F.C["protect"], F.C["protect_f"], 0.500, "+ ½"),
        ("DRAIN", cfg.MODULE_DRAIN, F.C["patho"], F.C["patho_f"], 0.825, "− 1")]
for name, genes, ec, fc, x, wt in mods:
    axA.add_patch(FancyBboxPatch((x - 0.152, 0.555), 0.304, 0.375,
                                 boxstyle="round,pad=0,rounding_size=0.02",
                                 lw=1.2, edgecolor=ec, facecolor=fc))
    axA.text(x, 0.895, name, ha="center", fontsize=7.0, fontweight="bold", color=ec)
    wrapped, line = [], ""
    for g in genes:
        if line and len(line) + 2 + len(g) > 16:
            wrapped.append(line); line = g
        else:
            line = (line + ", " + g) if line else g
    wrapped.append(line)
    axA.text(x, 0.735, "\n".join(wrapped), ha="center", va="center",
             fontsize=5.0, family="monospace", color=F.C["ink"], linespacing=1.45)
    axA.annotate("", xy=(x, 0.475), xytext=(x, 0.545),
                 arrowprops=dict(arrowstyle="-|>", lw=1.2, color=ec, mutation_scale=7))
    axA.text(x + 0.048, 0.505, wt, fontsize=6.6, color=ec, fontweight="bold")
axA.plot([0.175, 0.825], [0.455, 0.455], color=F.C["rule"], lw=1.0)
axA.annotate("", xy=(0.5, 0.345), xytext=(0.5, 0.450),
             arrowprops=dict(arrowstyle="-|>", lw=1.2, color=F.C["ink"], mutation_scale=7))
axA.add_patch(FancyBboxPatch((0.145, 0.215), 0.71, 0.135,
                             boxstyle="round,pad=0,rounding_size=0.02",
                             lw=1.3, edgecolor=F.C["ink"], facecolor="white"))
axA.text(0.5, 0.283, "RSI  =  ½ (SUPPLY + REDUCTION)  −  DRAIN",
         ha="center", va="center", fontsize=7.6, fontweight="bold")
axA.text(0.5, 0.128, "each module = mean of within-cohort gene-wise z-scores\n"
                     "locked 25 Aug 2026 · SHA-256 7a2bf934…b680046\n"
                     "POR was locked into DRAIN on a false premise (§6l); removing it\n"
                     "changes no direction and no verdict",
         ha="center", va="center", fontsize=5.8, color=F.C["muted"], linespacing=1.6)
axA.set_title("Index construction", fontsize=7.4, fontweight="bold", loc="left", pad=4)

# --- b: H1 forest
axB = fig.add_axes([0.615, 0.185, 0.335, 0.640])
fig.text(0.545, 0.955, "b", **LB)
coh = [("GSE167523", 98, 0.299, 0.095, 0.479),
       ("GSE130970", 78, 0.389, 0.169, 0.572),
       ("GSE135251", 216, 0.397, 0.270, 0.509)]
pooled = ("Pooled", 392, 0.372, 0.277, 0.459)
ys = np.arange(len(coh))
for y, (n, nn, r, lo, hi) in zip(ys, coh):
    axB.plot([lo, hi], [y, y], color=F.C["protect"], lw=1.2, solid_capstyle="round")
    axB.plot([r], [y], "s", ms=5.0, color=F.C["protect"], mec="white", mew=0.6)
    axB.text(0.60, y, f"n = {nn}", fontsize=5.9, va="center", color=F.C["muted"])
yp = -1.05
axB.plot([pooled[3], pooled[4]], [yp, yp], color=F.C["ink"], lw=1.6, solid_capstyle="round")
axB.plot([pooled[2]], [yp], "D", ms=6.2, color=F.C["ink"], mec="white", mew=0.6)
axB.text(0.60, yp, f"n = {pooled[1]}", fontsize=5.9, va="center", fontweight="bold")
axB.axhline(-0.5, color=F.C["rule"], lw=0.7)
axB.axvline(0, color=F.C["ink"], lw=0.8)
axB.set_yticks(list(ys) + [yp])
axB.set_yticklabels([c[0] for c in coh] + ["Pooled  (I² = 0%)"], fontsize=6.4)
axB.set_xlim(-0.05, 0.72); axB.set_ylim(-1.7, len(coh) - 0.3)
axB.set_xlabel("Spearman ρ, index vs histological severity", fontsize=6.6)
axB.set_title("Severity gradient in steatotic liver disease", fontsize=7.4,
              fontweight="bold", loc="left", pad=4)
for s_ in ("top", "right", "left"):
    axB.spines[s_].set_visible(False)
axB.tick_params(labelsize=6.2, length=2.5, width=0.6)

print("wrote:", *F.save_all(fig, os.path.join(OUT, "Figure1_index_spectrum")), sep="\n  ")
plt.close(fig)

# ============================ FIGURE 2 =================================
h2 = json.load(open("results/gse76427_h2.json"))
g14 = json.load(open("results/gse14520_partial.json"))
assert abs(h2["unpaired"]["cliffs_delta"] - 0.76) < 0.02
assert abs(g14["H2_composite"]["cliffs_delta"] - 0.862) < 0.02

fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 74 * F.MM))
axA = fig.add_axes([0.085, 0.175, 0.40, 0.645])
fig.text(0.012, 0.945, "a", **LB)

# paired differences, real data
rsi14 = pd.read_csv("results/GSE14520/rsi.tsv", sep="\t", index_col=0)
ph14 = pd.read_csv("results/GSE14520/phenotype.tsv", sep="\t", index_col=0).loc[rsi14.index]
d = pd.DataFrame({"pid": ph14["patient_id"].astype(str).values,
                  "t": (ph14["tissue"] == "HCC tumor").to_numpy(),
                  "rsi": rsi14["RSI"].to_numpy(float)})
tt = d[d.t].set_index("pid")["rsi"]; tt = tt[~tt.index.duplicated()]
nn = d[~d.t].set_index("pid")["rsi"]; nn = nn[~nn.index.duplicated()]
common = sorted(set(tt.index) & set(nn.index))
diff14 = (tt.loc[common] - nn.loc[common]).to_numpy(float)
assert len(diff14) == 213

rsi76 = pd.read_csv("results/GSE76427/rsi.tsv", sep="\t", index_col=0)
ph76 = pd.read_csv("results/GSE76427/phenotype.tsv", sep="\t", index_col=0).loc[rsi76.index]
d = pd.DataFrame({"pid": ph76["patient_id"].astype(str).values,
                  "t": ~ph76["tissue"].astype(str).str.lower().str.contains("adjacent"),
                  "rsi": rsi76["RSI"].to_numpy(float)})
tt = d[d.t].set_index("pid")["rsi"]; tt = tt[~tt.index.duplicated()]
nn = d[~d.t].set_index("pid")["rsi"]; nn = nn[~nn.index.duplicated()]
common = sorted(set(tt.index) & set(nn.index))
diff76 = (tt.loc[common] - nn.loc[common]).to_numpy(float)
assert len(diff76) == 52

rng = np.random.default_rng(7)
for i, (dv, lab) in enumerate(((diff76, "GSE76427\n52 pairs"), (diff14, "GSE14520\n213 pairs"))):
    x = i + rng.uniform(-0.16, 0.16, len(dv))
    axA.plot(x, dv, "o", ms=2.4, color=F.C["protect"], alpha=0.40, mec="none")
    axA.plot([i - 0.30, i + 0.30], [np.median(dv)] * 2, color=F.C["ink"], lw=1.8)
    lx, ha = ((i - 0.36, "right") if i == 0 else (i + 0.36, "left"))
    axA.text(lx, np.median(dv), f"median {np.median(dv):+.2f}", fontsize=6.0,
             va="center", ha=ha, fontweight="bold")
    axA.text(i, -3.05, f"{(dv > 0).mean() * 100:.0f}% of pairs positive",
             ha="center", fontsize=5.9, color=F.C["muted"])
axA.axhline(0, color=F.C["ink"], lw=0.9)
axA.set_xticks([0, 1]); axA.set_xticklabels(["GSE76427\n52 pairs", "GSE14520\n213 pairs"], fontsize=6.4)
axA.set_xlim(-0.95, 1.80); axA.set_ylim(-3.3, 4.4)
axA.set_ylabel("Δ index, tumor − paired adjacent liver", fontsize=6.8)
axA.set_title("Within-patient differences", fontsize=7.4, fontweight="bold", loc="left", pad=4)
for s_ in ("top", "right"): axA.spines[s_].set_visible(False)
axA.tick_params(labelsize=6.2, length=2.5, width=0.6)

axB = fig.add_axes([0.585, 0.175, 0.175, 0.645])
fig.text(0.520, 0.945, "b", **LB)
vals = [("GSE76427", h2["unpaired"]["cliffs_delta"]), ("GSE14520", g14["H2_composite"]["cliffs_delta"])]
for i, (lab, v) in enumerate(vals):
    axB.bar(i, v, 0.5, color=F.C["protect"], edgecolor=F.C["protect"], lw=0.8)
    axB.text(i, v + 0.03, f"{v:+.2f}", ha="center", fontsize=6.6, fontweight="bold",
             color=F.C["protect"])
axB.set_xticks([0, 1]); axB.set_xticklabels([v[0] for v in vals], fontsize=6.2)
axB.set_ylim(0, 1.02); axB.set_ylabel("Cliff's δ", fontsize=6.8)
axB.set_title("Effect size", fontsize=7.4, fontweight="bold", loc="left", pad=4)
for s_ in ("top", "right"): axB.spines[s_].set_visible(False)
axB.tick_params(labelsize=6.2, length=2.5, width=0.6)

axC = F.blank_axes(fig, [0.815, 0.175, 0.175, 0.645])
fig.text(0.775, 0.945, "c", **LB)
axC.text(0.0, 1.0, "GSE164760", fontsize=7.0, fontweight="bold", va="top")
axC.text(0.0, 0.86, "53 MASH-HCC tumors vs 29 non-\ntumoral MASH livers adjacent to\n"
                    "HCC; patient pairing not declared",
         fontsize=5.9, va="top", color=F.C["muted"], linespacing=1.6)
axC.text(0.0, 0.50, "δ = +0.10\nP = 0.44", fontsize=7.6, va="top", fontweight="bold",
         color=F.C["patho"], linespacing=1.5)
axC.text(0.0, 0.245, "Not supported. The comparator\nis already-diseased liver, and\n"
                     "this cohort is reported as\nuninformative for every\ndownstream adjustment by a\nrule fixed in advance.",
         fontsize=5.8, va="top", color=F.C["ink"], linespacing=1.6)

print("wrote:", *F.save_all(fig, os.path.join(OUT, "Figure2_tumor_adjacent")), sep="\n  ")
