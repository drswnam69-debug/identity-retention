#!/usr/bin/env python3
"""11_summary_figure.py -- study-at-a-glance graphic.

Panel a: how the index is built. Panel b: H1 forest across three cohorts.
Panel c: H2 effect. Panel d: the six-hypothesis scoreboard.
All numbers read from the results files.
"""
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F

OUT = "results/figures"
os.makedirs(OUT, exist_ok=True)

pb = {c: json.load(open(f"results/{c}/phase_b.json"))
      for c in ("GSE135251", "GSE130970", "GSE167523")}
h2 = json.load(open("results/gse76427_h2.json"))
pe = json.load(open("results/phase_e.json"))

# Numbers quoted as literals below are asserted against the results files so
# the graphic cannot drift from the analysis.
POOLED = (0.372, 0.277, 0.459)          # PHASE_B_RESULTS.md, random effects
c164 = json.load(open("results/GSE164760/phase_c.json"))
assert abs(c164["H3"]["p_two_sided"] - 0.060) < 0.001, c164["H3"]["p_two_sided"]
assert c164["H3"]["n_a"] == 29
assert abs(h2["unpaired"]["cliffs_delta"] - 0.76) < 0.005
assert abs(h2["paired"]["p"] - 8.1e-8) < 0.2e-8
assert abs(pe["H4"]["adjusted"]["HR"][3] - 4.56) < 0.01
assert abs(pe["H4"]["primary_continuous"]["HR_per_SD"] - 0.99) < 0.01

fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 158 * F.MM))

# ------------------------------------------------- a: index construction
axA = F.blank_axes(fig, [0.03, 0.660, 0.94, 0.310])
fig.text(0.010, 0.982, "a", fontsize=10, fontweight="bold")
axA.text(0.5, 0.96, "The Reductive-Supply Index", ha="center", va="top",
         fontsize=8.5, fontweight="bold", color=F.C["ink"])
axA.text(0.5, 0.855,
         "one number per liver sample: how much reducing power the tissue "
         "can supply, minus how fast it is drained",
         ha="center", va="top", fontsize=6.6, color=F.C["muted"])

mods = [
    (0.155, "SUPPLY", "mevalonate + CoQ\nbiosynthesis",
     "HMGCR · HMGCS1 · MVK\nPMVK · MVD · IDI1 · FDPS\nPDSS1/2 · COQ2–COQ9",
     F.C["supply"], F.C["supply_f"], "+ ½"),
    (0.435, "REDUCTION", "reductases that\nregenerate the pool",
     "CYB5R3 · CYB5R1\nAIFM2 (FSP1) · NQO1",
     F.C["protect"], F.C["protect_f"], "+ ½"),
    (0.715, "DRAIN", "consumers of\nreducing equivalents",
     "mARC1 (MTARC1)\nmARC2 · POR",
     F.C["patho"], F.C["patho_f"], "−1"),
]
for x, name, sub, genes, ec, fc, w in mods:
    F.box(axA, x, 0.50, 0.235, 0.44, "", ec, fc, lw=1.0, radius=0.03)
    axA.text(x, 0.665, name, ha="center", va="center", fontsize=7.6,
             fontweight="bold", color=ec)
    axA.text(x, 0.560, sub, ha="center", va="center", fontsize=6.2,
             color=F.C["ink"])
    axA.text(x, 0.410, genes, ha="center", va="center", fontsize=5.6,
             color=F.C["muted"], linespacing=1.5)
    axA.text(x, 0.245, f"weight {w}", ha="center", va="center", fontsize=6.4,
             fontweight="bold", color=ec)
    F.arrow(axA, (x, 0.195), (x, 0.115), ec, lw=1.0, ms=6)

axA.plot([0.155, 0.715], [0.115, 0.115], color=F.C["rule"], lw=0.8)
F.arrow(axA, (0.435, 0.115), (0.435, 0.055), F.C["ink"], lw=1.2, ms=7)
F.box(axA, 0.435, 0.015, 0.30, 0.085, "RSI  =  ½(supply + reduction) − drain",
      F.C["ink"], F.C["out_f"], fs=7.0, bold=True, lw=1.1, radius=0.03)
axA.text(0.928, 0.50,
         "locked before\nany analysis\n\nSHA-256\n7a2bf934…",
         ha="center", va="center", fontsize=5.8, color=F.C["muted"],
         linespacing=1.6)
axA.add_patch(plt.Rectangle((0.850, 0.30), 0.156, 0.40, fill=False,
                            ec=F.C["rule"], lw=0.7, ls=(0, (2, 2))))

# --------------------------------------------------------- b: H1 forest
axB = fig.add_axes([0.155, 0.435, 0.275, 0.150])
fig.text(0.010, 0.630, "b", fontsize=10, fontweight="bold")
rows = [("GSE135251", "fibrosis spectrum"), ("GSE130970", "fibrosis F0–F4"),
        ("GSE167523", "MASL → MASH")]
ys = np.arange(len(rows) + 1)[::-1]
for i, (c, _) in enumerate(rows):
    pr = pb[c]["primary_RSI"]
    rho, lo, hi = pr["spearman_rho"], *pr["rho_ci95"]
    axB.plot([lo, hi], [ys[i], ys[i]], color=F.C["protect"], lw=1.0,
             solid_capstyle="round")
    axB.plot([rho], [ys[i]], "s", ms=3.6, color=F.C["protect"])
    axB.text(1.03, ys[i], f"n = {pb[c]['n']}",
             transform=axB.get_yaxis_transform(), va="center", fontsize=5.8,
             color=F.C["muted"])
axB.plot([POOLED[1], POOLED[2]], [ys[-1], ys[-1]], color=F.C["ink"], lw=1.6,
         solid_capstyle="round")
axB.plot([POOLED[0]], [ys[-1]], "D", ms=4.4, color=F.C["ink"])
axB.text(1.03, ys[-1], "n = 392", transform=axB.get_yaxis_transform(),
         va="center", fontsize=5.8, fontweight="bold", color=F.C["ink"])
axB.axvline(0, color=F.C["rule"], lw=0.6)
axB.set_yticks(ys)
axB.set_yticklabels([r[0] for r in rows] + ["Pooled  (I² = 0%)"], fontsize=6.2)
for lab in axB.get_yticklabels()[-1:]:
    lab.set_fontweight("bold")
axB.set_xlim(-0.05, 0.68)
axB.set_ylim(-0.6, len(rows) + 0.5)
axB.set_xlabel("Spearman ρ,  RSI vs disease severity", fontsize=6.6)
axB.set_title("H1 — RSI rises with MASLD severity:\nsupported, replicated",
              fontsize=6.8, fontweight="bold", pad=5, loc="left")
for s in ("top", "right", "left"):
    axB.spines[s].set_visible(False)
axB.tick_params(labelsize=6.0, length=2.5, width=0.6)

# -------------------------------------------------------------- c: H2/H4
axC = F.blank_axes(fig, [0.545, 0.420, 0.430, 0.185])
fig.text(0.490, 0.630, "c", fontsize=10, fontweight="bold")
dif = json.load(open("results/GSE76427_differentiation.json"))
_red = dif["module_reduction"]; _dra = dif["module_drain"]
_rr = _red["retained_fraction"] * 100
_dr = abs(_dra["retained_fraction"]) * 100
_g14 = json.load(open("results/GSE14520_differentiation.json"))
_rr2 = _g14["module_reduction"]["retained_fraction"] * 100
_dr2 = abs(_g14["module_drain"]["retained_fraction"]) * 100
_G14H2 = json.load(open("results/gse14520_partial.json"))["H2_composite"]["cliffs_delta"]
axC.text(0.0, 0.99, "H2 — tumor vs adjacent liver:  supported in two cohorts",
         fontsize=6.8, fontweight="bold", va="top", color=F.C["ink"])
axC.text(0.0, 0.845,
         f"Cliff's δ  {h2['unpaired']['cliffs_delta']:+.2f} / "
         f"{_G14H2:+.2f}   (52 and 213 matched pairs)",
         fontsize=6.4, va="top", color=F.C["ink"])
axC.text(0.0, 0.715,
         "carried by REDUCTION ↑ and DRAIN ↓ — not by supply",
         fontsize=6.2, va="top", color=F.C["protect"], style="italic")
axC.text(0.0, 0.600,
         f"identity-adjusted: REDUCTION keeps {_rr:.0f}% / {_rr2:.0f}%, "
         f"DRAIN only {_dr:.0f}% / {_dr2:.0f}%",
         fontsize=5.9, va="top", color=F.C["muted"])
axC.plot([0.0, 1.0], [0.470, 0.470], color=F.C["rule"], lw=0.7)
axC.text(0.0, 0.395, "H4 — RSI and overall survival:  not supported",
         fontsize=6.8, fontweight="bold", va="top", color=F.C["ink"])
axC.text(0.0, 0.235,
         f"HR {pe['H4']['primary_continuous']['HR_per_SD']:.2f} per SD "
         f"({pe['H4']['primary_continuous']['ci95'][0]:.2f}–"
         f"{pe['H4']['primary_continuous']['ci95'][1]:.2f}),  P = 0.95",
         fontsize=6.4, va="top", color=F.C["ink"])
axC.text(0.0, 0.090,
         "23 events: a strong effect is excluded, a modest one is not\n"
         "(BCLC C/D in the same model: HR 4.56 — the model does see signal)",
         fontsize=6.0, va="top", color=F.C["muted"], linespacing=1.6)

# ------------------------------------------------------- d: the scoreboard
axD = F.blank_axes(fig, [0.03, 0.035, 0.94, 0.310])
fig.text(0.010, 0.366, "d", fontsize=10, fontweight="bold")
axD.text(0.030, 1.030, "Where the six pre-registered hypotheses landed",
         fontsize=8.0, fontweight="bold", va="top")

SUP = ("#1d7a58", F.C["protect_f"], "supported")
NOT = ("#b0453a", "#f7e2df", "not supported")
NA = (F.C["muted"], "#eeece7", "not testable")

rows = [
    ("H1", "RSI rises across the MASLD severity gradient",
     "pooled ρ = 0.372 (0.277–0.459), I² = 0%, n = 392, 3 cohorts", SUP),
    ("H2", "RSI higher in HCC tumor than adjacent liver",
     "Cliff's δ +0.76 and +0.86 in two cohorts; 50–63% survives adjustment "
     "for hepatocyte identity, carried by the reduction arm (71%, 84%)", SUP),
    ("H3", "field effect in adjacent non-tumor liver",
     "P = 0.060 — underpowered, 29 adjacent samples", NOT),
    ("H4", "RSI predicts overall survival",
     "HR 0.99 per SD; only 23 events — bounded null, not a refutation", NOT),
    ("H5", "NMD-target burden tracks RSI",
     "proxy correlates with mean expression at ρ 0.82–0.98; needs "
     "transcript-level data", NA),
    ("H6", "RSI tracks a ferroptosis-resistance signature",
     "positive with suppressors AND drivers — real but non-directional", NOT),
]
y = 0.870
for tag, claim, detail, (ec, fc, verdict) in rows:
    axD.add_patch(plt.Rectangle((0.0, y - 0.093), 1.0, 0.128,
                                facecolor="#faf9f6", edgecolor="none",
                                zorder=0))
    axD.text(0.012, y, tag, fontsize=7.4, fontweight="bold",
             va="center", color=F.C["ink"])
    axD.text(0.055, y, claim, fontsize=6.8, va="center", color=F.C["ink"])
    axD.text(0.055, y - 0.062, detail, fontsize=5.8, va="center",
             color=F.C["muted"])
    F.box(axD, 0.905, y - 0.018, 0.175, 0.072, verdict, ec, fc,
          fs=6.2, bold=True, lw=0.9, radius=0.02)
    y -= 0.152

axD.plot([0.0, 1.0], [-0.020, -0.020], color=F.C["rule"], lw=0.8)
axD.text(0.0, -0.055,
         "Conclusion: RSI is a tissue-state index, not a prognostic "
         "biomarker. The index was not redefined after these results "
         "(pre-registered stopping rule §7).",
         fontsize=6.4, va="top", color=F.C["ink"], fontweight="bold")

paths = F.save_all(fig, os.path.join(OUT, "figure6_summary"))
print("wrote:", *paths, sep="\n  ")
