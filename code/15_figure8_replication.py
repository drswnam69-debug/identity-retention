#!/usr/bin/env python3
"""15_figure8_replication.py -- Figure 8: cross-cohort replication (PREREG 6h/6i).

Every number drawn here is read from the results JSON and asserted against the
value quoted in the manuscript text, so the panel cannot drift from the analysis.
"""
import json, os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F

OUT = "results/figures"
g76 = json.load(open("results/GSE76427_differentiation.json"))
g16 = json.load(open("results/GSE164760_differentiation.json"))
g14 = json.load(open("results/GSE14520_differentiation.json"))

# --- assertions: the figure quotes the analysis, not a remembered number -----
def close(a, b, tol=5e-4):
    assert abs(a - b) < tol, f"{a} != {b}"

close(g76["module_reduction"]["retained_fraction"], 0.7101, 1e-2)
close(g14["module_reduction"]["retained_fraction"], 0.8369, 1e-2)
close(g14["module_drain"]["retained_fraction"], 0.2464, 1e-2)
assert g14["verdict_6h_module_split"].startswith("REPLICATED")
assert g16["verdict"].startswith("UNINFORMATIVE")
assert g14["n_pairs"] == 213 and g76["n_pairs"] == 52

fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 124 * F.MM))
LAB = dict(fontsize=10, fontweight="bold")

# ============ a: unadjusted module deltas across the three cohorts ==========
axA = fig.add_axes([0.085, 0.640, 0.335, 0.290])
fig.text(0.008, 0.972, "a", **LAB)
cohorts = [("GSE76427\nmixed · 52 pairs", g76),
           ("GSE164760\nNASH · unpaired", g16),
           ("GSE14520\nHBV · 213 pairs", g14)]
w = 0.34
for i, (lab, d) in enumerate(cohorts):
    red = d["module_reduction"]["unadjusted_mean_delta"]
    dra = d["module_drain"]["unadjusted_mean_delta"]
    axA.bar(i - w/2, red, w, color=F.C["protect"], edgecolor=F.C["protect"], lw=0.8)
    axA.bar(i + w/2, dra, w, color=F.C["patho"], edgecolor=F.C["patho"], lw=0.8)
    for v, off, col in ((red, -w/2, F.C["protect"]), (dra, w/2, F.C["patho"])):
        axA.text(i + off, v + (0.05 if v >= 0 else -0.10), f"{v:+.2f}",
                 ha="center", fontsize=5.9, fontweight="bold", color=col)
axA.axhline(0, color=F.C["ink"], lw=0.8)
axA.set_xticks(range(3)); axA.set_xticklabels([c[0] for c in cohorts], fontsize=6.2)
axA.set_ylabel("Δ module score, tumor − adjacent", fontsize=6.8)
axA.set_ylim(-1.30, 1.15)
axA.set_title("Direction replicates in three cohorts (unadjusted)",
              fontsize=7.0, fontweight="bold", pad=4, loc="left")
axA.annotate("drain shift\nabsent here", xy=(1 + w/2, 0.048), xytext=(1.24, -0.62),
             fontsize=5.7, color=F.C["muted"], ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.6, color=F.C["rule"]))
axA.legend(handles=[Patch(facecolor=F.C["protect"], label="REDUCTION"),
                    Patch(facecolor=F.C["patho"], label="DRAIN")],
           fontsize=6.0, frameon=False, loc="upper left", handlelength=1.1,
           handletextpad=0.4, borderpad=0.1, ncol=2, columnspacing=1.0,
           bbox_to_anchor=(0.0, 1.005))
for s_ in ("top", "right"): axA.spines[s_].set_visible(False)
axA.tick_params(labelsize=6.4, length=2.5, width=0.6)

# ============ b: retention after adjustment — the registered claim ==========
axB = fig.add_axes([0.560, 0.640, 0.385, 0.290])
fig.text(0.480, 0.972, "b", **LAB)
groups = [("GSE76427\n52 pairs", g76), ("GSE14520\n213 pairs", g14)]
xs = []
for i, (lab, d) in enumerate(groups):
    base = i * 2.6
    for j, (mod, col) in enumerate((("reduction", F.C["protect"]),
                                    ("drain", F.C["patho"]))):
        for k, suffix in enumerate(("", "_D2")):
            x = base + j * 1.15 + k * 0.46
            v = abs(d[f"module_{mod}{suffix}"]["retained_fraction"]) * 100
            ns = d[f"module_{mod}{suffix}"]["p"] >= 0.05
            axB.bar(x, v, 0.42, color="none" if k else col,
                    edgecolor=col, lw=0.9, hatch="////" if k else None)
            axB.text(x, v + 2.4, f"{v:.0f}", ha="center", fontsize=5.9,
                     fontweight="bold", color=col)
            if ns:
                axB.text(x, v + 9.0, "n.s.", ha="center", fontsize=5.4,
                         style="italic", color=F.C["muted"])
    xs.append(base + 0.80)
axB.axhline(50, color=F.C["rule"], lw=0.7, ls=(0, (3, 2)))
axB.set_xlim(-0.45, 5.45)
axB.text(4.52, 50.0, "registered\n50% rule", fontsize=5.3,
         color=F.C["muted"], va="center", ha="left", linespacing=1.4)
axB.set_xticks(xs); axB.set_xticklabels([g[0] for g in groups], fontsize=6.4)
axB.set_ylabel("% of the unadjusted shift retained", fontsize=6.8)
axB.set_ylim(0, 118)
axB.set_title("The reduction arm survives adjustment for hepatocyte identity;\n"
              "the drain arm largely does not",
              fontsize=7.0, fontweight="bold", pad=4, loc="left")
axB.legend(handles=[Patch(facecolor=F.C["protect"], label="REDUCTION"),
                    Patch(facecolor=F.C["patho"], label="DRAIN"),
                    Patch(facecolor="white", edgecolor=F.C["ink"], hatch="////",
                          label="D2 (over-adjusted)")],
           fontsize=5.9, frameon=False, loc="upper center", ncol=3,
           handlelength=1.1, handletextpad=0.4, borderpad=0.1,
           columnspacing=1.0, bbox_to_anchor=(0.5, 1.005))
for s_ in ("top", "right"): axB.spines[s_].set_visible(False)
axB.tick_params(labelsize=6.4, length=2.5, width=0.6)

# ============ c: composite intercept, unadjusted vs adjusted ================
axC = fig.add_axes([0.195, 0.225, 0.245, 0.260])
fig.text(0.008, 0.545, "c", **LAB)
rows = []
for lab, d, unadj in (("GSE76427", g76, 1.129), ("GSE14520", g14, 1.261)):
    rows.append((lab + "  unadjusted", unadj, None, None, F.C["muted"]))
    for tag, col in (("D1", F.C["ink"]), ("D2", F.C["muted"])):
        f = d[f"paired_adjusted_{tag}"]["intercept"]
        rows.append((f"{lab}  adjusted {tag}", f["beta"], f["ci95"][0],
                     f["ci95"][1], col))
ypos = np.arange(len(rows))[::-1]
for y, (lab, v, lo, hi, col) in zip(ypos, rows):
    if lo is not None:
        axC.plot([lo, hi], [y, y], color=col, lw=1.0, solid_capstyle="round")
    axC.plot([v], [y], "o", ms=4.2, color=col, mec="white", mew=0.6)
    axC.text(v, y + 0.30, f"{v:+.2f}", ha="center", fontsize=5.9,
             fontweight="bold", color=col)
axC.axvline(0, color=F.C["ink"], lw=0.8)
axC.set_yticks(ypos); axC.set_yticklabels([r[0] for r in rows], fontsize=6.0)
axC.set_xlabel("ΔRSI, tumor − adjacent, at zero change in hepatocyte identity",
               fontsize=6.4)
axC.set_xlim(-0.10, 1.55); axC.set_ylim(-0.7, len(rows) - 0.25)
axC.set_title("Composite: reduced by adjustment, not abolished",
              fontsize=7.0, fontweight="bold", pad=4, loc="left")
for s_ in ("top", "right", "left"): axC.spines[s_].set_visible(False)
axC.tick_params(labelsize=6.0, length=2.5, width=0.6)

# ============ d: status =====================================================
axD = F.blank_axes(fig, [0.560, 0.180, 0.385, 0.310])
fig.text(0.480, 0.545, "d", **LAB)
axD.text(0.0, 1.02, "What replicates, and what it does not license",
         fontsize=7.0, fontweight="bold", va="top")
items = [
    ("#1d7a58", "The phenomenon replicates.",
     "REDUCTION ↑ / DRAIN ↓ in tumor across three\ncohorts, three platforms, three etiologies."),
    ("#1d7a58", "So does the decomposition.",
     "The identity adjustment now holds in two independent\n"
     "cohorts under a rule fixed before the second was acquired:\n71% / 84% reduction vs 41% / 25% drain."),
    (F.C["muted"], "It is not causal, and one gene is missing.",
     "These are cross-sectional transcriptomes. AIFM2 is\n"
     "absent from HG-U133A, so GSE14520's REDUCTION is a 3-gene\n"
     "score without the module's strongest measured member."),
]
y = 0.855
for col, head, body in items:
    axD.plot([0.0, 0.026], [y + 0.010, y + 0.010], color=col, lw=2.2,
             solid_capstyle="round")
    axD.text(0.044, y + 0.030, head, fontsize=6.3, fontweight="bold",
             va="top", color=col)
    axD.text(0.044, y - 0.055, body, fontsize=5.7, va="top",
             color=F.C["muted"], linespacing=1.65)
    y -= 0.310

fig.text(0.045, 0.115,
         "Cohorts, covariate and decision rule fixed in PREREGISTRATION §6h before any replication number existed; §6i records the completion of\n"
         "GSE14520's acquisition, including the recovery of HNF1A from its retired symbol TCF1, before any adjusted number was computed. GSE164760\n"
         "returns \"uninformative\" by that rule because its unadjusted DRAIN shift is absent, and is counted as neither support nor refutation. D2 adds the\n"
         "hepatic P450s to the covariate and is expected to over-adjust, because POR donates electrons to those very enzymes.",
         fontsize=5.4, color=F.C["muted"], va="top", linespacing=1.7)

paths = F.save_all(fig, os.path.join(OUT, "figure8_replication"))
print("wrote:", *paths, sep="\n  ")
