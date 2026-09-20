#!/usr/bin/env python3
"""71_figures.py -- the three figures, drawn from PANEL_COHERENCE.json only."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_cands = [_here, _os.path.join(_here, "rsi", "code"), _os.path.join(_here, "code"),
          _os.path.join(_os.path.dirname(_here), "code")]
if _os.environ.get("IR_ROOT"):
    _cands.insert(0, _os.path.join(_os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        if _c not in _sys.path:
            _sys.path.insert(0, _c)
        break
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, \
    DATA as IR_DATA, CODE as IR_CODE, FIGURES as IR_FIGURES, DOCS as IR_DOCS
from figstyle import C, save_all                      # noqa: E402
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Rectangle               # noqa: E402

R = json.load(open(f"{IR_RESULTS}/PANEL_COHERENCE.json"))
OUT = f"{IR_FIGURES}/panel_coherence"
os.makedirs(OUT, exist_ok=True)

DS = list(R["datasets"])
SHORT = {"GSE14520 (array transcriptome)": "GSE14520\narray\n213 pairs",
         "GSE76427 (array transcriptome)": "GSE76427\narray\n52 pairs",
         "TCGA-LIHC (RNA sequencing)":     "TCGA-LIHC\nRNA-seq\n50 pairs",
         "Gao 2019 proteome":              "Gao 2019\nproteome\n159 pairs",
         "Jiang 2019 proteome":            "Jiang 2019\nproteome\n124 pairs"}
LEVEL = {n: R["datasets"][n]["level"] for n in DS}
CLASSES = R["panel_classes"]
GENES = [g for v in CLASSES.values() for g in v]
CLS_OF = {g: c for c, v in CLASSES.items() for g in v}
CLS_COLOR = {"secreted plasma protein": C["patho"],
             "metabolic enzyme":        C["supply"],
             "transcription factor":    C["protect"],
             "surface receptor":        C["hub"]}

# Diverging pair from the validated palette: blue and red poles, neutral gray
# midpoint. A hue at the midpoint would read as a value rather than as nothing.
DIV = LinearSegmentedColormap.from_list(
    "bl_gy_rd", ["#0d366b", "#256abf", "#86b6ef", "#f0efec",
                 "#eda48f", "#d03b3b", "#8f2424"])


def fig1():
    fig = plt.figure(figsize=(6.85, 4.6)); fig.patch.set_facecolor("white")
    L, B, W, H = 0.305, 0.175, 0.525, 0.70
    ax = fig.add_axes([L, B, W, H])
    M = np.full((len(GENES), len(DS)), np.nan)
    for j, n in enumerate(DS):
        pm = R["datasets"][n]["panel_members"]
        for i, g in enumerate(GENES):
            if g in pm:
                M[i, j] = pm[g]["z"]
    im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(vcenter=0, vmin=-1.7, vmax=1.7),
                   aspect="auto")
    for i in range(len(GENES)):
        for j in range(len(DS)):
            if np.isnan(M[i, j]):
                ax.text(j, i, "ns", ha="center", va="center", fontsize=6,
                        color=C["muted"]); continue
            ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(M[i, j]) > 0.95 else C["ink"])
    ax.set_xticks(range(len(DS)))
    ax.set_xticklabels([SHORT[n] for n in DS], fontsize=6.6, linespacing=1.35)
    ax.set_yticks(range(len(GENES)))
    ax.set_yticklabels(GENES, fontsize=6.8, style="italic")
    for t, g in zip(ax.get_yticklabels(), GENES):
        t.set_color(CLS_COLOR[CLS_OF[g]])
    ax.tick_params(length=0, pad=2)
    for s in ax.spines.values():
        s.set_color(C["rule"])
    # the line between transcriptomes and proteomes
    k = sum(1 for n in DS if LEVEL[n] == "transcriptome")
    ax.axvline(k - 0.5, color=C["ink"], lw=1.4)
    ax.text((k - 1) / 2, -0.95, "transcriptome", ha="center", fontsize=7.4,
            fontweight="bold", color=C["ink"])
    ax.text((k + len(DS) - 1) / 2, -0.95, "proteome", ha="center", fontsize=7.4,
            fontweight="bold", color=C["ink"])
    # class brackets
    start = 0
    for cls, gs in CLASSES.items():
        n = len(gs)
        ax.plot([-0.62, -0.62], [start - 0.38, start + n - 0.62],
                color=CLS_COLOR[cls], lw=2.4, clip_on=False,
                solid_capstyle="butt", transform=ax.get_yaxis_transform(
                    which="grid") if False else ax.transData)
        fig.text(0.175, B + H * (1 - (start + n / 2) / len(GENES)),
                 cls.replace(" ", "\n"), ha="right", va="center", fontsize=6.9,
                 color=CLS_COLOR[cls], linespacing=1.3)
        start += n
    cb = fig.add_axes([0.86, 0.30, 0.016, 0.42])
    bar = fig.colorbar(im, cax=cb)
    bar.set_label("paired shift, tumor minus adjacent (z)", fontsize=6.8, labelpad=3)
    bar.ax.tick_params(labelsize=6.6, length=2)
    bar.outline.set_edgecolor(C["rule"])
    fig.text(0.02, 0.965, "The panel holds together at transcript level and splits at protein level",
             fontsize=8.6, fontweight="bold", color=C["ink"], va="top")
    return save_all(fig, f"{OUT}/Figure1_panel_members")


def fig2():
    fig = plt.figure(figsize=(6.85, 3.05)); fig.patch.set_facecolor("white")
    # (a) the composite, and the verdict it produces
    ax = fig.add_axes([0.095, 0.245, 0.345, 0.545])
    y = np.arange(len(DS))[::-1]
    for yi, n in zip(y, DS):
        t = R["datasets"][n]["panel_composite_sample_test"]
        col = C["supply"] if LEVEL[n] == "transcriptome" else C["patho"]
        ax.barh(yi, t["mean_paired_delta"], height=0.5, color=col, alpha=0.85,
                edgecolor="none")
        xoff = 0.03 if t["mean_paired_delta"] > 0 else -0.03
        ax.text(t["mean_paired_delta"] + xoff, yi,
                f"P = {t['wilcoxon_p']:.0e}".replace("e-0", "e-"), fontsize=6.4,
                va="center", ha="left" if t["mean_paired_delta"] > 0 else "right",
                color=C["muted"])
    ax.axvline(0, color=C["ink"], lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([n.split(" (")[0].replace(" 2019 proteome", "") for n in DS],
                       fontsize=7)
    ax.set_xlim(-1.12, 0.46)
    ax.set_xlabel("composite panel score, paired shift (z)", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C["rule"])
    ax.set_title("(a)  Every dataset gives a significant verdict,\n"
                 "      and the two proteomes disagree on the sign",
                 fontsize=7.6, fontweight="bold", color=C["ink"], loc="left",
                 pad=5, linespacing=1.3)

    # (b) how much of the members' movement survives into the composite
    ax2 = fig.add_axes([0.585, 0.245, 0.375, 0.545])
    rng = np.random.default_rng(20260920)
    for j, n in enumerate(DS):
        sets = R["datasets"][n]["gene_set_coherence"]
        v = [s["cancellation"] for s in sets.values() if s["cancellation"] is not None]
        ax2.scatter(j + rng.uniform(-0.22, 0.22, len(v)), v, s=4,
                    facecolor=C["muted"], alpha=0.3, linewidth=0, zorder=2)
        pm = [m["z"] for m in R["datasets"][n]["panel_members"].values()]
        a = np.array(pm)
        panel = abs(a.mean()) / np.abs(a).mean()
        col = C["supply"] if LEVEL[n] == "transcriptome" else C["patho"]
        ax2.plot([j], [panel], "D", color=col, ms=6.5, zorder=4,
                 markeredgecolor="white", markeredgewidth=0.8)
        ax2.text(j, panel + 0.05, f"{panel:.2f}", ha="center", fontsize=6.6,
                 color=col, fontweight="bold")
    k = sum(1 for n in DS if LEVEL[n] == "transcriptome")
    ax2.axvline(k - 0.5, color=C["ink"], lw=1.0)
    ax2.set_xticks(range(len(DS)))
    ax2.set_xticklabels([n.split(" (")[0].replace(" 2019 proteome", "") for n in DS],
                        fontsize=6.6, rotation=20, ha="right")
    ax2.set_ylim(-0.03, 1.05)
    ax2.set_ylabel("share of member movement that\nsurvives into the composite",
                   fontsize=7, labelpad=2, linespacing=1.3)
    ax2.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax2.spines[s].set_color(C["rule"])
    ax2.set_title("(b)  The panel (diamond) against 114 published\n"
                  "      liver signatures scored in the same data (dots)",
                  fontsize=7.6, fontweight="bold", color=C["ink"], loc="left",
                  pad=5, linespacing=1.3)
    return save_all(fig, f"{OUT}/Figure2_composite_and_coherence")


def fig3():
    fig = plt.figure(figsize=(6.85, 2.7)); fig.patch.set_facecolor("white")
    rows = [("D1 secreted plasma proteins", "secreted plasma protein", C["patho"]),
            ("D1 metabolic enzymes", "metabolic enzyme", C["supply"]),
            ("plasma proteins outside the panel", None, C["hub"]),
            ("erythrocyte proteins", None, C["out"])]
    ax = fig.add_axes([0.305, 0.265, 0.355, 0.575])
    for i, (lab, cls, col) in enumerate(rows[::-1]):
        for j, n in enumerate(DS):
            d = R["datasets"][n]
            if cls:
                v = d["panel_by_class"].get(cls, {}).get("mean_z")
            else:
                v = (d["off_panel_plasma_mean_z"] if "outside" in lab
                     else d["erythrocyte_mean_z"])
            if v is None:
                continue
            tx = LEVEL[n] == "transcriptome"
            ax.plot([v], [i + (0.11 if tx else -0.11)], "o" if tx else "s",
                    color=col, ms=4.6, markeredgecolor="white",
                    markeredgewidth=0.6, alpha=0.5 if tx else 1.0)
    ax.axvline(0, color=C["ink"], lw=0.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=6.9)
    for t, r in zip(ax.get_yticklabels(), rows[::-1]):
        t.set_color(r[2])
    ax.set_xlabel("mean paired shift (z)", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C["rule"])
    ax.set_ylim(-0.6, len(rows) - 0.35)
    ax.set_title("(a)  Blood content does not explain the split", fontsize=7.6,
                 fontweight="bold", color=C["ink"], loc="left", pad=5)
    ax.plot([], [], "o", color=C["muted"], ms=4.6, alpha=0.5,
            label="transcriptome (3 datasets)")
    ax.plot([], [], "s", color=C["muted"], ms=4.6, label="proteome (2 datasets)")
    ax.legend(fontsize=6.4, frameon=False, ncol=2, handletextpad=0.35,
              columnspacing=1.1, loc="upper center",
              bbox_to_anchor=(0.5, -0.20))

    ax2 = fig.add_axes([0.765, 0.265, 0.195, 0.575])
    ns = R["jiang_normalization_sensitivity"]
    labs = ["raw iBAQ", "median\ncentered", "quantile\nnormalized"]
    x = np.arange(3)
    ax2.plot(x, [ns[k]["secreted_mean_z"] for k in ("raw", "median", "quantile")],
             "s-", color=C["patho"], ms=5, lw=1.4, markeredgecolor="white",
             markeredgewidth=0.6, label="secreted")
    ax2.plot(x, [ns[k]["metabolic_mean_z"] for k in ("raw", "median", "quantile")],
             "s-", color=C["supply"], ms=5, lw=1.4, markeredgecolor="white",
             markeredgewidth=0.6, label="metabolic")
    ax2.axhline(0, color=C["ink"], lw=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels(labs, fontsize=6.4, linespacing=1.25)
    ax2.set_xlim(-0.45, 2.45)
    ax2.set_ylabel("mean paired shift (z)", fontsize=7, labelpad=2)
    ax2.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax2.spines[s].set_color(C["rule"])
    ax2.legend(fontsize=6.4, frameon=False, loc="center right", handletextpad=0.4)
    ax2.set_title("(b)  Nor does the processing\n      (Jiang, three ways)",
                  fontsize=7.6, fontweight="bold", color=C["ink"], loc="left",
                  pad=5, linespacing=1.25)
    return save_all(fig, f"{OUT}/Figure3_controls")


for f in (fig1, fig2, fig3):
    print("\n".join("  " + p for p in f()))
