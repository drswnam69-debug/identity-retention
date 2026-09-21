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

def pfmt(p):
    """One spelling of a P value, used in every panel and matching the text."""
    if p is None:
        return ""
    if p >= 1e-3:
        return f"P = {p:.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"P = {m} \u00d7 10$^{{{int(e)}}}$"


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
# Four classes need four hues. The first three are the validated categorical
# trio; the fourth was added and the set revalidated:
#   node scripts/validate_palette.js "#4a3aa7,#199e70,#d95926,#256abf" --mode light
#   -> PASS on all six checks (worst adjacent CVD dE 9.4 deutan, normal dE 26.5).
RECEPTOR = "#256abf"
CLS_COLOR = {"secreted plasma protein": C["patho"],
             "metabolic enzyme":        C["supply"],
             "transcription factor":    C["protect"],
             "surface receptor":        RECEPTOR}

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
    NP = [[None] * len(DS) for _ in GENES]
    FULL = [R["datasets"][n]["n_pairs"] for n in DS]
    for j, n in enumerate(DS):
        pm = R["datasets"][n]["panel_members"]
        for i, g in enumerate(GENES):
            if g in pm:
                M[i, j] = pm[g]["z"]
                NP[i][j] = pm[g]["n_pairs"]
    NP = np.array(NP, dtype=object)
    im = ax.imshow(M, cmap=DIV, norm=TwoSlopeNorm(vcenter=0, vmin=-1.7, vmax=1.7),
                   aspect="auto")
    for i in range(len(GENES)):
        for j in range(len(DS)):
            if np.isnan(M[i, j]):
                ax.text(j, i, "n.q.", ha="center", va="center", fontsize=5.8,
                        color=C["muted"]); continue
            ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=6,
                    color="white" if abs(M[i, j]) > 0.95 else C["ink"])
            if NP[i, j] is not None and NP[i, j] < FULL[j]:
                ax.plot([j + 0.355], [i - 0.32], marker="o", ms=1.9,
                        color="white" if abs(M[i, j]) > 0.95 else C["ink"],
                        markeredgewidth=0)
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
    fig.text(0.02, 0.965, "The panel mostly holds together at transcript level and splits at protein level",
             fontsize=8.6, fontweight="bold", color=C["ink"], va="top")
    return save_all(fig, f"{OUT}/Figure1_panel_members")


def fig2():
    """Level is carried by fill, never by hue: the hues in this paper mean
    protein class, and reusing them for transcriptome and proteome would make
    the same color mean two things in two figures."""
    fig = plt.figure(figsize=(6.85, 3.15)); fig.patch.set_facecolor("white")
    TXC, PRC = C["hub"], "white"

    def style(n):
        return dict(color=TXC if LEVEL[n] == "transcriptome" else PRC,
                    edgecolor=TXC,
                    linewidth=0.0 if LEVEL[n] == "transcriptome" else 0.9)

    # (a) the composite, and the verdict it produces
    ax = fig.add_axes([0.105, 0.300, 0.240, 0.495])
    y = np.arange(len(DS))[::-1]
    for yi, n in zip(y, DS):
        t = R["datasets"][n]["panel_composite_sample_test"]
        d = t["mean_paired_delta"]
        ax.barh(yi, d, height=0.5, **style(n))
        # every P label in one column, so the longest bar cannot push its own
        # label into the axis labels
        ax.text(0.15, yi, pfmt(t["wilcoxon_p"]), fontsize=6.2, va="center",
                ha="left", color=C["muted"])
    # The same Jiang composite under the other five pipelines built from the
    # same deposited file. Leaving them unlabeled would let this panel show a
    # null value while its title spoke only of the significant one.
    ji = DS.index("Jiang 2019 proteome"); yj = y[ji]
    ns = R["jiang_normalization_sensitivity"]
    alt = []
    for k in ("raw", "quantile", "median"):
        for key in ("composite_sample_test", "composite_sample_test_unfiltered"):
            alt.append((ns[k][key]["mean_paired_delta"], ns[k][key]["wilcoxon_p"]))
    ya = yj - 0.85
    for v, _pv in alt:
        ax.plot([v], [ya], marker="o", ms=3.0, markerfacecolor="white",
                markeredgecolor=TXC, markeredgewidth=0.8, zorder=5)
    lo, hi = min(a[0] for a in alt), max(a[0] for a in alt)
    plo = max(a[1] for a in alt); phi = min(a[1] for a in alt)
    ax.plot([lo, hi], [ya, ya], color=TXC, lw=0.7, zorder=4)
    ax.text(0.20, ya, f"six pipelines, one file\nP = {phi:.3f} to {plo:.2f}",
            fontsize=5.8, color=C["muted"], va="center", ha="left",
            linespacing=1.3, clip_on=False)
    ax.axvline(0, color=C["ink"], lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([n.split(" (")[0].replace(" 2019 proteome", "") for n in DS],
                       fontsize=7)
    ax.set_ylim(-1.55, len(DS) - 0.45)
    ax.set_xlim(-0.92, 1.28)
    ax.set_xlabel("composite panel score, paired shift (z)", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(C["rule"])
    ax.set_title("(a)  One verdict per dataset, opposite in\n"
                 "      sign in the two proteomes",
                 fontsize=7.0, fontweight="bold", color=C["ink"], loc="left",
                 pad=5, linespacing=1.3)

    # (b) surviving share: the panel against the 114 sets scored at BOTH levels
    shared = [r["set"] for r in R["coherence_by_level"]["per_set"]]
    nshared = len(shared)
    ax2 = fig.add_axes([0.450, 0.300, 0.265, 0.495])
    rng = np.random.default_rng(20260920)
    pw = R["panel_within_shared_sets"]["per_dataset"]
    for j, n in enumerate(DS):
        sets = R["datasets"][n]["gene_set_coherence"]
        v = [sets[k]["cancellation"] for k in shared
             if sets[k]["cancellation"] is not None]
        ax2.scatter(j + rng.uniform(-0.22, 0.22, len(v)), v, s=4,
                    facecolor=C["muted"], alpha=0.30, linewidth=0, zorder=2)
        panel = pw[n]["panel_surviving_share"]
        ax2.plot([j], [panel], "D", ms=6.5, zorder=4,
                 markerfacecolor=TXC if LEVEL[n] == "transcriptome" else PRC,
                 markeredgecolor=TXC, markeredgewidth=0.9)
        # keep the label off the metabolic square when the two nearly coincide
        above = panel < 0.85
        ax2.text(j, panel + (0.055 if above else -0.065), f"{panel:.2f}",
                 ha="center", va="bottom" if above else "top", fontsize=6.6,
                 color=C["ink"], fontweight="bold")
        # the metabolic class alone, which is what does not move
        met = R["class_within_shared_sets"][n]["metabolic enzyme"]["surviving_share"]
        ax2.plot([j], [met], "s", ms=5.2, zorder=5,
                 markerfacecolor=CLS_COLOR["metabolic enzyme"]
                 if LEVEL[n] == "transcriptome" else "white",
                 markeredgecolor=CLS_COLOR["metabolic enzyme"], markeredgewidth=1.0)
    k = sum(1 for n in DS if LEVEL[n] == "transcriptome")
    ax2.axvline(k - 0.5, color=C["ink"], lw=1.0)
    ax2.set_xticks(range(len(DS)))
    ax2.set_xticklabels([n.split(" (")[0].replace(" 2019 proteome", "") for n in DS],
                        fontsize=6.6, rotation=20, ha="right")
    ax2.set_ylim(-0.03, 1.12)
    ax2.set_ylabel("share of member movement that\nsurvives into their own mean",
                   fontsize=7, labelpad=2, linespacing=1.3)
    ax2.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax2.spines[sp].set_color(C["rule"])
    ax2.set_title(f"(b)  The panel (diamond) and its metabolic class\n"
                  f"      (square) against the {nshared} signatures (dots)",
                  fontsize=7.0, fontweight="bold", color=C["ink"], loc="left",
                  pad=5, linespacing=1.3)

    # (c) how far each set falls between the two levels
    ax3 = fig.add_axes([0.822, 0.300, 0.163, 0.495])
    pw2 = R["panel_within_shared_sets"]
    drops = np.array([r["transcriptome_cancellation"] - r["proteome_cancellation"]
                      for r in R["coherence_by_level"]["per_set"]])
    ax3.hist(drops, bins=22, color=C["hub_f"], edgecolor=C["hub"], linewidth=0.5)
    pd_ = pw2["panel_drop"]
    ax3.axvline(pd_, color=C["block"], lw=1.6, zorder=5)
    ax3.set_xlim(min(drops.min(), -0.05) - 0.06, pd_ + 0.13)
    ax3.text(pd_ - 0.03, ax3.get_ylim()[1] * 0.96, f"panel\n{pd_:.2f}",
             fontsize=6.3, color=C["block"], ha="right", va="top",
             fontweight="bold", linespacing=1.25)
    ax3.set_xlabel("fall in surviving share,\ntranscriptome to proteome",
                   fontsize=7, labelpad=2, linespacing=1.25)
    ax3.set_ylabel(f"signatures (of {nshared})", fontsize=7, labelpad=2)
    ax3.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax3.spines[sp].set_color(C["rule"])
    ax3.set_title("(c)  No signature\n      falls as far", fontsize=7.0,
                  fontweight="bold", color=C["ink"], loc="left", pad=5,
                  linespacing=1.3)

    fig.text(0.105, 0.045, "Filled marker = transcriptome     Open marker = proteome",
             fontsize=6.4, color=C["muted"])
    return save_all(fig, f"{OUT}/Figure2_composite_and_coherence")


def fig3():
    fig = plt.figure(figsize=(6.85, 3.35)); fig.patch.set_facecolor("white")
    rows = [("panel: secreted plasma proteins", "secreted plasma protein",
             CLS_COLOR["secreted plasma protein"], None),
            ("panel: transcription factors", "transcription factor",
             CLS_COLOR["transcription factor"], None),
            ("panel: metabolic enzymes", "metabolic enzyme",
             CLS_COLOR["metabolic enzyme"], None),
            ("panel: surface receptor (ASGR1)", "surface receptor",
             CLS_COLOR["surface receptor"], None),
            ("liver-made plasma proteins,\noutside the panel", None, C["muted"],
             "off_panel_plasma_mean_z"),
            ("immunoglobulins,\nnot made by the liver", None, C["muted"],
             "immunoglobulin_mean_z"),
            ("erythrocyte proteins", None, C["muted"], "erythrocyte_mean_z")]
    ax = fig.add_axes([0.230, 0.300, 0.300, 0.520])
    for i, (lab, cls, col, key) in enumerate(rows[::-1]):
        for n in DS:
            d = R["datasets"][n]
            v = d["panel_by_class"].get(cls, {}).get("mean_z") if cls else d[key]
            if v is None:
                continue
            tx = LEVEL[n] == "transcriptome"
            ax.plot([v], [i + (0.13 if tx else -0.13)], "o" if tx else "s",
                    ms=4.6, markerfacecolor=col if tx else "white",
                    markeredgecolor=col, markeredgewidth=0.8,
                    alpha=0.55 if tx else 1.0)
    ax.axvline(0, color=C["ink"], lw=0.8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=6.6, linespacing=1.25)
    for t, r in zip(ax.get_yticklabels(), rows[::-1]):
        t.set_color(r[2])
    ax.set_xlabel("mean paired shift (z)", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(C["rule"])
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_title("(a)  Neither blood compartment explains the split",
                 fontsize=7.0, fontweight="bold", color=C["ink"], loc="left", pad=5)
    ax.plot([], [], "o", color=C["muted"], ms=4.6, alpha=0.55,
            label="transcriptome (3 datasets)")
    ax.plot([], [], "s", markerfacecolor="white", markeredgecolor=C["muted"],
            ms=4.6, label="proteome (2 datasets)")
    ax.legend(fontsize=6.4, frameon=False, ncol=2, handletextpad=0.35,
              columnspacing=1.1, loc="upper center", bbox_to_anchor=(0.5, -0.205))

    ns = R["jiang_normalization_sensitivity"]
    order = ("raw", "quantile", "median")
    labs = ["raw", "quantile", "median"]
    x = np.arange(3)

    ax2 = fig.add_axes([0.620, 0.300, 0.125, 0.520])
    # every member, not only the class mean: the claim being made is that the
    # two classes do not overlap, which is a statement about members
    for cls, key in (("secreted plasma protein", "secreted"),
                     ("metabolic enzyme", "metabolic")):
        for xi, k in enumerate(order):
            vals = [ns[k]["members"][g] for g in CLASSES[cls]
                    if g in ns[k]["members"]]
            ax2.scatter([xi + 0.12] * len(vals), vals, s=6,
                        facecolor=CLS_COLOR[cls], alpha=0.45, linewidth=0, zorder=2)
    ax2.plot(x, [ns[k]["secreted_mean_z"] for k in order], "o-",
             color=CLS_COLOR["secreted plasma protein"], ms=5, lw=1.4,
             markeredgecolor="white", markeredgewidth=0.6, label="secreted", zorder=3)
    ax2.plot(x, [ns[k]["metabolic_mean_z"] for k in order], "o-",
             color=CLS_COLOR["metabolic enzyme"], ms=5, lw=1.4,
             markeredgecolor="white", markeredgewidth=0.6, label="metabolic", zorder=3)
    ax2.axhline(0, color=C["ink"], lw=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels(labs, fontsize=6.3, rotation=28, ha="right")
    ax2.set_xlim(-0.45, 2.45); ax2.set_ylim(-1.75, 1.35)
    ax2.set_ylabel("mean paired shift (z)", fontsize=7, labelpad=2)
    ax2.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax2.spines[sp].set_color(C["rule"])
    ax2.text(0.05, ns["raw"]["secreted_mean_z"] + 0.52, "secreted", fontsize=6.2,
             color=CLS_COLOR["secreted plasma protein"], fontweight="bold")
    ax2.text(0.05, ns["raw"]["metabolic_mean_z"] - 0.58, "metabolic", fontsize=6.2,
             color=CLS_COLOR["metabolic enzyme"], fontweight="bold")
    ax2.set_title("(b)  The split is stable\n      across processing",
                  fontsize=7.0, fontweight="bold", color=C["ink"], loc="left",
                  pad=5, linespacing=1.25)

    ax3 = fig.add_axes([0.845, 0.300, 0.125, 0.520])
    vals = [ns[k]["composite_sample_test"]["mean_paired_delta"] for k in order]
    ps = [ns[k]["composite_sample_test"]["wilcoxon_p"] for k in order]
    ax3.bar(x, vals, width=0.55, color=C["hub_f"], edgecolor=C["hub"], linewidth=0.9)
    for xi, v, pv in zip(x, vals, ps):
        ax3.text(xi, v + 0.003, f"{pv:.2g}", ha="center", va="bottom",
                 fontsize=6.0, color=C["muted"])
    ax3.text(0.02, 0.985, "Wilcoxon P", transform=ax3.transAxes, fontsize=5.9,
             color=C["muted"], va="top")
    ax3.axhline(0, color=C["ink"], lw=0.8)
    ax3.set_xticks(x); ax3.set_xticklabels(labs, fontsize=6.3, rotation=28, ha="right")
    ax3.set_xlim(-0.55, 2.55); ax3.set_ylim(0, 0.125)
    ax3.set_ylabel("composite panel score (z)", fontsize=7, labelpad=2)
    ax3.tick_params(labelsize=6.8, length=2.4, pad=1.6)
    for sp in ("top", "right"):
        ax3.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax3.spines[sp].set_color(C["rule"])
    ax3.set_title("(c)  The composite\n      is not",
                  fontsize=7.0, fontweight="bold", color=C["ink"], loc="left",
                  pad=5, linespacing=1.25)
    return save_all(fig, f"{OUT}/Figure3_controls")


for f in (fig1, fig2, fig3):
    print("\n".join("  " + p for p in f()))
