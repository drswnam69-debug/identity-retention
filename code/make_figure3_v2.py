"""Figure 3 for the v2 master manuscript.

Panel order follows the v2 argument: the module dissociation leads, the
composite index follows. Every value is read from the archived result JSONs
(GSE76427_differentiation / _composition, GSE14520_differentiation /
_composition) and none is recomputed here.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PIL import Image

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 7.2,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
})

MM = 1 / 25.4
WIDTH_MM = 174

C = {
    "red": "#199e70",     # REDUCTION, CoQ-reducing arm
    "red_f": "#d7f0e6",
    "dra": "#d95926",     # DRAIN, N-reductive arm
    "dra_f": "#fbe4d8",
    "ink": "#0b0b0b",
    "muted": "#898781",
    "rule": "#c3c2b7",
}

# Read from the archive, not from an upload directory. This script used to
# read its input from outside the deposited tree, so the figure could not be
# regenerated from the deposit alone.
D = json.load(open("/home/claude/rsi/results/FIG4_INPUTS.json"))
COH = ["GSE76427", "GSE14520"]
LAB = {"GSE76427": "GSE76427\n(52 pairs)", "GSE14520": "GSE14520\n(213 pairs)"}


def panel_label(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=10, fontweight="bold",
             color=C["ink"], ha="left", va="top")


fig = plt.figure(figsize=(WIDTH_MM * MM, 132 * MM))
gs = gridspec.GridSpec(
    2, 2, figure=fig,
    left=0.145, right=0.985, top=0.915, bottom=0.085,
    hspace=0.55, wspace=0.26,
    width_ratios=[1.0, 1.30], height_ratios=[1.0, 1.05])

# ----------------------------------------------------------------- panel (a)
# Identity-retention fraction, joint D1 + C1, both modules, both cohorts.
axa = fig.add_subplot(gs[0, 0])
rows, colors, labels = [], [], []
for coh in COH:
    rows.append(D[coh]["module_reduction_joint"]["ret"] * 100)
    colors.append(C["red"])
    labels.append("REDUCTION")
    rows.append(abs(D[coh]["module_drain_joint"]["ret"]) * 100)
    colors.append(C["dra"])
    labels.append("DRAIN")

ypos = [3.4, 2.6, 1.4, 0.6]
axa.barh(ypos, rows, height=0.62, color=colors, edgecolor="none", zorder=3)
axa.axvline(50, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=4)
for y, v in zip(ypos, rows):
    axa.text(v + 1.8, y, f"{v:.1f}%", va="center", ha="left",
             fontsize=7.2, color=C["ink"], zorder=5)
axa.set_yticks(ypos)
axa.set_yticklabels(labels, fontsize=7.2)
axa.set_xlim(0, 108)
axa.set_ylim(0.05, 4.45)
axa.set_xlabel("Identity-retention fraction, joint model (%)")
axa.text(50, 4.30, "50% threshold", fontsize=7.2, color=C["ink"],
         ha="center", va="bottom")
axa.text(2, 3.90, "GSE76427, 52 pairs", fontsize=7.2, color=C["muted"],
         ha="left", va="center")
axa.text(2, 1.90, "GSE14520, 213 pairs", fontsize=7.2, color=C["muted"],
         ha="left", va="center")
axa.axhline(2.0, lw=0.6, color=C["rule"])
for s in ("top", "right"):
    axa.spines[s].set_visible(False)

# ----------------------------------------------------------------- panel (b)
# Absolute module shift at three successive stages.
gsb = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 1], wspace=0.20)
stages = ["Unadjusted", "D1", "D1 + C1"]
for j, coh in enumerate(COH):
    ax = fig.add_subplot(gsb[0, j])
    d = D[coh]
    series = [
        ("REDUCTION", C["red"], [
            (d["module_reduction"]["unadj"], None),
            (d["module_reduction"]["adj"], d["module_reduction"]["ci"]),
            (d["module_reduction_joint"]["est"], d["module_reduction_joint"]["ci"]),
        ]),
        ("DRAIN", C["dra"], [
            (d["module_drain"]["unadj"], None),
            (d["module_drain"]["adj"], d["module_drain"]["ci"]),
            (d["module_drain_joint"]["est"], d["module_drain_joint"]["ci"]),
        ]),
    ]
    x = [0, 1, 2]
    for name, col, pts in series:
        ys = [p[0] for p in pts]
        ax.plot(x, ys, "-", lw=1.0, color=col, zorder=3)
        for xi, (yi, ci) in zip(x, pts):
            if ci is None:
                ax.plot([xi], [yi], "o", ms=4.2, mfc="white", mec=col,
                        mew=1.1, zorder=5)
            else:
                ax.plot([xi, xi], ci, "-", lw=1.0, color=col, zorder=4)
                ax.plot([xi], [yi], "o", ms=4.2, color=col, zorder=5)
    ax.axhline(0, lw=0.6, color=C["ink"])
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=7.2)
    ax.set_xlim(-0.35, 2.35)
    ax.set_ylim(-1.15, 1.15)
    ax.set_title(LAB[coh].replace("\n", " "), fontsize=7.2, pad=3,
                 color=C["muted"])
    if j == 0:
        ax.set_ylabel("Module shift, tumor minus adjacent\n(z units)")
    else:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if j == 0:
        ax.text(-0.25, 1.06, "REDUCTION", fontsize=7.2, color=C["red"],
                ha="left", va="center", fontweight="bold")
        ax.text(-0.25, -1.06, "DRAIN", fontsize=7.2, color=C["dra"],
                ha="left", va="center", fontweight="bold")

# ----------------------------------------------------------------- panel (c)
# Composite index across the full adjustment cascade.
gsc = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.10)
models = [
    ("Unadjusted", "unadj"),
    ("D1", "composite_D1"),
    ("D2 (over-adj.)", "composite_D2"),
    ("C1 alone", "composite_C1"),
    ("D1 + C1", "composite_joint"),
    ("D1 + C2", "composite_joint_C2"),
]
for j, coh in enumerate(COH):
    ax = fig.add_subplot(gsc[0, j])
    d = D[coh]
    ys = list(range(len(models)))[::-1]
    for y, (name, key) in zip(ys, models):
        if key == "unadj":
            est, ci = d["composite_unadj_median"], None
        elif key in ("composite_D1", "composite_D2"):
            est, ci = d[key]["beta"], d[key]["ci95"]
        else:
            est, ci = d[key]["est"], d[key]["ci"]
        emph = key in ("composite_joint",)
        col = C["ink"] if emph else C["muted"]
        if ci is None:
            ax.plot([est], [y], "o", ms=4.6, mfc="white", mec=C["ink"],
                    mew=1.1, zorder=5)
        else:
            ax.plot(ci, [y, y], "-", lw=1.4 if emph else 1.0, color=col,
                    zorder=4)
            ax.plot([est], [y], "o", ms=4.6 if emph else 4.0, color=col,
                    zorder=5)
    ax.axvline(0, lw=0.6, color=C["ink"])
    ax.set_yticks(ys)
    ax.set_yticklabels([m[0] for m in models], fontsize=7.2)
    ax.set_ylim(-0.6, len(models) - 0.4)
    ax.set_xlim(-0.15, 1.75)
    ax.set_xlabel("Composite index shift at zero covariate change (z units)")
    ax.set_title(LAB[coh].replace("\n", " "), fontsize=7.2, pad=3,
                 color=C["muted"])
    if j == 1:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

panel_label(fig, 0.006, 0.985, "a")
panel_label(fig, 0.425, 0.985, "b")
panel_label(fig, 0.006, 0.495, "c")

stem = "/home/claude/Figure3_dissociation"
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png")
if im.mode != "RGB":
    bg = Image.new("RGB", im.size, "white")
    bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
    im = bg
im.save(stem + ".png", dpi=(300, 300))
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size)
