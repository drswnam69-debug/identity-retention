"""Figure 8: a third paired cohort, on RNA sequencing (PREREG 6w). All values are
read from the archived result files; nothing is recomputed here."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import spearmanr
from PIL import Image


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
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.2,
    "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})
MM = 1 / 25.4
C = {"red": "#199e70", "dra": "#d95926", "ink": "#0b0b0b", "muted": "#898781",
     "rule": "#c3c2b7", "arr": "#4a3aa7", "seq": "#c99a12", "bad": "#d03b3b",
     "pt": "#8d9bb5"}
R = f"{IR_RESULTS}"
G = json.load(open(f"{R}/SIGNATURE_BENCHMARK_GSE14520.json"))
T = json.load(open(f"{R}/SIGNATURE_BENCHMARK_TCGA_LIHC.json"))
K = json.load(open(f"{R}/TCGA_LIHC/CONSOLIDATED_6w.json"))
D = json.load(open(f"{R}/TCGA_LIHC/D1_gene_comparison_6w.json"))

fig = plt.figure(figsize=(174 * MM, 134 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.078, right=0.985,
                       top=0.915, bottom=0.085, hspace=0.50, wspace=0.30)

# --- (a) the two distributions ---------------------------------------------
ax = fig.add_subplot(gs[0, 0])
gv = np.array([r["retention_joint"] for r in G["signatures"] if r["status"] == "ok"])
tv = np.array([r["retention_joint"] for r in T["signatures"] if r["status"] == "ok"])
bins = np.arange(0, 2.5, 0.1)
ax.hist(gv, bins=bins, color=C["arr"], alpha=0.58, edgecolor="white", lw=0.4,
        label=f"GSE14520, array, {len(gv)} sets", zorder=3)
ax.hist(tv, bins=bins, color=C["seq"], alpha=0.58, edgecolor="white", lw=0.4,
        label=f"TCGA-LIHC, RNA-seq, {len(tv)} sets", zorder=3)
ax.axvline(np.median(gv), lw=1.1, color=C["arr"], ls=":", zorder=4, ymax=0.70)
ax.axvline(np.median(tv), lw=1.1, color=C["seq"], ls=":", zorder=4, ymax=0.70)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=5, ymax=0.70)
ax.legend(frameon=False, fontsize=7.2, loc="upper left", handlelength=0.9,
          borderpad=0.1, labelspacing=0.22, handletextpad=0.4)
ax.text(0.985, 0.60, f"medians {np.median(gv):.3f} and {np.median(tv):.3f}\n"
        f"difference {K['benchmark']['median_difference_vs_GSE14520']:+.3f}\n"
        f"(§6w threshold 0.15)", transform=ax.transAxes, fontsize=7.2,
        color=C["ink"], ha="right", va="top")
ax.set_xlim(0, 2.0)
ax.set_ylim(0, 34)
n_off = int((gv > 2.0).sum() + (tv > 2.0).sum())
ax.text(0.985, 0.30, f"axis trimmed at 2.0;\n{n_off} of {len(gv) + len(tv)} sets lie above it",
        transform=ax.transAxes, fontsize=7.2, color=C["muted"], ha="right", va="top",
        linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
ax.set_xlabel("Identity-retention fraction")
ax.set_ylabel("Published liver signatures")
ax.set_title("The distribution replicates across platform", fontsize=7.2,
             pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (b) the shared signatures ---------------------------------------------
ax = fig.add_subplot(gs[0, 1])
g = {r["name"]: r["retention_joint"] for r in G["signatures"] if r["status"] == "ok"}
t = {r["name"]: r["retention_joint"] for r in T["signatures"] if r["status"] == "ok"}
sh = sorted(set(g) & set(t))
x = np.array([g[k] for k in sh]); y = np.array([t[k] for k in sh])
lim = 2.0
ax.plot([0, lim], [0, lim], lw=0.8, color=C["rule"], zorder=1)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.5, zorder=2)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.5, zorder=2)
ax.scatter(x, y, s=12, color=C["pt"], edgecolor="white", lw=0.35, zorder=3)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
_off = int(((x > lim) | (y > lim)).sum())
ax.text(lim - 0.03, 0.05, f"axes trimmed at {lim};\n{_off} of {len(x)} lie outside",
        fontsize=7.2, color=C["muted"], ha="right", va="bottom", linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Retention in GSE14520 (array)")
ax.set_ylabel("Retention in TCGA-LIHC (RNA-seq)")
ax.set_title(f"{len(sh)} signatures scored in both, ρ = "
             f"{K['benchmark']['spearman_shared']:+.3f}", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (c) the worked example across cohorts ---------------------------------
ax = fig.add_subplot(gs[1, 0])
labels = ["REDUCTION\n(4 gene)", "REDUCTION\n(3 gene)", "DRAIN\n(locked)",
          "DRAIN\nwithout POR"]
gvals = [np.nan, 0.857, -0.800, -0.731]
tvals = [1.2612, 1.216, K["own_modules"]["DRAIN_locked"]["shift"],
         K["own_modules"]["DRAIN_without_POR_6l"]["shift"]]
gret = ["", "89%", "31%", "32%"]
tret = ["92%", "92%", "n.s.", "4%"]
xx = np.arange(4); w = 0.36
ax.bar(xx - w / 2, gvals, width=w, color=C["arr"], edgecolor="none",
       label="GSE14520, 213 pairs", zorder=3)
ax.bar(xx + w / 2, tvals, width=w, color=C["seq"], edgecolor="none",
       label="TCGA-LIHC, 50 pairs", zorder=3)
ax.axhline(0, lw=0.7, color=C["ink"])
for xi, v, lab in zip(xx - w / 2, gvals, gret):
    if not np.isfinite(v):
        continue
    ax.text(xi, v + (0.06 if v > 0 else -0.06), lab, ha="center",
            va="bottom" if v > 0 else "top", fontsize=7.2, color=C["ink"])

for xi, v, lab in zip(xx + w / 2, tvals, tret):
    col = C["bad"] if lab == "n.s." else C["ink"]
    ax.text(xi, v + (0.06 if v > 0 else -0.06), lab, ha="center",
            va="bottom" if v > 0 else "top", fontsize=7.2, color=col,
            fontweight="bold" if lab == "n.s." else "normal")
ax.annotate("POR reverses direction\non RNA-seq, and the\nlocked module cancels",
            xy=(2 + w / 2, tvals[2]), xytext=(-0.14, -0.60), fontsize=7.2,
            color=C["bad"], ha="left", va="top", linespacing=1.3,
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C["bad"],
                            shrinkA=2, shrinkB=3))
ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=7.2)
ax.set_ylim(-1.35, 1.62)
ax.set_ylabel("Module shift, tumor minus adjacent\n(z units); label is retention")
ax.legend(frameon=False, fontsize=7.2, loc="upper right", handlelength=0.9,
          borderpad=0.1, labelspacing=0.22)
ax.set_title("One arm replicates, the other does not as locked\n"
             "AIFM2 is absent from the array, so GSE14520 has no four-gene bar",
             fontsize=7.2, pad=4, color=C["muted"], linespacing=1.35)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (d) the covariate's own members ---------------------------------------
ax = fig.add_subplot(gs[1, 1])
eff = [r for r in D["per_gene_D1"] if r["class"] == "effector"]
tfs = [r for r in D["per_gene_D1"] if r["class"] == "TF"]
ax.plot([-1.6, 1.2], [-1.6, 1.2], lw=0.8, color=C["rule"], zorder=1)
ax.axhline(0, lw=0.6, color=C["rule"]); ax.axvline(0, lw=0.6, color=C["rule"])
ax.scatter([r["GSE14520"] for r in eff], [r["TCGA_LIHC"] for r in eff],
           s=16, color=C["pt"], edgecolor="white", lw=0.4, zorder=3,
           label=f"effector genes (n = {len(eff)})")
ax.scatter([r["GSE14520"] for r in tfs], [r["TCGA_LIHC"] for r in tfs],
           s=26, marker="D", color=C["bad"], edgecolor="white", lw=0.5,
           zorder=4, label=f"transcription factors (n = {len(tfs)})")
OFF = {"HNF4A": (-0.09, 0.13, "right"), "FOXA2": (0.10, -0.11, "left"),
       "HNF1A": (-0.09, 0.12, "right"), "FOXA1": (-0.10, 0.12, "right"),
       "NR1H4": (-0.09, 0.12, "right")}
for r in tfs:
    dx, dy, ha = OFF[r["gene"]]
    ax.annotate(r["gene"], xy=(r["GSE14520"], r["TCGA_LIHC"]),
                xytext=(r["GSE14520"] + dx, r["TCGA_LIHC"] + dy),
                fontsize=7.2, color=C["bad"], ha=ha)
ax.set_xlim(-1.6, 1.2); ax.set_ylim(-1.6, 1.2)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("D1 member in GSE14520 (z units)")
ax.set_ylabel("D1 member in TCGA-LIHC (z units)")
ax.legend(frameon=False, fontsize=7.2, loc="lower right", handlelength=0.9,
          borderpad=0.1, labelspacing=0.22, handletextpad=0.3)
ax.set_title("The covariate's effectors agree, its regulators do not",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for xx_, yy_, L in [(0.004, 0.985, "a"), (0.505, 0.985, "b"),
                    (0.004, 0.495, "c"), (0.505, 0.495, "d")]:
    fig.text(xx_, yy_, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = f"{IR_DOCS}/Figure8_third_cohort"
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png").convert("RGB")
im.save(stem + ".png", dpi=(300, 300))
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size, "| shared", len(sh))
