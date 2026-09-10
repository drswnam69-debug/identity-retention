"""Figure 1: what the identity-retention fraction is, and how it is read off
the data. Panel (a) is a schematic. Panel (b) is real GSE14520 data for
the two arms of the worked example, drawn from the archived pairing."""
import importlib.util, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
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
C = {"ink": "#0b0b0b", "muted": "#898781", "rule": "#c3c2b7", "pt": "#8d9bb5",
     "id": "#4a3aa7", "co": "#c99a12", "red": "#199e70", "dra": "#d95926",
     "bad": "#d03b3b"}

fig = plt.figure(figsize=(174 * MM, 78 * MM))
gs = gridspec.GridSpec(1, 2, figure=fig, left=0.055, right=0.965, top=0.88,
                       bottom=0.145, wspace=0.30, width_ratios=[1.0, 1.12])

# ---------------------------------------------------------------- (a) schematic
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(0, 1); ax.set_ylim(-0.16, 1.06); ax.axis("off")
ax.annotate("", xy=(0.97, 0.10), xytext=(0.13, 0.10),
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=C["ink"]))
ax.annotate("", xy=(0.13, 0.95), xytext=(0.13, 0.10),
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=C["ink"]))
ax.text(0.58, -0.155, "more non-parenchymal content\n(what purity adjustment targets)",
        ha="center", va="bottom", fontsize=7.2, color=C["co"])
ax.text(0.045, 0.52, "more loss of tissue identity\n(what the covariate stands for)",
        ha="center", va="center", fontsize=7.2, color=C["id"], rotation=90)
ax.scatter([0.13], [0.10], s=34, color=C["ink"], zorder=5)
ax.text(0.115, 0.055, "adjacent\ntissue", fontsize=7.2, color=C["ink"],
        ha="right", va="top", linespacing=1.25)
paths = [((0.86, 0.14), C["co"], "purity change only",
          "identity adjustment leaves it intact"),
         ((0.20, 0.80), C["id"], "dedifferentiation only",
          "identity adjustment removes it"),
         ((0.82, 0.52), C["bad"], "a real tumor", "both, in an unknown ratio")]
for (x, y), col, lab, sub in paths:
    ax.add_patch(FancyArrowPatch((0.13, 0.10), (x, y), arrowstyle="-|>",
                                 mutation_scale=9, lw=1.7, color=col,
                                 connectionstyle="arc3,rad=0.0", zorder=4))
    ax.scatter([x], [y], s=30, color=col, zorder=5, edgecolor="white", lw=0.6)
va = "bottom"
ax.text(0.905, 0.155, "purity change only\nthe identity covariate\ndoes not move",
        fontsize=7.2, color=C["co"], ha="right", va="bottom", linespacing=1.3)
ax.text(0.245, 0.955, "dedifferentiation only\nthe identity covariate\nmoves with it",
        fontsize=7.2, color=C["id"], ha="left", va="top", linespacing=1.3)
ax.text(0.845, 0.52, "a real tumor:\nboth, in an\nunknown ratio",
        fontsize=7.2, color=C["bad"], ha="left", va="center", linespacing=1.3)
ax.add_patch(FancyArrowPatch((0.13, 0.10), (0.42, 0.52), arrowstyle="-|>",
                             mutation_scale=8, lw=1.5, color=C["muted"],
                             linestyle=(0, (3.2, 2.0)), zorder=3))
ax.text(0.235, 0.62, "a random covariate that falls as far\nremoves as much (\u00a76ae)",
        fontsize=7.2, color=C["muted"], ha="left", va="bottom", linespacing=1.3)
ax.text(0.13, 1.02, "The two axes are not the same axis, but the vertical one is\n"
        "the design of the measure rather than a demonstrated cause",
        fontsize=7.2, color=C["muted"], ha="left", va="bottom", linespacing=1.3)

# ---------------------------------------------------------------- (b) real data
sys.path.insert(0, f"{IR_ROOT}/code")
os.chdir(f"{IR_ROOT}")
def L(n, f):
    s = importlib.util.spec_from_file_location(n, f"code/{f}")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
_da = L("d", "12_differentiation_adjust.py"); _bm = L("b", "30_signature_benchmark.py")
import rsi_config as cfg

ols_ci, D1, zmean = _da.ols_ci, _da.D1, _bm.zmean
rsi = pd.read_csv("results/GSE14520/rsi.tsv", sep="\t", index_col=0)
ph = pd.read_csv("results/GSE14520/phenotype.tsv", sep="\t", index_col=0).loc[rsi.index]
expr = pd.read_csv("data/GSE14520_symbols.tsv.gz", sep="\t", index_col=0)[rsi.index]
is_t = (~ph["tissue"].astype(str).str.lower().str.contains("adjacent")).to_numpy()
fr = pd.DataFrame({"pid": ph["patient_id"].astype(str).values, "t": is_t}, index=rsi.index)
def paired(s):
    f = fr.assign(v=s.to_numpy())
    tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
    nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
    k = sorted(set(tt.index) & set(nn.index))
    return (tt.loc[k] - nn.loc[k]).to_numpy(float)
dD1 = paired(zmean(expr, D1)[0])
os.chdir(IR_DOCS)
ax = fig.add_subplot(gs[0, 1])
NOUT = 0
X = np.column_stack([np.ones(len(dD1)), dD1])
xs = np.linspace(-3.1, 0.75, 40)
for mod, col, lab in (("DRAIN", C["dra"], "DRAIN"), ("REDUCTION", C["red"], "REDUCTION")):
    os.chdir(f"{IR_ROOT}")
    dv = paired(zmean(expr, getattr(cfg, f"MODULE_{mod}"))[0])
    os.chdir(IR_DOCS)
    fit = ols_ci(dv, X, ["intercept", "dD1"])
    a0 = fit["intercept"]["beta"]; b1 = fit["dD1"]["beta"]
    NOUT += int((np.abs(dv) > 2.45).sum())
    ax.scatter(dD1, dv, s=8, color=col, alpha=0.32, edgecolor="none", zorder=3)
    ax.plot(xs, a0 + b1 * xs, lw=1.5, color=col, zorder=4)
    ax.scatter([0], [a0], s=48, marker="D", color=col, edgecolor="white", lw=0.8,
               zorder=6)
    un = dv.mean()
    ax.annotate("", xy=(0.62, un), xytext=(0.62, 0), zorder=5,
                arrowprops=dict(arrowstyle="<->", lw=0.9, color=col, shrinkA=0, shrinkB=0))
    ax.text(0.72, un / 2, f"unadjusted\n{un:+.2f}".replace("-", "\u2212"),
            fontsize=7.2, color=col,
            va="center", ha="left", linespacing=1.25)
    # The falling arm's label sat on top of its own regression line, which
    # struck the text through. Place it clear of the line and give both labels
    # a light box so neither can be overprinted.
    ly = 2.05 if a0 > 0 else -1.95
    lx = -2.45 if a0 > 0 else -0.95
    ax.annotate(f"{lab}\nintercept {a0:+.2f}\nretains {abs(a0)/abs(un)*100:.0f}%".replace("-", "\u2212"),
                xy=(0, a0), xytext=(lx, ly), fontsize=7.2, color=col,
                ha="center", va="center", fontweight="bold", linespacing=1.35,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none",
                          alpha=0.85), zorder=7)
ax.axhline(0, lw=0.7, color=C["ink"])
ax.axvline(0, ls=(0, (3, 2)), lw=0.9, color=C["ink"], zorder=2)
ax.text(0.12, -2.42, "no identity lost", fontsize=7.2, color=C["ink"],
        ha="left", va="bottom")
ax.set_xlim(-3.25, 2.35)
ax.set_ylim(-2.45, 2.45)
ax.set_xlabel("ΔD1, paired change in hepatocyte identity")
ax.text(2.30, -1.15, f"axes trimmed at −3.25 and at \u00b12.45;\n"
        f"{int((dD1 < -3.25).sum())} pairs lie to the left and {NOUT} points\n"
        f"above or below, all of them in every fit",
        fontsize=7.2, color=C["muted"], ha="right", va="top", linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.85))
ax.set_ylabel("Δ module score, tumor minus adjacent (z units)")
ax.set_title(f"The intercept as a share of the shift (GSE14520, {len(dD1)} pairs)",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.text(0.004, 0.985, "a", fontsize=10, fontweight="bold", ha="left", va="top")
fig.text(0.487, 0.985, "b", fontsize=10, fontweight="bold", ha="left", va="top")

stem = f"{IR_DOCS}/Figure1_concept"
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png").convert("RGB")
im.save(stem + ".png", dpi=(300, 300))
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size)
