"""Figure 5 for v4: the identity-retention measure does not transfer between
transcriptome and proteome (PREREG 6s). Every value is read from the archived
result files or from the verified per-gene diagnostics; nothing is recomputed.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
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
     "rule": "#c3c2b7", "tx": "#4a3aa7", "pr": "#256abf", "bad": "#d03b3b"}

G = json.load(open(f"{IR_RESULTS}/PROTEIN_VALIDATION_6s_Gao2019.json"))
J = json.load(open(f"{IR_RESULTS}/PROTEIN_VALIDATION_6s_Jiang2019.json"))

SETS = ["GSE14520\ntranscript\n213 pairs", "Gao 2019\nprotein\n159 pairs",
        "Jiang 2019\nprotein\n124 pairs"]

# Liberation Sans has no superscript minus (U+207B) or superscript zero
# (U+2070); written as literal characters they print as empty boxes, which is
# how "10 to the minus ten" came to render as "10 1 0". Mathtext renders the
# exponent from matplotlib's own math font, which carries both.
CYB5R3 = [(+0.555, r"$9.1 \times 10^{-10}$"),
          (-0.387, r"$3.0 \times 10^{-4}$"),
          (+0.521, r"$1.2 \times 10^{-6}$")]
RED = [+0.858, G["modules"]["REDUCTION"]["unadjusted_mean_delta"],
       J["modules"]["REDUCTION"]["unadjusted_mean_delta"]]
DRA = [-0.800, G["modules"]["DRAIN"]["unadjusted_mean_delta"],
       J["modules"]["DRAIN"]["unadjusted_mean_delta"]]
RED_P = ["8.9e-31", "0.36", "1.9e-12"]
DRA_P = ["1.9e-24", "2.3e-21", "8.3e-16"]

fig = plt.figure(figsize=(174 * MM, 144 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.095, right=0.985,
                       top=0.92, bottom=0.10, hspace=0.55, wspace=0.34)

# --- (a) CYB5R3 direction reversal --------------------------------------
ax = fig.add_subplot(gs[0, 0])
x = np.arange(3)
vals = [v for v, _ in CYB5R3]
cols = [C["tx"], C["bad"], C["pr"]]
ax.bar(x, vals, width=0.55, color=cols, edgecolor="none", zorder=3)
ax.axhline(0, lw=0.7, color=C["ink"])
for xi, (v, p) in zip(x, CYB5R3):
    off = 0.04 if v > 0 else -0.04
    ax.text(xi, v + off, f"{v:+.3f}\nP = {p}", ha="center",
            va="bottom" if v > 0 else "top", fontsize=7.2, color=C["ink"])
ax.set_xticks(x)
ax.set_xticklabels(SETS, fontsize=7.2)
ax.set_ylim(-0.72, 0.88)
ax.set_ylabel("CYB5R3, tumor minus adjacent\n(z units)")
ax.set_title("The donor to both arms reverses direction", fontsize=7.2,
             pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (b) the two modules across datasets ---------------------------------
ax = fig.add_subplot(gs[0, 1])
w = 0.36
ax.bar(x - w / 2, RED, width=w, color=C["red"], edgecolor="none",
       label="REDUCTION", zorder=3)
ax.bar(x + w / 2, DRA, width=w, color=C["dra"], edgecolor="none",
       label="DRAIN", zorder=3)
ax.axhline(0, lw=0.7, color=C["ink"])
for xi, v, p in zip(x - w / 2, RED, RED_P):
    ax.text(xi, v + (0.04 if v > 0 else -0.04), f"{v:+.2f}", ha="center",
            va="bottom" if v > 0 else "top", fontsize=7.2, color=C["ink"])
for xi, v, p in zip(x + w / 2, DRA, DRA_P):
    ax.text(xi, v - 0.04, f"{v:+.2f}", ha="center", va="top",
            fontsize=7.2, color=C["ink"])
ax.annotate("no shift\n(P = 0.36)", xy=(1 - w / 2, RED[1]),
            xytext=(1 - w / 2 - 0.42, 0.52), fontsize=7.2, color=C["bad"],
            ha="center",
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C["bad"],
                            shrinkA=1, shrinkB=2))
ax.set_xticks(x)
ax.set_xticklabels(SETS, fontsize=7.2)
ax.set_ylim(-1.25, 1.15)
ax.set_ylabel("Module shift, tumor minus adjacent\n(z units)")
ax.legend(frameon=False, fontsize=7.2, loc="upper right",
          handlelength=1.1, borderpad=0.2)
ax.set_title("The N-reductive fall is the robust result in these three", fontsize=7.2,
             pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (c) DRAIN retention where the premise holds -------------------------
ax = fig.add_subplot(gs[1, 0])
labels = ["Transcript\nD1 full\n(22 genes)",
          "Transcript\nD1 restricted\n(19 measured)",
          "Gao protein\nD1\n(19 genes)"]
ret = [24.7, 22.0, G["modules"]["DRAIN"]["D1"]["retention"] * 100]
cols = [C["tx"], C["tx"], C["bad"]]
ax.bar(np.arange(3), ret, width=0.55, color=cols, edgecolor="none", zorder=3)
ax.axhline(50, ls=(0, (3, 2)), lw=0.8, color=C["ink"])
for xi, v in zip(np.arange(3), ret):
    ax.text(xi, v + 1.8, f"{v:.1f}%", ha="center", va="bottom",
            fontsize=7.2, color=C["ink"])
ax.text(0.55, 53.5, "50% threshold", fontsize=7.2, color=C["ink"],
        ha="right", va="bottom")
ax.set_xticks(np.arange(3))
ax.set_xticklabels(labels, fontsize=7.2)
ax.set_ylim(0, 104)
ax.set_ylabel("DRAIN identity-retention fraction (%)")
ax.set_title("Restricting the covariate does not explain the difference",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (d) why the second proteome is excluded -----------------------------
ax = fig.add_subplot(gs[1, 1])
d1 = [("ARG1", -0.883), ("OTC", -0.608), ("CPS1", -0.802),
      ("FGA", +0.425), ("HNF4A", +0.429), ("ALB", +0.731),
      ("SERPINA1", +0.797)]
d1 = sorted(d1, key=lambda t: t[1])
y = np.arange(len(d1))
vals = [v for _, v in d1]
TF = {"HNF4A"}
cols = [C["muted"] if v < 0 else ("#7b4ea8" if n in TF else C["bad"])
        for n, v in d1]
ax.barh(y, vals, height=0.62, color=cols, edgecolor="none", zorder=3)
ax.axvline(0, lw=0.7, color=C["ink"])
for yi, v in zip(y, vals):
    ax.text(v + (0.03 if v > 0 else -0.03), yi, f"{v:+.3f}",
            va="center", ha="left" if v > 0 else "right",
            fontsize=7.2, color=C["ink"])
ax.set_yticks(y)
ax.set_yticklabels([n for n, _ in d1], fontsize=7.2)
ax.set_xlim(-1.55, 1.55)
ax.set_xlabel("Jiang 2019, tumor minus adjacent (z units)")
ax.text(1.50, 0.35, "urea-cycle\nenzymes fall", fontsize=7.2,
        color=C["muted"], ha="right", va="center", linespacing=1.3)
ax.text(-1.50, 5.75, "secreted plasma\nproteins rise", fontsize=7.2,
        color=C["bad"], ha="left", va="center", linespacing=1.3)
ax.text(-1.50, 4.10, "and so does the\ntranscription factor", fontsize=7.2,
        color="#7b4ea8", ha="left", va="center", linespacing=1.3)
ax.set_title("The hepatocyte-identity premise fails in this matrix",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for xx, yy, L in [(0.004, 0.985, "a"), (0.505, 0.985, "b"),
                  (0.004, 0.495, "c"), (0.505, 0.495, "d")]:
    fig.text(xx, yy, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = f"{IR_DOCS}/Figure5_level_transfer"
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
