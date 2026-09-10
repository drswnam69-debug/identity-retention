"""Figure 12: what the covariate contributes, and what its size contributes (PREREG 6ae).

The paper's most consequential control had no display item. Three panels:
(a) the median retention each pool of covariates leaves, against D1;
(b) the mechanism, one point per random draw: how far the covariate itself moves
    against how much it removes;
(c) per-signature agreement between D1 and the matched random covariates.
"""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np, pandas as pd
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
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.2,
    "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5})
MM = 1 / 25.4
C = {"ink": "#0b0b0b", "muted": "#898781", "pt": "#8d9bb5", "id": "#4a3aa7",
     "co": "#c99a12", "red": "#199e70", "dra": "#d95926", "bad": "#d03b3b"}
R = f"{IR_ROOT}"
NC = json.load(open(f"{R}/results/NEGATIVE_CONTROL_6ae.json"))
REF = json.load(open(f"{R}/results/NEGATIVE_CONTROL_6ae_reference.json"))
D = pd.read_csv(f"{R}/results/NEGATIVE_CONTROL_6ae_draws.csv")
# S9 ships at the archive root; fall back to the author's build tree so the
# same script serves both.
_s9 = f"{IR_ROOT}/SupplementaryTable_S9_negative_control.csv"
if not _os.path.exists(_s9):
    _s9 = (f"{IR_DOCS}/SF2build/SupplementaryFile2/"
           "SupplementaryTable_S9_negative_control.csv")
S9 = pd.read_csv(_s9)
nc = NC["negative_control"]

fig = plt.figure(figsize=(174 * MM, 62 * MM))
gs = gridspec.GridSpec(1, 3, figure=fig, left=0.132, right=0.982, top=0.83,
                       bottom=0.215, wspace=0.50)

# --- (a) what each pool of covariates leaves --------------------------------
ax = fig.add_subplot(gs[0, 0])
POOLS = [("no_tumor_shift", "random,\nno tumor shift"),
         ("unrestricted", "random,\nunrestricted"),
         ("matched", "random, matched\nto D1's fall"),
         ("rises_in_tumor", "random,\nrises in tumor")]
ys = np.arange(len(POOLS))
for i, (k, lab) in enumerate(POOLS):
    v = D[D["pool"] == k]["median_retention"].to_numpy()
    lo, hi = np.percentile(v, [2.5, 97.5])
    ax.plot([lo, hi], [i, i], lw=2.6, color=C["pt"], solid_capstyle="butt", zorder=3)
    ax.scatter([np.median(v)], [i], s=26, color=C["ink"], zorder=4)
ax.axvline(REF["D1"]["median"], lw=1.2, color=C["id"], zorder=5)
ax.text(REF["D1"]["median"] + 0.012, len(POOLS) - 0.42,
        f"D1 as locked\n{REF['D1']['median']:.3f}", color=C["id"], fontsize=7.2,
        ha="left", va="top", linespacing=1.3, fontweight="bold")
ax.axvline(1.0, ls=(0, (1, 2)), lw=0.7, color=C["muted"], zorder=2)
ax.set_yticks(ys)
ax.set_yticklabels([l for _, l in POOLS], fontsize=7.2)
ax.set_ylim(-0.7, len(POOLS) - 0.3)
ax.set_xlim(0.25, 1.12)
ax.set_xlabel("Median retention across the 119 signatures")
# The exact count says more than the rounded percentage, and says it in the
# direction the text now takes: D1 sits at the bottom edge of its own null.
_n_ge = int(round((1 - nc["fraction_random_at_or_below_D1_median"]) * nc["n_draws"]))
ax.set_title(f"Only {_n_ge} of {nc['n_draws']} matched random covariates\n"
             f"leave as much standing as D1 does",
             fontsize=7.2, pad=4, color=C["muted"], linespacing=1.35)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

# --- (b) the mechanism ------------------------------------------------------
ax = axb = fig.add_subplot(gs[0, 1])
COL = {"no_tumor_shift": C["pt"], "unrestricted": C["muted"],
       "matched": C["co"], "rises_in_tumor": C["red"]}
for k, _ in POOLS:
    d = D[D["pool"] == k]
    ax.scatter(d["covariate_mean_paired_shift"], d["median_retention"], s=8,
               color=COL[k], alpha=0.55, edgecolor="none", zorder=3)
ax.scatter([REF["D1"]["shift"]], [REF["D1"]["median"]], s=52, marker="D",
           color=C["id"], edgecolor="white", lw=0.7, zorder=6)
ax.text(REF["D1"]["shift"] + 0.05, REF["D1"]["median"] + 0.035, "D1",
        color=C["id"], fontsize=7.2, fontweight="bold", ha="left", va="bottom")
ax.scatter([REF["C1_alone"]["shift"]], [REF["C1_alone"]["median"]], s=52, marker="D",
           color=C["co"], edgecolor="white", lw=0.7, zorder=6)
ax.text(REF["C1_alone"]["shift"] - 0.05, REF["C1_alone"]["median"] + 0.03, "C1",
        color=C["co"], fontsize=7.2, fontweight="bold", ha="right", va="bottom")
ax.axhline(1.0, ls=(0, (1, 2)), lw=0.7, color=C["muted"], zorder=2)
ax.axvline(0.0, lw=0.7, color=C["ink"], zorder=2)
ax.set_xlabel("The covariate's own paired shift (z units)")
ax.set_ylabel("Median retention it leaves")
ax.set_title("What a covariate removes is set by how far it moves",
             fontsize=7.2, pad=4, color=C["muted"])
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

# --- (c) per-signature agreement --------------------------------------------
ax = fig.add_subplot(gs[0, 2])
x = S9["retention_joint_D1"].to_numpy()
y = S9["retention_random_matched_median"].to_numpy()
ax.scatter(x, y, s=11, color=C["pt"], edgecolor="white", lw=0.35, zorder=3)
lim = (0, 2.0)
ax.plot(lim, lim, lw=0.8, color=C["muted"], zorder=2)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], zorder=2)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], zorder=2)
n_out = int(((x > lim[1]) | (y > lim[1])).sum())
ax.set_xlim(*lim); ax.set_ylim(*lim)
ax.text(1.97, 0.06, f"axes trimmed at 2.0;\n{n_out} signatures lie outside",
        fontsize=7.2, color=C["muted"], ha="right", va="bottom", linespacing=1.3,
        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.82))
ax.text(0.05, 1.94, f"median $r$ = {nc['median_r_with_D1_retention']:.3f}",
        fontsize=7.2, color=C["ink"], ha="left", va="top")
ax.set_xlabel("Retention under D1")
ax.set_ylabel("Retention under a matched\nrandom covariate")
ax.set_title("The ordering is preserved, so what is\nordered is not identity-specific",
             fontsize=7.2, pad=4, color=C["muted"], linespacing=1.35)
ax.set_aspect("equal")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

# Panel b carries four pools in four colors and had no key at all, which made
# it unreadable on paper. The key goes in the panel, not the caption.
import matplotlib.lines as mlines  # noqa: E402
_h = [mlines.Line2D([], [], marker="o", ls="none", ms=3.6, color=COL[k], label=v)
      for k, v in (("matched", "matched to D1"), ("unrestricted", "unrestricted"),
                   ("no_tumor_shift", "no tumor shift"), ("rises_in_tumor", "rises in tumor"))]
axb.legend(handles=_h, frameon=False, fontsize=6.4, loc="lower left",
           handlelength=0.6, borderpad=0.05, labelspacing=0.18,
           handletextpad=0.35, borderaxespad=0.15)

for xx, yy, L in [(0.004, 0.985, "a"), (0.386, 0.985, "b"), (0.700, 0.985, "c")]:
    fig.text(xx, yy, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = (f"{IR_DOCS}/Figure12_negative_control"
        if _os.path.isdir(IR_DOCS) and _os.access(IR_DOCS, _os.W_OK)
        else f"{IR_FIGURES}/Figure12_negative_control")
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png").convert("RGB")
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size)
