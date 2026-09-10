"""Figure 7: how precise the retention fraction is, and whether the estimator
recovers a known value (PREREG 6v). Every value is read from the archived result
files; nothing is recomputed here.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from PIL import Image

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.2,
    "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})
MM = 1 / 25.4
C = {"red": "#199e70", "dra": "#d95926", "ink": "#0b0b0b", "muted": "#898781",
     "rule": "#c3c2b7", "pt": "#8d9bb5", "bad": "#d03b3b", "big": "#4a3aa7",
     "small": "#c99a12"}

R = "/home/claude/rsi/results"
A = json.load(open(f"{R}/RETENTION_INTERVALS_6v_GSE14520.json"))
A2 = json.load(open(f"{R}/RETENTION_INTERVALS_6v_GSE76427.json"))
S = json.load(open(f"{R}/ESTIMATOR_SIMULATION_6v.json"))

fig = plt.figure(figsize=(174 * MM, 132 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.075, right=0.985,
                       top=0.92, bottom=0.085, hspace=0.46, wspace=0.30)

# --- (a) caterpillar --------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
sig = [r for r in A["signatures"] if not r["ratio_unstable"]]
sig.sort(key=lambda r: r["retention_joint"])
v = np.array([r["retention_joint"] for r in sig])
lo = np.array([r["ci95_retention_joint"][0] for r in sig])
hi = np.array([r["ci95_retention_joint"][1] for r in sig])
x = np.arange(len(sig))
supported = hi < 0.5
ax.vlines(x[~supported], lo[~supported], hi[~supported], lw=0.55,
          color=C["rule"], zorder=2)
ax.vlines(x[supported], lo[supported], hi[supported], lw=0.8,
          color=C["bad"], alpha=0.75, zorder=3)
ax.plot(x, v, ".", ms=2.0, color=C["ink"], zorder=4)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=5)
ax.axhline(1.0, ls=(0, (1, 2)), lw=0.7, color=C["muted"], zorder=5)
for mod, col in (("DRAIN", C["dra"]), ("REDUCTION", C["red"])):
    m = A["own_modules"][mod]
    c = m["ci95_retention_joint"]
    ax.axhline(m["retention_joint"], lw=1.0, color=col, zorder=6)
    # The label used to sit on its own reference line, with the caterpillar of
    # point estimates running through the letters. Offset it from the line and
    # give it a white backing so nothing overprints it.
    dy = -0.11 if m["retention_joint"] > 0.5 else 0.11
    ax.text(len(sig) * 0.995, m["retention_joint"] + dy, f"{mod}", color=col,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9),
            zorder=7,
            fontsize=7.2, fontweight="bold", ha="right", va="center")
ax.text(1, 2.30, f"{int(supported.sum())} of {A['n_below_50_point']} "
        f"below 50% keep\nthe whole interval below 50%", fontsize=7.2,
        color=C["bad"], ha="left", va="top")
ax.set_xlim(-2, len(sig) + 1)
ax.set_ylim(-0.05, 2.45)
n_clip = int(sum(1 for r in sig if r.get("ci95_retention_joint", [0, 0])[1] > 2.45))
ax.text(len(sig) * 0.99, 0.02,
        f"axis trimmed at 2.45; {n_clip} intervals extend above it",
        fontsize=7.2, color=C["muted"], ha="right", va="bottom")
ax.set_xlabel("Signatures, ranked by retention")
ax.set_ylabel("Identity-retention fraction\nwith 95% bootstrap interval")
ax.set_title(f"GSE14520, {A['n_pairs']} pairs", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (b) interval width by cohort size -------------------------------------
ax = fig.add_subplot(gs[0, 1])
w1 = np.array([r["ci95_retention_joint"][1] - r["ci95_retention_joint"][0]
               for r in A["signatures"] if not r["ratio_unstable"]])
w2 = np.array([r["ci95_retention_joint"][1] - r["ci95_retention_joint"][0]
               for r in A2["signatures"] if not r["ratio_unstable"]])
bins = np.arange(0, 2.2, 0.1)
ax.hist(w2, bins=bins, color=C["small"], alpha=0.62, edgecolor="white", lw=0.4,
        label=f"GSE76427, {A2['n_pairs']} pairs", zorder=3)
ax.hist(w1, bins=bins, color=C["big"], alpha=0.62, edgecolor="white", lw=0.4,
        label=f"GSE14520, {A['n_pairs']} pairs", zorder=3)
for arr, col in ((w2, C["small"]), (w1, C["big"])):
    ax.axvline(np.median(arr), lw=1.1, color=col, ls=":", zorder=4)
ax.legend(frameon=False, fontsize=7.2, loc="upper right", handlelength=0.9,
          borderpad=0.1, labelspacing=0.22, handletextpad=0.4)
ax.text(0.985, 0.60, f"median width\n{np.median(w1):.3f} and {np.median(w2):.3f}\n"
        f"supported below 50%\n{A['n_below_50_interval']} and "
        f"{A2['n_below_50_interval']}",
        transform=ax.transAxes, fontsize=7.2, color=C["ink"], ha="right",
        va="top")
ax.set_xlabel("Width of the 95% interval")
ax.set_ylabel("Signatures")
ax.set_title("Precision depends on the number of pairs", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (c) simulation, recovery ----------------------------------------------
ax = fig.add_subplot(gs[1, 0])
cells = [c for c in S["cells"] if c["mu"] == 0.50 and c["comp_share"] == 0.0]
sigmas = sorted({c["sigma"] for c in cells})
cols = {0.25: "#199e70", 0.50: "#4a3aa7", 1.00: "#d03b3b"}
ax.plot([0, 1.1], [0, 1.1], lw=0.9, color=C["rule"], zorder=1)
for s in sigmas:
    cc = sorted([c for c in cells if c["sigma"] == s], key=lambda c: c["tau"])
    t = [c["tau"] for c in cc]; m = [c["median_estimate"] for c in cc]
    q1 = [c["iqr_estimate"][0] for c in cc]; q3 = [c["iqr_estimate"][1] for c in cc]
    ax.fill_between(t, q1, q3, color=cols[s], alpha=0.13, lw=0, zorder=2)
    ax.plot(t, m, "-o", lw=1.1, ms=3.2, color=cols[s], zorder=3,
            label=f"σ = {s:.2f}")
ax.annotate("upward bias only where\ntrue retention is near zero\nand noise is high",
            xy=(0.10, 0.154), xytext=(0.30, 0.10), fontsize=7.2, color=C["bad"],
            ha="left", va="center",
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C["bad"],
                            shrinkA=2, shrinkB=3))
ax.legend(frameon=False, fontsize=7.2, loc="upper left", handlelength=1.1,
          borderpad=0.15, labelspacing=0.25)
ax.set_xlim(0, 1.1); ax.set_ylim(0, 1.15)
ax.set_xlabel("True retention fraction")
ax.set_ylabel("Estimated retention\n(median and interquartile range)")
ax.set_title(f"{S['n_rep']} replicates per point, on the real covariates",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (d) coverage -----------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axhspan(0.92, 0.97, color=C["rule"], alpha=0.30, lw=0, zorder=1)
ax.axhline(0.95, lw=0.8, color=C["ink"], ls=(0, (3, 2)), zorder=2)
jit = {0.25: -0.022, 0.50: 0.0, 1.00: 0.022}
for c in S["cells"]:
    ax.plot(c["tau"] + jit[c["sigma"]], c["coverage_95"], "o", ms=3.6,
            color=cols[c["sigma"]], alpha=0.85, zorder=3,
            markeredgecolor="white", markeredgewidth=0.35)
ax.text(1.09, 0.9015, "shaded: the 0.92 to 0.97 band fixed in §6v",
        fontsize=7.2, color=C["muted"], ha="right", va="bottom")
ax.set_xlim(0.0, 1.10); ax.set_ylim(0.895, 1.0)
ax.set_xticks([0.10, 0.25, 0.50, 0.75, 1.00])
ax.set_xlabel("True retention fraction")
ax.set_ylabel("Coverage of the nominal 95% interval")
ax.set_title(f"achieved {S['coverage_min']:.3f} to {S['coverage_max']:.3f} "
             f"across {len(S['cells'])} cells", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for xx, yy, L in [(0.004, 0.985, "a"), (0.505, 0.985, "b"),
                  (0.004, 0.495, "c"), (0.505, 0.495, "d")]:
    fig.text(xx, yy, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = "/home/claude/Figure7_uncertainty_and_recovery"
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png").convert("RGB")
im.save(stem + ".png", dpi=(300, 300))
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size, "| supported", int(supported.sum()),
      "| widths", round(float(np.median(w1)),3), round(float(np.median(w2)),3))
