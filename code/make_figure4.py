"""Figure 4 for the v3 manuscript: identity retention across published liver
signatures (PREREG 6t). Every value is read from
results/SIGNATURE_BENCHMARK_GSE14520.json; nothing is recomputed here.
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
C = {"red": "#199e70", "dra": "#d95926", "ink": "#0b0b0b",
     "muted": "#898781", "rule": "#c3c2b7", "bar": "#b8c4d9",
     "bar_lo": "#8d9bb5"}

D = json.load(open("/home/claude/rsi/results/SIGNATURE_BENCHMARK_GSE14520.json"))
ok = [r for r in D["signatures"] if r["status"] == "ok"]
vals = np.array([r["retention_joint"] for r in ok], float)
names = [r["name"] for r in ok]
S = D["summary_joint"]
RED = D["own_modules"]["REDUCTION"]["retention_joint"]
DRA = abs(D["own_modules"]["DRAIN"]["retention_joint"])

fig = plt.figure(figsize=(174 * MM, 140 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.078, right=0.985,
                       top=0.93, bottom=0.085, hspace=0.48, wspace=0.42,
                       width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0])

# (a) histogram
ax = fig.add_subplot(gs[0, 0])
bins = np.arange(0, 2.5, 0.1)
# One bar color. The two-tone fill carried no key and the dashed 0.5 line
# already marks the threshold.
cols = [C["bar"] for b in bins[:-1]]
n, _, patches = ax.hist(vals, bins=bins, color=C["bar"], edgecolor="white", lw=0.4)
for p, c in zip(patches, cols):
    p.set_facecolor(c)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"])
top = max(n) * 1.02
for x, col, lab in [(DRA, C["dra"], "DRAIN"), (RED, C["red"], "REDUCTION")]:
    ax.plot([x, x], [0, top], lw=1.6, color=col, zorder=5)
    ax.text(x, top * 1.03, lab, color=col, fontsize=7.2, fontweight="bold",
            ha="center", va="bottom")
ax.axvline(S["median"], lw=1.0, color=C["ink"], ls=":")
ax.axvline(1.0, ls=(0, (1, 2)), lw=0.7, color=C["muted"], zorder=2)
from matplotlib.ticker import MaxNLocator
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.text(S["median"], top * 0.55, f"  median {S['median']:.3f}", fontsize=7.2,
        color=C["ink"], ha="left", va="center")
ax.set_xlim(0, 2.45)
ax.set_ylim(0, top * 1.18)
ax.set_xlabel("Identity-retention fraction, joint model")
ax.set_ylabel("Published liver signatures")
ax.text(2.40, top * 0.95, f"{S['n_below_50pct']} of {len(vals)}\nbelow 50%",
        fontsize=7.2, color=C["ink"], ha="right", va="top")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (b) ranked, ten lowest named
ax = fig.add_subplot(gs[0, 1])
order = np.argsort(vals)
sv = vals[order]
x = np.arange(len(sv))
ax.plot(x, sv, "-", lw=1.1, color=C["muted"])
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"])
ax.scatter(x[:10], sv[:10], s=9, color=C["bar_lo"], zorder=4)
for y, col, lab in [(DRA, C["dra"], "DRAIN"), (RED, C["red"], "REDUCTION")]:
    ax.axhline(y, lw=1.2, color=col)
    ax.text(len(sv) * 0.99, y, f" {lab} ", color=col, fontsize=7.2,
            fontweight="bold", ha="right",
            va="bottom" if lab == "REDUCTION" else "top")
ax.set_xlim(-2, len(sv) + 1)
ax.set_ylim(-0.05, 2.45)
ax.set_xlabel("Signatures, ranked by retention")
ax.set_ylabel("Identity-retention fraction")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (c) sensitivity
ax = fig.add_subplot(gs[1, 0])
pairs = [(r["retention_joint"],
          r["sensitivity_no_covariate_genes"]["retention_joint"])
         for r in ok
         if "sensitivity_no_covariate_genes" in r
         and r["sensitivity_no_covariate_genes"].get("status") == "ok"]
a = np.array([p[0] for p in pairs])
b = np.array([p[1] for p in pairs])
lim = max(a.max(), b.max()) * 1.05
ax.plot([0, lim], [0, lim], lw=0.8, color=C["rule"], zorder=1)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.5)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.5)
ax.scatter(a, b, s=13, color=C["bar_lo"], edgecolor="white", lw=0.4, zorder=3)
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
ax.set_xlabel("Retention, signature as published")
ax.set_ylabel("Retention, after removing genes\nshared with D1 or C1")
ax.set_title(f"{len(pairs)} of 60 sharing signatures, evaluable after trimming\n"
             f"median {np.median(a):.3f} to {np.median(b):.3f}",
             fontsize=7.2, pad=4, color=C["muted"])
ax.set_aspect("equal", adjustable="box")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (d) the ten lowest, named
axd = fig.add_subplot(gs[1, 1])
low_i = order[:10][::-1]
low_v = vals[low_i]
low_n = [names[i].replace("_", " ").title() for i in low_i]
yy = np.arange(10)
axd.barh(yy, low_v, height=0.66, color=C["bar_lo"], edgecolor="none", zorder=3)
axd.axvline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"])
axd.axvline(DRA, lw=1.3, color=C["dra"], zorder=4)
axd.text(DRA, 9.9, " DRAIN", color=C["dra"], fontsize=7.2,
         fontweight="bold", ha="left", va="top")
for y, v in zip(yy, low_v):
    if v > 0.12:
        axd.text(v - 0.012, y, f"{v:.3f}", va="center", ha="right",
                 fontsize=7.2, color="white", zorder=5)
    else:
        axd.text(v + 0.012, y, f"{v:.3f}", va="center", ha="left",
                 fontsize=7.2, color=C["ink"])
axd.set_yticks(yy)
axd.set_yticklabels(low_n, fontsize=7.2)
axd.set_xlim(0, 0.62)
axd.set_ylim(-0.7, 10.2)
axd.set_xlabel("Identity-retention fraction")
axd.set_title("The ten lowest-retaining signatures", fontsize=7.2, pad=4,
              color=C["muted"])
for sp in ("top", "right"):
    axd.spines[sp].set_visible(False)

for xx, yy, L in [(0.004, 0.985, "a"), (0.505, 0.985, "b"),
                  (0.004, 0.495, "c"), (0.505, 0.495, "d")]:
    fig.text(xx, yy, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = "/home/claude/Figure4_signature_benchmark"
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
print("written", im.size, "| n sets", len(vals), "| sensitivity pairs", len(pairs))
