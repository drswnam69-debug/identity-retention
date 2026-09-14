#!/usr/bin/env python3
"""23_graphical_abstract.py -- mandatory graphical abstract for Hepatology International.

Journal specification (link.springer.com/journal/12072/submission-guidelines):
single-panel image, 300 dpi or greater, .tif/.jpg, Arial 12-16 pt, RGB,
reads top to bottom or left to right, own title, bullet labels, emphasis on the
new findings and a future vision.

Liberation Sans is metric-compatible with Arial and is substituted here; the
production file should be re-rendered with true Arial if the publisher requires
the exact face.
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F

OUT = "results/figures"
C = {c: json.load(open(f"results/{c}_composition.json")) for c in ("GSE76427", "GSE14520")}
assert abs(C["GSE14520"]["composite_joint"]["estimate"] - 0.814) < 5e-3
assert abs(C["GSE76427"]["composite_joint"]["estimate"] - 0.451) < 5e-3

INK, MUT, RULE = "#12181a", "#5d6b69", "#c9d4d1"
SUP, RED, DRA = F.C["supply"], F.C["protect"], F.C["patho"]

fig = plt.figure(figsize=(7.6, 6.6), dpi=300, facecolor="white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, ec, fc, r=1.4, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                lw=lw, edgecolor=ec, facecolor=fc, zorder=2))

def arrow(p0, p1, color, lw=1.8, ms=11):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                                 lw=lw, color=color, zorder=3,
                                 shrinkA=1.5, shrinkB=1.5))

# ---------------- title ----------------
ax.text(50, 96.4, "Reductive supply in liver cancer:", ha="center", va="center",
        fontsize=16, fontweight="bold", color=INK)
ax.text(50, 92.2, "how much survives adjustment for what the tumor is made of",
        ha="center", va="center", fontsize=12, color=MUT)
ax.plot([6, 94], [89.2, 89.2], color=RULE, lw=1.2)

# ---------------- 1. the question ----------------
ax.text(6, 85.6, "1", fontsize=12, fontweight="bold", color=DRA, va="center")
ax.text(10, 85.6, "One electron pool, two demands", fontsize=13,
        fontweight="bold", color=INK, va="center")

box(8, 71.0, 24, 11.5, SUP, "#e9e6f8")
ax.text(20, 79.4, "SUPPLY", ha="center", fontsize=12, fontweight="bold", color=SUP)
ax.text(20, 75.0, "mevalonate → CoQ\nbiosynthesis", ha="center", va="center",
        fontsize=12, color=INK, linespacing=1.35)

box(38, 71.0, 24, 11.5, RED, "#dcf1e9")
ax.text(50, 79.4, "REDUCTION", ha="center", fontsize=12, fontweight="bold", color=RED)
ax.text(50, 75.0, "CYB5R3 · FSP1\nkeeps CoQ reduced", ha="center", va="center",
        fontsize=12, color=INK, linespacing=1.35)

box(68, 71.0, 24, 11.5, DRA, "#fbe6da")
ax.text(80, 79.4, "DRAIN", ha="center", fontsize=12, fontweight="bold", color=DRA)
ax.text(80, 75.0, "mARC1 · mARC2\nN-reduction (+ POR, §6l)", ha="center", va="center",
        fontsize=12, color=INK, linespacing=1.35)

ax.text(50, 67.6, "index  =  ½ (SUPPLY + REDUCTION)  −  DRAIN", ha="center",
        va="center", fontsize=13, fontweight="bold", color=INK)
ax.text(50, 64.2, "locked with a published hash before any cohort was analyzed",
        ha="center", va="center", fontsize=12, color=MUT)

ax.plot([6, 94], [61.4, 61.4], color=RULE, lw=1.0)

# ---------------- 2. what was done ----------------
ax.text(6, 58.0, "2", fontsize=12, fontweight="bold", color=DRA, va="center")
ax.text(10, 58.0, "1,086 liver transcriptomes · 6 public cohorts · 265 matched pairs",
        fontsize=13, fontweight="bold", color=INK, va="center")

# ---------------- 3. the cascade ----------------
ys = 42.0
stages = [("unadjusted", "+1.13", "+1.26", MUT, "#f2f4f3"),
          ("− hepatocyte\nidentity", "+0.57", "+0.79", INK, "#eef1f0"),
          ("− identity AND\ncomposition", "+0.45", "+0.81", RED, "#dcf1e9")]
xs = [10, 39, 68]
for (lab, v76, v14, col, fc), x in zip(stages, xs):
    box(x, ys, 22, 12.0, col, fc, lw=1.6 if col == RED else 1.2)
    ax.text(x + 11, ys + 9.4, lab, ha="center", va="center", fontsize=12,
            color=MUT, linespacing=1.25)
    ax.text(x + 11, ys + 5.0, f"{v14}", ha="center", va="center",
            fontsize=16, fontweight="bold", color=col)
    ax.text(x + 11, ys + 1.5, f"{v76} · 52 pairs", ha="center",
            va="center", fontsize=12, color=MUT)
for x in xs[:-1]:
    arrow((x + 22.6, ys + 6), (x + 28.4, ys + 6), RULE, lw=2.0, ms=12)
ax.text(50, ys - 3.8, "the rise at zero change in the covariate — "
        "large figure 213 pairs, small figure 52 pairs", ha="center",
        va="center", fontsize=12, color=INK)

ax.plot([6, 94], [34.6, 34.6], color=RULE, lw=1.0)

# ---------------- 4. findings ----------------
ax.text(6, 31.2, "3", fontsize=12, fontweight="bold", color=DRA, va="center")
ax.text(10, 31.2, "What is new", fontsize=13, fontweight="bold", color=INK, va="center")

bullets = [
    (RED, "Cell composition explains none of the rise", "adjusting for leukocyte, endothelial and stromal content leaves 101–116%"),
    (INK, "Loss of hepatocyte identity explains about half", "the drain arm is mostly identity loss; the reduction arm mostly is not — 89% vs 31%"),
    (DRA, "Half the pre-registered hypotheses failed", "survival prediction among them — reported as specified, not reframed"),
]
y = 26.6
for col, head, sub in bullets:
    ax.plot([9.4, 11.6], [y + 0.55, y + 0.55], color=col, lw=2.6,
            solid_capstyle="round")
    ax.text(13.0, y + 0.9, head, fontsize=12.5, fontweight="bold", color=col, va="center")
    ax.text(13.0, y - 2.2, sub, fontsize=12, color=MUT, va="center")
    y -= 7.2

ax.plot([6, 94], [8.4, 8.4], color=RULE, lw=1.0)

# ---------------- 5. future vision ----------------
ax.text(6.4, 4.6, "→", fontsize=14, fontweight="bold", color=RED, va="center")
ax.text(10.4, 4.6, "Next", fontsize=12.5, fontweight="bold", color=RED, va="center")
ax.text(19.5, 4.6, "treat the reduction arm, not the drain, as the tractable node — "
        "and test it in a\npaired cohort carrying AIFM2, at protein level, "
        "before any clinical claim.",
        fontsize=12, color=INK, va="center", linespacing=1.55)

for ext in ("tif", "jpg", "png"):
    path = os.path.join(OUT, f"GraphicalAbstract.{ext}")
    kw = dict(dpi=300, facecolor="white")
    if ext == "tif":
        kw["pil_kwargs"] = {"compression": "tiff_lzw"}
    if ext == "jpg":
        kw["pil_kwargs"] = {"quality": 95}
    fig.savefig(path, **kw)
    print("wrote", path)
