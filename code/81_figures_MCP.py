#!/usr/bin/env python3
"""81_figures_MCP.py

Three figures for manuscript_MCP.md, drawn from results/MANUSCRIPT_NUMBERS.json
and results/PANEL_COHERENCE.json only. Each is written as PNG at 300 dpi, TIFF
with LZW compression and SVG.

Encoding held constant across all three figures:
  hue      protein class      metabolic enzymes indigo, plasma proteins orange
  fill     measurement level  filled transcriptome, open proteome
  shape    cohort             circle, square, triangle where two or three are shown
Colors are the validated four-hue categorical palette; the two used here are
adjacent in it and clear the colorblind separation check, and every series also
carries a legend entry and a direct label, so identity is never color alone.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/figures/mcp"
os.makedirs(OUT, exist_ok=True)

METAB = "#4a3aa7"     # metabolic enzymes
PLASMA = "#d95926"    # secreted plasma proteins
NEUTRAL = "#8a8a84"   # the comparison-set distribution
INK = "#1c1c1a"
MUTED = "#6b6b66"
GRID = "#e3e3df"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.titlesize": 8.5, "legend.fontsize": 7,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


def save(fig, stem):
    png = f"{OUT}/{stem}.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(f"{OUT}/{stem}.svg", bbox_inches="tight")
    Image.open(png).convert("RGB").save(f"{OUT}/{stem}.tiff", compression="tiff_lzw")
    plt.close(fig)
    print(f"  {stem}: png, tiff, svg")


def tidy(ax, xlabel=None, ylabel=None, xgrid=False, ygrid=False):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if xgrid:
        ax.xaxis.grid(True, color=GRID, lw=0.5)
    if ygrid:
        ax.yaxis.grid(True, color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


# --------------------------------------------------------------- figure 1
def figure1(num, pc):
    order = ["GSE14520 (array transcriptome)", "GSE76427 (array transcriptome)",
             "TCGA-LIHC (RNA sequencing)", "Gao 2019 proteome", "Jiang 2019 proteome"]
    short = ["GSE14520", "GSE76427", "TCGA-LIHC", "Gao 2019", "Jiang 2019"]
    pdrop = num["comparison_sets"]["panel_drop"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1),
                             gridspec_kw={"width_ratios": [1.35, 1]})

    ax = axes[0]
    rng = np.random.default_rng(0)
    for i, name in enumerate(order):
        vals = [v["cancellation"] for v in
                pc["datasets"][name]["gene_set_coherence"].values()
                if v and v.get("cancellation") is not None]
        x = i + (rng.random(len(vals)) - 0.5) * 0.5
        ax.scatter(x, vals, s=3.2, color=NEUTRAL, alpha=0.5, linewidths=0, zorder=2)
        ax.plot([i - 0.33, i + 0.33], [np.median(vals)] * 2,
                color=MUTED, lw=1.6, zorder=3, solid_capstyle="round")
        p = pdrop["per_dataset"][name]["panel_surviving_share"]
        proteome = "proteome" in name
        ax.scatter([i], [p], s=64, marker="D", zorder=5,
                   facecolor=SURFACE if proteome else METAB,
                   edgecolor=METAB, linewidths=1.5)
        ax.annotate(f"{p:.2f}", (i, p), textcoords="offset points",
                    xytext=(0, 9 if proteome else 9), ha="center",
                    fontsize=7, color=INK)
    ax.axvline(2.5, color=GRID, lw=1.2, zorder=1)
    ax.text(1.0, 1.06, "transcriptome", ha="center", fontsize=7, color=MUTED)
    ax.text(3.5, 1.06, "proteome", ha="center", fontsize=7, color=MUTED)
    ax.set_xticks(range(5)); ax.set_xticklabels(short, rotation=20, ha="right")
    ax.set_ylim(-0.02, 1.13); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    tidy(ax, ylabel="Surviving share", ygrid=True)
    ax.set_title("a   The panel against 114 published liver signatures", loc="left")
    ax.legend(handles=[
        Line2D([], [], marker="D", ls="", markerfacecolor=METAB,
               markeredgecolor=METAB, markersize=7, label="22-gene panel"),
        Line2D([], [], marker="D", ls="", markerfacecolor=SURFACE,
               markeredgecolor=METAB, markersize=7, label="panel, proteome"),
        Line2D([], [], marker="o", ls="", color=NEUTRAL, markersize=4,
               label="114 signatures"),
        Line2D([], [], color=MUTED, lw=1.6, label="their median")],
        loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
        frameon=False, handletextpad=0.5, columnspacing=1.4, borderpad=0.2)

    ax = axes[1]
    drops = [r["transcriptome_cancellation"] - r["proteome_cancellation"]
             for r in pc["coherence_by_level"]["per_set"]]
    ax.hist(drops, bins=26, color=NEUTRAL, alpha=0.75, edgecolor=SURFACE, lw=0.4)
    ax.axvline(pdrop["panel_drop"], color=METAB, lw=1.8)
    ax.annotate(f"panel  {pdrop['panel_drop']:.3f}",
                (pdrop["panel_drop"], ax.get_ylim()[1] * 0.93),
                xytext=(-6, 0), textcoords="offset points",
                ha="right", fontsize=7.5, color=METAB)
    ax.annotate(f"median {pdrop['set_drop_median']:.3f}\nlargest {pdrop['set_drop_max']:.3f}",
                (pdrop["set_drop_median"], ax.get_ylim()[1] * 0.55),
                xytext=(10, 0), textcoords="offset points",
                ha="left", fontsize=7, color=MUTED)
    ax.set_xlim(min(drops) - 0.05, pdrop["panel_drop"] + 0.09)
    tidy(ax, xlabel="Loss of surviving share between levels",
         ylabel="Number of signatures", ygrid=True)
    ax.set_title("b   No signature loses as much as the panel", loc="left")
    fig.tight_layout(w_pad=2.0)
    save(fig, "Figure1_panel_comes_apart")


# --------------------------------------------------------------- figure 2
def figure2(num):
    rows = num["datasets"]
    groups = [("Liver transcriptomes", [r for r in rows if r["level"] == "transcriptome"]),
              ("Hepatocellular carcinoma proteomes",
               [r for r in rows if r["level"] == "proteome"
                and r["tumor_type"] == "HCC"]),
              ("Cholangiocarcinoma proteomes (exploratory)",
               [r for r in rows if r["tumor_type"] != "HCC"])]
    labels, ys, band = [], [], []
    y = 0
    for title, rs in groups:
        band.append((title, y, y + len(rs) - 1))
        for r in rs:
            labels.append(r["dataset"].replace(" (array transcriptome)", "")
                          .replace(" (RNA sequencing)", "").replace(" proteome", "")
                          .replace(" (held out)", "").replace(" (exploratory)", ""))
            ys.append((y, r)); y += 1
        y += 0.8

    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    ax.axvline(0, color=MUTED, lw=0.8, zorder=2)
    for yy, r in ys:
        proteome = r["level"] == "proteome"
        for val, share, color, dy in (
                (r["metabolic_mean_z"], r["metabolic_share"], METAB, -0.16),
                (r["secreted_mean_z"], r["secreted_share"], PLASMA, 0.16)):
            ax.plot([0, val], [-(yy + dy)] * 2, color=color, lw=1.6,
                    alpha=0.5, solid_capstyle="round", zorder=3)
            ax.scatter([val], [-(yy + dy)], s=46,
                       facecolor=SURFACE if proteome else color,
                       edgecolor=color, linewidths=1.4, zorder=4)
            ax.annotate(f"{val:+.2f}  ({share:.2f})",
                        (val, -(yy + dy)),
                        xytext=(7 if val > 0 else -7, 0),
                        textcoords="offset points", fontsize=6.6,
                        ha="left" if val > 0 else "right", va="center",
                        color=color)
    ax.set_yticks([-(yy) for yy, _ in ys])
    ax.set_yticklabels([f"{l}   n={r['n_pairs']}" for l, (_, r) in zip(labels, ys)])
    for title, a, b in band:
        ax.annotate(title, (-2.32, -(a) + 0.60), fontsize=7, color=MUTED,
                    ha="left", annotation_clip=False)
    ax.set_xlim(-2.35, 1.35)
    tidy(ax, xlabel="Mean paired shift, tumor minus adjacent\n"
                    "(SD of the pooled spread of both arms; "
                    "surviving share in parentheses)", xgrid=True)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", markerfacecolor=METAB,
               markeredgecolor=METAB, markersize=6.5,
               label="metabolic enzymes (transcriptome)"),
        Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE,
               markeredgecolor=METAB, markersize=6.5,
               label="metabolic enzymes (proteome)"),
        Line2D([], [], marker="o", ls="", markerfacecolor=PLASMA,
               markeredgecolor=PLASMA, markersize=6.5,
               label="plasma proteins (transcriptome)"),
        Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE,
               markeredgecolor=PLASMA, markersize=6.5,
               label="plasma proteins (proteome)")],
        loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
        frameon=False, handletextpad=0.5, columnspacing=1.6, borderpad=0.2)
    ax.set_ylim(-(y - 1.8) - 0.75, 1.30)
    fig.tight_layout()
    save(fig, "Figure2_eight_datasets")


# --------------------------------------------------------------- figure 3
def figure3(num):
    coh = num["cholangiocarcinoma"]["cohorts"]
    names = ["MSKCC-ICC", "UKF-ICC"]
    marks = {"MSKCC-ICC": "o", "UKF-ICC": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4),
                             gridspec_kw={"width_ratios": [1.28, 1]})

    ax = axes[0]
    order = ["SERPINA1", "ALB", "AHSG", "APOH", "TF", "TTR",
             "FGB", "FGG", "F2", "FGA"]
    ax.axvline(0, color=MUTED, lw=0.8, zorder=2)
    dodge = {"MSKCC-ICC": 0.15, "UKF-ICC": -0.15}
    for i, g in enumerate(order):
        for c in names:
            v = coh[c]["G3_plasma"]["members"].get(g)
            if v is None:
                continue
            ax.scatter([v], [-i + dodge[c]], s=42, marker=marks[c],
                       facecolor=SURFACE, edgecolor=PLASMA, linewidths=1.4,
                       zorder=4)
    ax.set_yticks([-i for i in range(len(order))]); ax.set_yticklabels(order, style="italic")
    ax.set_xlim(-0.85, 1.45)
    tidy(ax, xlabel="Mean paired shift", xgrid=True)
    ax.set_title("a   Plasma proteins in cholangiocarcinoma", loc="left")
    ax.annotate("rise", (1.35, -0.0), fontsize=7, color=PLASMA, ha="right")
    ax.annotate("fall", (-0.78, -9.0), fontsize=7, color=PLASMA, ha="left")

    ax = axes[1]
    morder = ["OTC", "PCK1", "CPS1", "ARG1"]
    ax.axvline(0, color=MUTED, lw=0.8, zorder=2)
    dodge2 = {"MSKCC-ICC": 0.15, "UKF-ICC": -0.15}
    for i, g in enumerate(morder):
        for c in names:
            v = coh[c]["G1_module"]["members"].get(g)
            if v is None:
                continue
            ax.scatter([v], [-i + dodge2[c]], s=42, marker=marks[c],
                       facecolor=SURFACE, edgecolor=METAB, linewidths=1.4,
                       zorder=4)
    for c, dy in zip(names, (-4.7, -5.4)):
        ax.annotate(f"{c}: mean {coh[c]['G1_module']['mean_z']:+.3f},"
                    f" surviving share {coh[c]['G1_module']['coherence']['cancellation']:.3f}",
                    (-2.35, dy), fontsize=6.8, color=METAB, ha="left",
                    annotation_clip=False)
    ax.set_yticks([-i for i in range(len(morder))])
    ax.set_yticklabels(morder, style="italic")
    ax.set_xlim(-2.4, 0.35); ax.set_ylim(-5.8, 0.6)
    tidy(ax, xlabel="Mean paired shift", xgrid=True)
    ax.set_title("b   The metabolic module in the same cohorts", loc="left")
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE,
               markeredgecolor=MUTED, markersize=6.5, label="MSKCC-ICC, 64 pairs"),
        Line2D([], [], marker="s", ls="", markerfacecolor=SURFACE,
               markeredgecolor=MUTED, markersize=6.5, label="UKF-ICC, 61 pairs"),
        Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE,
               markeredgecolor=PLASMA, markersize=6.5, label="plasma protein"),
        Line2D([], [], marker="o", ls="", markerfacecolor=SURFACE,
               markeredgecolor=METAB, markersize=6.5, label="metabolic enzyme")],
        loc="lower center", bbox_to_anchor=(0.5, -0.09), ncol=4,
        frameon=False, handletextpad=0.5, columnspacing=1.6)
    fig.tight_layout(w_pad=2.2, rect=(0, 0.04, 1, 1))
    save(fig, "Figure3_cholangiocarcinoma")


def main() -> int:
    num = json.load(open(f"{ROOT}/results/MANUSCRIPT_NUMBERS.json"))
    pc = json.load(open(f"{ROOT}/results/PANEL_COHERENCE.json"))
    print("drawing:")
    figure1(num, pc)
    figure2(num)
    figure3(num)
    print(f"\nwrote to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
