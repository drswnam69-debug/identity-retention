"""Figure 9: the measure in three tissues (PREREG 6x, 6y). Every value is read
from archived result files; nothing is recomputed here."""
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
C = {"liv": "#4a3aa7", "lun": "#199e70", "kid": "#d95926", "ink": "#0b0b0b",
     "muted": "#898781", "rule": "#c3c2b7", "id": "#4a3aa7", "co": "#c99a12",
     "bad": "#d03b3b", "pt": "#8d9bb5"}
R = "/home/claude/rsi/results"
B = {k: json.load(open(f"{R}/SIGNATURE_BENCHMARK_{v}.json"))
     for k, v in [("liver", "GSE14520"), ("lung", "TCGA_LUAD"), ("kidney", "TCGA_KIRC")]}
PK = json.load(open(f"{R}/TCGA_KIRC/PREMISE_6x.json"))
PL = json.load(open(f"{R}/TCGA_LUAD/PREMISE_6y.json"))
PH = json.load(open(f"{R}/TCGA_LIHC/PREMISE_6w.json"))
U = {k: json.load(open(f"{R}/COVARIATE_COMPARISON_6u_{v}.json"))
     for k, v in [("liver GSE14520", "GSE14520"), ("liver TCGA", "TCGA_LIHC"),
                  ("lung", "TCGA_LUAD"), ("kidney", "TCGA_KIRC")]}
V = {k: json.load(open(f"{R}/RETENTION_INTERVALS_6v_{v}.json"))
     for k, v in [("liver GSE14520", "GSE14520"), ("liver GSE76427", "GSE76427"),
                  ("liver TCGA", "TCGA_LIHC"), ("lung", "TCGA_LUAD"),
                  ("kidney", "TCGA_KIRC")]}
COL = {"liver": C["liv"], "lung": C["lun"], "kidney": C["kid"]}

fig = plt.figure(figsize=(174 * MM, 134 * MM))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.082, right=0.985,
                       top=0.915, bottom=0.088, hspace=0.50, wspace=0.32)

# --- (a) three distributions ------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
vals = {t_: np.array([r["retention_joint"] for r in B[t_]["signatures"]
                      if r["status"] == "ok"]) for t_ in ("liver", "lung", "kidney")}
rng = np.random.default_rng(7)
for i, t_ in enumerate(("liver", "lung", "kidney")):
    v = vals[t_]
    ax.scatter(v, np.full(len(v), i) + rng.normal(0, 0.085, len(v)), s=9,
               color=COL[t_], alpha=0.45, edgecolor="none", zorder=3)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([q1, q3], [i - 0.26] * 2, lw=3.0, color=COL[t_], solid_capstyle="butt",
            zorder=4)
    ax.plot([med, med], [i - 0.34, i - 0.18], lw=1.6, color=C["ink"], zorder=5)
    COHORT = {"liver": "GSE14520", "lung": "TCGA-LUAD", "kidney": "TCGA-KIRC"}
    ax.text(1.98, i + 0.20,
            f"{t_} ({COHORT[t_]}), {len(v)} sets   median {med:.3f}",
            fontsize=7.2, color=COL[t_], ha="right", va="bottom")
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=2)
ax.axvline(1.0, ls=(0, (1, 2)), lw=0.7, color=C["muted"], zorder=2)
ax.set_yticks([])
ax.set_xlim(-0.03, 2.0)
ax.set_ylim(-0.80, 2.62)
n_off = sum(int((vals[t_] > 2.0).sum()) for t_ in vals)
ax.text(0.02, -0.76, f"axis trimmed at 2.0; {n_off} sets lie above it",
        fontsize=7.2, color=C["muted"], ha="left", va="bottom")
ax.set_xlabel("Identity-retention fraction")
ax.set_title("Dispersion generalizes; the level is tissue-specific", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)

# --- (b) the premise in each tissue -----------------------------------------
ax = fig.add_subplot(gs[0, 1])
P = [("liver TCGA-LIHC\nD1, 50 pairs", PH, C["liv"]),
     ("lung TCGA-LUAD\nL1, 58 pairs", PL, C["lun"]),
     ("kidney TCGA-KIRC\nK1, 72 pairs", PK, C["kid"])]
x = np.arange(3)
for i, (lab, p, col) in enumerate(P):
    key = "per_gene_D1" if "per_gene_D1" in p else "per_gene_identity"
    eff = [g["paired_delta"] for g in p[key]
           if g.get("class", "effector") == "effector"]
    reg = [g["paired_delta"] for g in p[key] if g.get("class") == "regulator"]
    if not reg:   # liver file has no class field; split by the known regulator list
        TF = {"HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4"}
        eff = [g["paired_delta"] for g in p[key] if g["gene"] not in TF]
        reg = [g["paired_delta"] for g in p[key] if g["gene"] in TF]
    ax.scatter(np.full(len(eff), i - 0.13) + np.random.default_rng(1).normal(0, .022, len(eff)),
               eff, s=11, color=col, alpha=0.55, edgecolor="none", zorder=3)
    ax.scatter(np.full(len(reg), i + 0.16) + np.random.default_rng(2).normal(0, .018, len(reg)),
               reg, s=22, marker="D", color=col, edgecolor="white", lw=0.4, zorder=4)
    ax.plot([i - 0.24, i - 0.02], [np.mean(eff)] * 2, lw=1.5, color=C["ink"], zorder=5)
    ax.plot([i + 0.05, i + 0.27], [np.mean(reg)] * 2, lw=1.5, color=C["ink"], zorder=5)
ax.axhline(0, lw=0.7, color=C["ink"])
ax.set_xticks(x)
ax.set_xticklabels([p[0] for p in P], fontsize=7.2)
ax.set_ylabel("Identity gene, tumor minus adjacent\n(z units)")
ax.text(0.02, 0.02, "circles: effectors\ndiamonds: regulators\nbars: group means",
        transform=ax.transAxes, fontsize=7.2, color=C["muted"], ha="left",
        va="bottom", linespacing=1.35,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
ax.set_title("The premise holds in all three, differently", fontsize=7.2, pad=4,
             color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (c) identity against composition, four cohorts, plus the 6ac sensitivity
ax = fig.add_subplot(gs[1, 0])
order = ["liver GSE14520", "liver TCGA", "lung", "kidney"]
mi = [U[k]["retention_identity_only"]["median"] for k in order]
mc = [U[k]["retention_composition_only"]["median"] for k in order]
NOT_LUNG = {"MCCLUNG_COCAINE_REWARD_5D", "MCCLUNG_COCAIN_REWARD_4WK",
            "MCCLUNG_CREB1_TARGETS_DN", "MCCLUNG_CREB1_TARGETS_UP",
            "MCCLUNG_DELTA_FOSB_TARGETS_8WK", "REN_ALVEOLAR_RHABDOMYOSARCOMA_UP"}
keep = [r for r in U["lung"]["signatures"]
        if r["status"] == "ok" and r["name"] not in NOT_LUNG]
mi.append(float(np.median([r["retention_D1"] for r in keep])))
mc.append(float(np.median([r["retention_C1"] for r in keep])))
xi = np.arange(5); w = 0.36
hatch = [None, None, None, None, "///"]
for j in range(5):
    ax.bar(xi[j] - w / 2, mi[j], width=w, color=C["id"], edgecolor="white", lw=0.6,
           hatch=hatch[j], zorder=3, label="identity covariate" if j == 0 else None)
    ax.bar(xi[j] + w / 2, mc[j], width=w, color=C["co"], edgecolor="white", lw=0.6,
           hatch=hatch[j], zorder=3, label="composition covariate" if j == 0 else None)
ax.axhline(1.0, ls=(0, (1, 2)), lw=0.8, color=C["muted"], zorder=4)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=4)
for j in range(5):
    ax.text(xi[j], max(mi[j], mc[j]) + 0.04, f"{mi[j]:.2f} | {mc[j]:.2f}",
            ha="center", va="bottom", fontsize=7.2)
ax.annotate("composition removes nothing", xy=(0 + w / 2, mc[0]),
            xytext=(1.55, 1.50), fontsize=7.2, color=C["bad"], ha="left", va="center",
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C["bad"], shrinkA=2, shrinkB=3))
ax.annotate("the lung ordering reverses only\non the locked panel (\u00a76ac)",
            xy=(2 + w / 2, mc[2]), xytext=(4.55, 1.12), fontsize=7.2, color=C["bad"],
            ha="right", va="center", linespacing=1.3,
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color=C["bad"], shrinkA=2, shrinkB=3))
ax.set_xticks(xi)
ax.set_xticklabels(["liver\nGSE14520", "liver\nLIHC", "lung\nLUAD",
                    "kidney\nKIRC", "lung, 30 sets\n\u00a76ac"], fontsize=7.2)
ax.set_ylim(0, 1.62)
ax.set_ylabel("Median retention, single-covariate model")
ax.legend(frameon=False, fontsize=7.2, loc="upper left", handlelength=0.9,
          borderpad=0.1, labelspacing=0.25, bbox_to_anchor=(-0.02, 1.02))
ax.set_title("What each covariate removes is tissue-dependent", fontsize=7.2,
             pad=4, color=C["muted"])
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

# --- (d) precision against the number of pairs ------------------------------
ax = fig.add_subplot(gs[1, 1])
pts = [(V[k]["n_pairs"], V[k]["median_ci_width"], V[k]["n_below_50_interval"],
        V[k]["n_below_50_point"], k) for k in V]
cols = {"liver GSE14520": C["liv"], "liver GSE76427": C["liv"],
        "liver TCGA": C["liv"], "lung": C["lun"], "kidney": C["kid"]}
OFF = {"lung": (12, 0.07, "left", "bottom"), "kidney": (12, 0.05, "left", "bottom"),
       "liver GSE76427": (13, 0.05, "left", "bottom"),
       "liver TCGA": (13, -0.07, "left", "top"),
       "liver GSE14520": (0, -0.07, "center", "top")}
for n, wdt, sup, pt, k in pts:
    dx, dy, ha, va = OFF[k]
    ax.scatter([n], [wdt], s=46, color=cols[k], edgecolor="white", lw=0.6, zorder=4)
    NAME = {"liver GSE14520": "liver GSE14520",
            "liver GSE76427": "liver GSE76427",
            "liver TCGA": "liver TCGA-LIHC",
            "lung": "lung TCGA-LUAD", "kidney": "kidney TCGA-KIRC"}
    ax.annotate(f"{NAME[k]}\n{sup} of {pt}",
                xy=(n, wdt), xytext=(n + dx, wdt + dy), fontsize=7.2,
                color=C["ink"], ha=ha, va=va)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=2)
ax.text(243, 0.52, "interval width 0.50", fontsize=7.2, color=C["ink"],
        ha="right", va="bottom")
ax.set_xlim(30, 245)
ax.set_ylim(0, 1.80)
ax.set_xlabel("Paired tumor and adjacent cases (pairs)")
ax.set_ylabel("Median width of the 95% interval")
ax.set_title("Only the largest cohort supports per-signature values\n"
             "labels: below 50%, how many keep the whole interval there",
             fontsize=7.2, pad=4, color=C["muted"], linespacing=1.4)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for xx, yy, L in [(0.004, 0.985, "a"), (0.505, 0.985, "b"),
                  (0.004, 0.495, "c"), (0.505, 0.495, "d")]:
    fig.text(xx, yy, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = "/home/claude/Figure9_three_tissues"
fig.savefig(stem + ".png", dpi=300, facecolor="white")
fig.savefig(stem + ".svg", facecolor="white")
plt.close(fig)
im = Image.open(stem + ".png").convert("RGB")
im.save(stem + ".png", dpi=(300, 300))
im.save(stem + ".tiff", compression="tiff_lzw", dpi=(300, 300))
print("written", im.size)
