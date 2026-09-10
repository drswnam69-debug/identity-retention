"""Figure 6: the identity-retention fraction is not a re-description of
composition (tumor purity) adjustment (PREREG 6u). Every value is read from
results/COVARIATE_COMPARISON_6u_GSE14520.json; nothing is recomputed here, with
the single exception of the two paired covariate vectors in panel (a), which are
rebuilt by the same code path and checked against the archived correlation.
"""
import json
import importlib.util, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
from PIL import Image

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7.2,
    "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})
MM = 1 / 25.4
C = {"red": "#199e70", "dra": "#d95926", "ink": "#0b0b0b", "muted": "#898781",
     "rule": "#c3c2b7", "id": "#4a3aa7", "co": "#c99a12", "bad": "#d03b3b",
     "pt": "#8d9bb5"}

R = "/home/claude/rsi"
D = json.load(open(f"{R}/results/COVARIATE_COMPARISON_6u_GSE14520.json"))
ok = D["signatures"]
rd = np.array([r["retention_D1"] for r in ok], float)
rc = np.array([r["retention_C1"] for r in ok], float)
A = D["covariate_association"]
AG = D["agreement"]

# --- rebuild the two paired covariate vectors by the study's own code path ---
sys.path.insert(0, f"{R}/code")
os.chdir(R)
def _load(n, f):
    s = importlib.util.spec_from_file_location(n, f"{R}/code/{f}")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
_da = _load("d", "12_differentiation_adjust.py")
_ca = _load("c", "18_composition_adjust.py")
_bm = _load("b", "30_signature_benchmark.py")
rsi = pd.read_csv(f"{R}/results/GSE14520/rsi.tsv", sep="\t", index_col=0)
ph = pd.read_csv(f"{R}/results/GSE14520/phenotype.tsv", sep="\t",
                 index_col=0).loc[rsi.index]
expr = pd.read_csv(f"{R}/data/GSE14520_symbols.tsv.gz", sep="\t",
                   index_col=0)[rsi.index]
is_t = (~ph["tissue"].astype(str).str.lower().str.contains("adjacent")).to_numpy()
pid = ph["patient_id"].astype(str)
fr = pd.DataFrame({"pid": pid.values, "t": is_t}, index=rsi.index)
def paired(s):
    f = fr.assign(v=s.to_numpy())
    tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
    nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
    k = sorted(set(tt.index) & set(nn.index))
    return (tt.loc[k] - nn.loc[k]).to_numpy(float)
dD1 = paired(_bm.zmean(expr, _da.D1)[0])
dC1 = paired(_bm.zmean(expr, _ca.C1)[0])
r_check = float(np.corrcoef(dD1, dC1)[0, 1])
assert abs(r_check - A["pearson_r"]) < 1e-3, (r_check, A["pearson_r"])
os.chdir("/home/claude")

fig = plt.figure(figsize=(174 * MM, 66 * MM))
gs = gridspec.GridSpec(1, 3, figure=fig, left=0.062, right=0.988,
                       top=0.86, bottom=0.185, wspace=0.40)

# --- (a) the two covariates are largely independent ------------------------
ax = fig.add_subplot(gs[0, 0])
ax.axhline(0, lw=0.6, color=C["rule"]); ax.axvline(0, lw=0.6, color=C["rule"])
ax.scatter(dC1, dD1, s=10, color=C["pt"], edgecolor="white", lw=0.35, zorder=3)
b, a0 = np.polyfit(dC1, dD1, 1)
xs = np.linspace(dC1.min(), dC1.max(), 20)
ax.plot(xs, a0 + b * xs, lw=1.1, color=C["ink"], zorder=4)
ax.set_xlabel("ΔC1, non-parenchymal content")
ax.set_ylabel("ΔD1, hepatocyte identity")
_lab = (f"$r$ = {A['pearson_r']:+.3f}\n$R^2$ = {A['r2_dD1_on_dC1']:.3f}\n"
        f"{D['n_pairs']} pairs").replace("-", "\u2212")
ax.text(0.03, 0.05, _lab,
        transform=ax.transAxes, fontsize=7.2, color=C["ink"], va="bottom")
ax.set_title("The covariates carry separate information", fontsize=7.2,
             pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (b) the two retention distributions -----------------------------------
ax = fig.add_subplot(gs[0, 1])
bins = np.arange(0, 2.5, 0.1)
ax.hist(rc, bins=bins, color=C["co"], alpha=0.62, edgecolor="white", lw=0.4,
        label="composition-adjusted retention", zorder=3)
ax.hist(rd, bins=bins, color=C["id"], alpha=0.62, edgecolor="white", lw=0.4,
        label="identity-adjusted retention", zorder=3)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.8, color=C["ink"], zorder=4)
ax.set_xlim(0, 2.0)
n_hi = int(sum(1 for v in list(rd) + list(rc) if v > 2.0))
ax.text(1.98, ax.get_ylim()[1] * 0.62, f"axis trimmed at 2.0;\n{n_hi} values lie above it",
        fontsize=7.2, color=C["muted"], ha="right", va="top", linespacing=1.3)
ax.set_xlabel("Retention fraction, single-covariate model")
ax.set_ylabel("Published liver signatures")
ax.legend(frameon=False, fontsize=7.2, loc="upper left", handlelength=0.85,
          borderpad=0.1, labelspacing=0.22, handletextpad=0.4,
          bbox_to_anchor=(-0.015, 1.02))
ax.text(0.985, 0.985,
        f"below 50%\ncomposition {D['retention_composition_only']['n_below_50pct']}"
        f"\nidentity {D['retention_identity_only']['n_below_50pct']}",
        transform=ax.transAxes, fontsize=7.2, color=C["ink"], ha="right",
        va="top")
ax.set_title(f"medians {D['retention_composition_only']['median']:.3f} "
             f"and {D['retention_identity_only']['median']:.3f}",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# --- (c) per signature, one against the other ------------------------------
ax = fig.add_subplot(gs[0, 2])
lim = 1.75
ax.plot([0, lim], [0, lim], lw=0.8, color=C["rule"], zorder=1)
ax.axvline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.55, zorder=2)
ax.axhline(0.5, ls=(0, (3, 2)), lw=0.7, color=C["ink"], alpha=0.55, zorder=2)
disc = (rc >= 0.90) & (rd < 0.50)
ax.scatter(rc[~disc], rd[~disc], s=11, color=C["pt"], edgecolor="white",
           lw=0.35, zorder=3)
ax.scatter(rc[disc], rd[disc], s=13, color=C["bad"], edgecolor="white",
           lw=0.35, zorder=4)
for mod, col, lab in [("REDUCTION", C["red"], "REDUCTION"),
                      ("DRAIN", C["dra"], "DRAIN")]:
    m = D["own_modules"][mod]
    ax.scatter([m["retention_C1"]], [m["retention_D1"]], s=42, marker="D",
               color=col, edgecolor="white", lw=0.7, zorder=6)
    dx, dy = ((-0.10, 0.16) if mod == "REDUCTION" else (-0.42, -0.02))
    ax.annotate(lab, xy=(m["retention_C1"], m["retention_D1"]),
                xytext=(m["retention_C1"] + dx, m["retention_D1"] + dy),
                fontsize=7.2, fontweight="bold", color=col, ha="right",
                va="center",
                arrowprops=dict(arrowstyle="-", lw=0.7, color=col,
                                shrinkA=1, shrinkB=4))
ax.text(0.02, 1.76, f"{int(disc.sum())} signatures kept by\ncomposition "
        f"adjustment,\nremoved by identity", fontsize=7.2, color=C["bad"],
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=1.5),
        ha="left", va="top")
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Retention after composition adjustment")
ax.set_ylabel("Retention after identity adjustment")
ax.set_title(f"ρ = {AG['spearman_rho']:+.3f}".replace("-", "\u2212") + "; "
             f"{AG['n_reclassified_at_50pct']} of {D['n_evaluable']} switch side",
             fontsize=7.2, pad=4, color=C["muted"])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

for xx, L in [(0.004, "a"), (0.337, "b"), (0.670, "c")]:
    fig.text(xx, 0.985, L, fontsize=10, fontweight="bold", color=C["ink"],
             ha="left", va="top")

stem = "/home/claude/Figure6_covariate_comparison"
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
print("written", im.size, "| discordant", int(disc.sum()),
      "| reclass", AG["n_reclassified_at_50pct"])
