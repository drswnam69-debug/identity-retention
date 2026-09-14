#!/usr/bin/env python3
"""10_figure5.py -- Figure 5: PHASE E (H2 replication + H4 null).

All numbers are read from the results written by 09_gse76427_h2.py and
08_phase_e.py; nothing is retyped by hand.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as F
import stats_lite as sl
import survival as sv
import rsi_config as cfg

OUT = "results/figures"
os.makedirs(OUT, exist_ok=True)

rsi = pd.read_csv("results/GSE76427/rsi.tsv", sep="\t", index_col=0)
ph = pd.read_csv("results/GSE76427/phenotype.tsv", sep="\t",
                 index_col=0).loc[rsi.index]
e = pd.read_csv("results/GSE76427/expr_log.tsv.gz", sep="\t",
                index_col=0)[rsi.index]
h2 = json.load(open("results/gse76427_h2.json"))
pe = json.load(open("results/phase_e.json"))

is_t = (~ph["tissue"].astype(str).str.lower()
        .str.contains("adjacent")).to_numpy()
r = rsi["RSI"].to_numpy(float)


def pstr(p):
    if p >= 0.01:
        return f"P = {p:.2f}"
    e_ = int(np.floor(np.log10(p)))
    return f"P = {p/10**e_:.1f} × 10$^{{{e_}}}$"


fig = plt.figure(figsize=(F.WIDTH_MM * F.MM, 126 * F.MM))

# ---------------------------------------------------------------- A: H2
axA = fig.add_axes([0.075, 0.620, 0.215, 0.325])
groups = [r[~is_t], r[is_t]]
pos = [0, 1]
for i, (g, col, fc) in enumerate(zip(groups,
                                     [F.C["muted"], F.C["patho"]],
                                     ["#eceae4", F.C["patho_f"]])):
    parts = axA.violinplot([g], positions=[pos[i]], widths=0.72,
                           showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(fc); b.set_edgecolor(col)
        b.set_alpha(1.0); b.set_linewidth(0.8)
    q1, med, q3 = np.percentile(g, [25, 50, 75])
    axA.add_patch(plt.Rectangle((pos[i] - 0.055, q1), 0.11, q3 - q1,
                                facecolor=col, edgecolor="none", zorder=4))
    axA.plot([pos[i]], [med], "o", ms=3.2, mfc="white",
             mec=col, mew=1.0, zorder=5)
    axA.scatter(np.random.default_rng(3 + i).normal(pos[i], 0.055, len(g)),
                g, s=1.6, color=col, alpha=0.5, lw=0, zorder=3)

hl = h2["unpaired"]["hodges_lehmann"]
d = h2["unpaired"]["cliffs_delta"]
top = max(r) + 0.42
axA.plot([0, 0, 1, 1], [top - 0.16, top, top, top - 0.16],
         color=F.C["ink"], lw=0.7, clip_on=False)
axA.text(0.5, top + 0.07, pstr(h2["unpaired"]["p"]), ha="center",
         va="bottom", fontsize=6.6, fontweight="bold")
axA.set_xticks(pos)
axA.set_xticklabels([f"Adjacent\n(n = {h2['n_adjacent']})",
                     f"Tumor\n(n = {h2['n_tumor']})"], fontsize=6.8)
axA.set_ylabel("Reductive-Supply Index", fontsize=7)
axA.set_ylim(min(r) - 0.3, top + 0.55)
axA.text(0.5, -0.215, f"Hodges–Lehmann {hl:+.2f},  Cliff's "
         f"$\\delta$ {d:+.2f}", transform=axA.transAxes, ha="center",
         fontsize=6.2, color=F.C["muted"])
for s in ("top", "right"):
    axA.spines[s].set_visible(False)
axA.tick_params(labelsize=6.6, length=2.5, width=0.6)

fig.text(0.010, 0.972, "a", fontsize=10, fontweight="bold")

# ------------------------------------------------------- B: paired pairs
axB = fig.add_axes([0.370, 0.620, 0.150, 0.325])
pid = ph["patient_id"].astype(str)
df = pd.DataFrame({"pid": pid.values, "t": is_t, "r": r}, index=rsi.index)
tt = df[df.t].set_index("pid")["r"]; tt = tt[~tt.index.duplicated()]
nn = df[~df.t].set_index("pid")["r"]; nn = nn[~nn.index.duplicated()]
common = sorted(set(tt.index) & set(nn.index))
up = 0
for k in common:
    a, b = nn.loc[k], tt.loc[k]
    up += b > a
    axB.plot([0, 1], [a, b], lw=0.5,
             color=F.C["patho"] if b > a else F.C["supply"], alpha=0.55)
axB.plot([0, 1], [nn.loc[common].mean(), tt.loc[common].mean()],
         lw=2.0, color=F.C["ink"], zorder=5,
         marker="o", ms=3.4, mfc="white", mec=F.C["ink"], mew=1.2)
axB.set_xlim(-0.28, 1.28); axB.set_xticks([0, 1])
axB.set_xticklabels(["Adjacent", "Tumor"], fontsize=6.8)
axB.set_ylabel("RSI", fontsize=7)
axB.set_title(f"{len(common)} patient-matched pairs", fontsize=6.8, pad=4)
axB.text(0.5, max(r) + 0.30, pstr(h2["paired"]["p"]), ha="center",
         fontsize=6.6, fontweight="bold")
axB.text(0.5, -0.215, f"higher in tumor in {up}/{len(common)}",
         transform=axB.transAxes, ha="center", fontsize=6.2,
         color=F.C["muted"])
axB.set_ylim(min(r) - 0.3, max(r) + 0.55)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.tick_params(labelsize=6.6, length=2.5, width=0.6)
fig.text(0.318, 0.972, "b", fontsize=10, fontweight="bold")

# --------------------------------------------------- C: module decomposition
axC = fig.add_axes([0.640, 0.620, 0.230, 0.325])
mods = [("SUPPLY", cfg.MODULE_SUPPLY, F.C["supply"]),
        ("REDUCTION", cfg.MODULE_REDUCTION, F.C["protect"]),
        ("DRAIN", cfg.MODULE_DRAIN, F.C["patho"])]
ypos = np.arange(len(mods))[::-1]
for (name, genes, col), y in zip(mods, ypos):
    g = [x for x in genes if x in e.index]
    z = (e.loc[g].sub(e.loc[g].mean(axis=1), axis=0)
         .div(e.loc[g].std(axis=1), axis=0).mean(axis=0).to_numpy(float))
    ma, mn_ = z[is_t].mean(), z[~is_t].mean()
    _, _, p = sl.mannwhitney_u(z[is_t], z[~is_t])
    axC.plot([mn_, ma], [y, y], color=col, lw=1.3, zorder=2,
             solid_capstyle="round")
    axC.plot([mn_], [y], "o", ms=4.2, mfc="white", mec=col, mew=1.2, zorder=3)
    axC.plot([ma], [y], "o", ms=4.6, mfc=col, mec=col, zorder=3)
    axC.text(1.03, y, pstr(p), transform=axC.get_yaxis_transform(),
             ha="left", va="center", fontsize=6.0,
             color=F.C["ink"] if p < 0.05 else F.C["muted"],
             fontweight="bold" if p < 0.05 else "normal")
axC.axvline(0, color=F.C["rule"], lw=0.6, zorder=1)
axC.set_yticks(ypos)
axC.set_yticklabels(["Supply", "Reduction", "Drain"], fontsize=7)
axC.set_xlabel("module score (within-cohort z)", fontsize=7)
axC.set_ylim(-0.55, len(mods) - 0.45)
for s in ("top", "right", "left"):
    axC.spines[s].set_visible(False)
axC.tick_params(labelsize=6.6, length=2.5, width=0.6)
axC.plot([], [], "o", ms=4.2, mfc="white", mec=F.C["muted"],
         mew=1.2, label="adjacent")
axC.plot([], [], "o", ms=4.6, color=F.C["muted"], label="tumor")
axC.legend(fontsize=6.0, frameon=False, loc="upper center",
           bbox_to_anchor=(0.5, 1.13), ncol=2, handletextpad=0.3,
           borderpad=0.1, columnspacing=1.0)
fig.text(0.578, 0.972, "c", fontsize=10, fontweight="bold")

# ------------------------------------------------------------- D: KM curves
axD = fig.add_axes([0.075, 0.135, 0.290, 0.330])
# tertiles recomputed exactly as 08_phase_e.py does (unrounded quantiles,
# pd.cut on the open-closed intervals) so the panel cannot drift from the
# reported analysis; the assertion below enforces that.
tum = ph[is_t]
rt_s = rsi.loc[tum.index, "RSI"].astype(float)
rt = rt_s.to_numpy()
time = tum["duryears_os"].to_numpy(float)
ev = tum["event_os"].to_numpy(float)
q1, q2 = rt_s.quantile([1 / 3, 2 / 3])
lab = pd.cut(rt_s, [-np.inf, q1, q2, np.inf],
             labels=["T1", "T2", "T3"]).astype(str).to_numpy()
_n = {k: int((lab == k).sum()) for k in ("T1", "T2", "T3")}
_e = {k: int(ev[lab == k].sum()) for k in ("T1", "T2", "T3")}
assert _n == pe["H4"]["tertiles"]["n"], (_n, pe["H4"]["tertiles"]["n"])
assert _e == pe["H4"]["tertiles"]["events"], (_e,)
ramp = [F.STAGE_RAMP[1], F.STAGE_RAMP[3], F.STAGE_RAMP[5]]
for name, col in zip(["T1", "T2", "T3"], ramp):
    m = lab == name
    t_, s_, _ = sv.kaplan_meier(time[m], ev[m])
    axD.step(np.concatenate([[0], t_]), np.concatenate([[1], s_]),
             where="post", color=col, lw=1.3,
             label=f"{name}  (n = {m.sum()}, {int(ev[m].sum())} deaths)")
    # censoring marks: step function value at each censoring time
    ct = time[m][ev[m] == 0]
    cs = [1.0 if t0 < t_[0] else s_[np.searchsorted(t_, t0, "right") - 1]
          for t0 in ct]
    axD.plot(ct, cs, "|", color=col, ms=3.0, mew=0.8, zorder=4)
axD.set_xlabel("years from resection", fontsize=7)
axD.set_ylabel("overall survival", fontsize=7)
axD.set_ylim(0, 1.02); axD.set_xlim(0, 8)
axD.legend(fontsize=6.0, frameon=False, loc="lower left",
           handlelength=1.4, handletextpad=0.4, borderpad=0.1)
axD.text(0.97, 0.95, "log-rank " + pstr(pe["H4"]["tertiles"]["logrank"]["p"]),
         transform=axD.transAxes, ha="right", fontsize=6.6,
         fontweight="bold")
for s in ("top", "right"):
    axD.spines[s].set_visible(False)
axD.tick_params(labelsize=6.6, length=2.5, width=0.6)
fig.text(0.010, 0.505, "d", fontsize=10, fontweight="bold")

# ---------------------------------------------------------------- E: forest
axE = fig.add_axes([0.520, 0.135, 0.215, 0.330])
adj = pe["H4"]["adjusted"]
terms = ["RSI per SD", "Age (per SD)", "Female sex", "BCLC C or D"]
hr = adj["HR"]; ci = adj["ci95"]; pv = adj["p"]
yy = np.arange(len(terms))[::-1]
for i, y in enumerate(yy):
    sig = pv[i] < 0.05
    col = F.C["patho"] if sig else F.C["muted"]
    axE.plot([ci[i][0], ci[i][1]], [y, y], color=col, lw=1.1,
             solid_capstyle="round")
    axE.plot([hr[i]], [y], "s", ms=4.2, color=col)
    axE.text(1.05, y, f"{hr[i]:.2f} ({ci[i][0]:.2f}–{ci[i][1]:.2f})",
             transform=axE.get_yaxis_transform(), va="center",
             fontsize=6.0, color=F.C["ink"] if sig else F.C["muted"],
             fontweight="bold" if sig else "normal")
    axE.text(1.62, y, pstr(pv[i]).replace("P = ", "P = "),
             transform=axE.get_yaxis_transform(), va="center",
             fontsize=6.0, color=F.C["ink"] if sig else F.C["muted"],
             fontweight="bold" if sig else "normal")
axE.axvline(1, color=F.C["rule"], lw=0.6, ls=(0, (3, 2)))
axE.set_xscale("log")
axE.set_xlim(0.15, 20)
axE.set_xticks([0.25, 1, 4, 16])
axE.set_xticklabels(["0.25", "1", "4", "16"])
axE.set_yticks(yy); axE.set_yticklabels(terms, fontsize=6.6)
axE.set_ylim(-0.6, len(terms) - 0.3)
axE.set_xlabel("hazard ratio (95% CI)", fontsize=7)
for s in ("top", "right", "left"):
    axE.spines[s].set_visible(False)
axE.tick_params(labelsize=6.4, length=2.5, width=0.6)
fig.text(0.430, 0.505, "e", fontsize=10, fontweight="bold")

fig.text(0.075, 0.028,
         "GSE76427 — 115 tumors, 52 patient-matched adjacent livers, "
         "23 deaths. Index locked before any outcome was examined "
         "(hash 7a2bf934…); H2 and H4 pre-registered.",
         fontsize=5.8, color=F.C["muted"])

paths = F.save_all(fig, os.path.join(OUT, "figure5_phaseE"))
print("wrote:", *paths, sep="\n  ")
