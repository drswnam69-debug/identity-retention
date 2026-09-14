"""Figure 11: what the source proteogenomic study settles (PREREG 6aa)."""
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
plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.2,"svg.fonttype":"none",
    "pdf.fonttype":42,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
    "xtick.major.size":2.5,"ytick.major.size":2.5})
MM=1/25.4
C={"ink":"#0b0b0b","muted":"#898781","rule":"#c3c2b7","pt":"#8d9bb5","id":"#4a3aa7",
   "co":"#c99a12","red":"#199e70","dra":"#d95926","bad":"#d03b3b"}
R=f"{IR_ROOT}"
G=json.load(open(f"{R}/results/GAO_MRNA_PROTEIN_6aa.json"))
H=json.load(open(f"{R}/results/GAO_HE_PURITY_6aa.json"))
corr=pd.read_csv(f"{R}/data/gao_mrna_protein_corr.tsv",sep="\t")
corr["spearman"]=pd.to_numeric(corr["spearman"],errors="coerce")
allc=corr["spearman"].dropna().to_numpy()

fig=plt.figure(figsize=(174*MM,60*MM))
gs=gridspec.GridSpec(1,3,figure=fig,left=0.062,right=0.988,top=0.845,bottom=0.20,wspace=0.36)

# (a) mRNA-protein concordance of the study's own genes, against all genes
ax=fig.add_subplot(gs[0,0])
ax.hist(allc,bins=np.arange(-0.25,1.0,0.025),color=C["pt"],edgecolor="none",zorder=2)
own=[("MTARC1",C["dra"]),("POR",C["dra"]),("MTARC2",C["dra"]),
     ("CYB5R3",C["red"]),("AIFM2",C["red"]),("CYB5R1",C["red"])]
ymax=ax.get_ylim()[1]
for g,col in own:
    ax.plot([G["own_genes"][g]]*2,[0,ymax*0.86],lw=1.3,color=col,zorder=4)
ax.axvline(np.median(allc),ls=(0,(3,2)),lw=0.9,color=C["ink"],zorder=3)
ax.text(0.985,0.985,f"all genes\nmedian {np.median(allc):.3f}",
        transform=ax.transAxes,fontsize=7.2,color=C["ink"],ha="right",va="top",
        linespacing=1.25)
ax.text(0.03,0.985,"DRAIN\n  MTARC1 0.831\n  POR 0.785\n  MTARC2 0.762",
        transform=ax.transAxes,fontsize=7.2,color=C["dra"],va="top",linespacing=1.3)
ax.text(0.03,0.640,"REDUCTION\n  CYB5R3 0.639\n  AIFM2 0.580\n  CYB5R1 0.541",
        transform=ax.transAxes,fontsize=7.2,color=C["red"],va="top",linespacing=1.3)
ax.set_ylim(0,ymax*1.12)
ax.set_xlabel("mRNA to protein Spearman ρ, same 159 patients")
ax.set_ylabel(f"Genes (n = {len(allc)})")
ax.set_title("DRAIN genes agree best between levels",fontsize=7.2,pad=4,color=C["muted"])
for s in ("top","right"): ax.spines[s].set_visible(False)

# (b) the registered test
ax=fig.add_subplot(gs[0,1])
S=pd.DataFrame(G["signatures"])
ax.scatter(S["mean_mrna_protein_rho"],S["discrepancy"],s=13,color=C["pt"],
           edgecolor="white",lw=0.35,zorder=3)
b,a0=np.polyfit(S["mean_mrna_protein_rho"],S["discrepancy"],1)
xs=np.linspace(S["mean_mrna_protein_rho"].min(),S["mean_mrna_protein_rho"].max(),20)
ax.plot(xs,a0+b*xs,lw=1.1,color=C["ink"],zorder=4)
ax.set_xlabel("Mean mRNA to protein ρ of the signature")
ax.set_ylabel("|retention, transcript \u2212 protein|")
ax.text(0.03,0.97,f"ρ = {G['spearman_rho']:+.3f}  (P = {G['spearman_p']:.2f})\n"
        f"§6aa predicted ρ ≤ −0.30\nprediction NOT supported",
        transform=ax.transAxes,fontsize=7.2,color=C["bad"],va="top")
ax.set_title(f"{G['n_signatures_both_levels']} signatures scored at both levels",
             fontsize=7.2,pad=4,color=C["muted"])
for s in ("top","right"): ax.spines[s].set_visible(False)

# (c) histological purity against the two covariates
ax=fig.add_subplot(gs[0,2])
xs=["C1\ncomposition","D1\nidentity"]
vals=[H["purity_vs_C1"]["rho"],H["purity_vs_D1"]["rho"]]
ps=[H["purity_vs_C1"]["p"],H["purity_vs_D1"]["p"]]
cols=[C["co"],C["id"]]
ax.bar(np.arange(2),vals,width=0.5,color=cols,edgecolor="none",zorder=3)
ax.axhline(0,lw=0.7,color=C["ink"])
for i,(v,pv) in enumerate(zip(vals,ps)):
    ax.text(i,v+(0.035 if v>0 else -0.035),
            (f"\u03c1 = {v:+.3f}\n$P$ = {pv:.2g}" if pv >= 1e-3 else
             "\u03c1 = %+.3f\n$P$ = %.1f \u00d7 10$^{%d}$"
             % (v, float(f"{pv:.1e}".split("e")[0]), int(f"{pv:.1e}".split("e")[1]))
             ).replace("-", "\u2212").replace("$^{\u2212", "$^{-"),
            ha="center",va="bottom" if v>0 else "top",fontsize=7.2,color=C["ink"])
ax.set_xticks(np.arange(2)); ax.set_xticklabels(xs,fontsize=7.2)
ax.set_ylim(-0.80,0.46)
ax.set_ylabel("Spearman ρ with tumor purity\nread from H&E by a pathologist")
ax.text(0.02,0.985,f"{H['n_pairs_with_purity']} cases\npurity 0.50 to 0.99",
        transform=ax.transAxes,fontsize=7.2,color=C["muted"],ha="left",va="top",
        linespacing=1.25)
ax.set_title("Purity tracks the composition covariate,\nnot the identity covariate",
             fontsize=7.2, pad=4, color=C["muted"], linespacing=1.3)
for s in ("top","right"): ax.spines[s].set_visible(False)

for xx,L in [(0.004,"a"),(0.340,"b"),(0.672,"c")]:
    fig.text(xx,0.985,L,fontsize=10,fontweight="bold",color=C["ink"],ha="left",va="top")
stem=f"{IR_DOCS}/Figure11_gao_proteogenomic"
fig.savefig(stem+".png",dpi=300,facecolor="white"); fig.savefig(stem+".svg",facecolor="white")
plt.close(fig)
im=Image.open(stem+".png").convert("RGB")
im.save(stem+".png",dpi=(300,300)); im.save(stem+".tiff",compression="tiff_lzw",dpi=(300,300))
print("written",im.size)
