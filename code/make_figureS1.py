"""Supplementary Figure S1: how the comparator panels were enumerated.
Every count is read from results/ENUMERATION_FLOW.json."""
import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image
plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.6,
    "svg.fonttype":"none","pdf.fonttype":42})
MM=1/25.4
C={"ink":"#0b0b0b","muted":"#898781","rule":"#c3c2b7","box":"#eef1f6",
   "liv":"#4a3aa7","lun":"#199e70","kid":"#d95926","out":"#f4e6e0","bad":"#d03b3b"}
F=json.load(open("/home/claude/rsi/results/ENUMERATION_FLOW.json"))

fig=plt.figure(figsize=(174*MM,178*MM))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

def box(x,y,w,h,text,fc=C["box"],ec=C["ink"],fs=7.6,bold=False,tc=None):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,
        boxstyle="round,pad=0,rounding_size=0.012",lw=0.9,edgecolor=ec,facecolor=fc))
    ax.text(x,y,text,ha="center",va="center",fontsize=fs,color=tc or C["ink"],
            fontweight="bold" if bold else "normal",linespacing=1.4)

def arrow(x0,y0,x1,y1,col=None):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",
        mutation_scale=8,lw=1.0,color=col or C["ink"]))

ax.text(0.5,0.975,"Enumeration of the comparator panels",ha="center",
        fontsize=9.2,fontweight="bold")
ax.text(0.5,0.945,"Rules fixed in §6t and applied unchanged to each tissue in §6x and §6y",
        ha="center",fontsize=7.2,color=C["muted"])

box(0.5,0.885,0.56,0.055,
    f"MSigDB human C2:CGP, release v2026.1.Hs\n{F['c2_cgp_total']:,} gene sets",bold=True)
arrow(0.5,0.857,0.5,0.828)
box(0.5,0.800,0.70,0.048,
    "Keyword rule, fixed before any set was scored: the standard name contains\n"
    "LIVER or HEPAT · LUNG, PULMONARY or ALVEOLAR · KIDNEY or RENAL",fs=7.0)

COLS={"liver":(0.17,C["liv"],"Liver"),"lung":(0.48,C["lun"],"Lung"),
      "kidney":(0.79,C["kid"],"Kidney")}
for t,(x,col,lab) in COLS.items():
    d=F[t]
    arrow(0.5,0.776,x,0.742,col)
    box(x,0.712,0.26,0.050,f"{lab}\n{d['keyword']} sets enumerated",
        fc="white",ec=col,bold=True,tc=col)
    def drop(y,n,why):
        ax.text(x+0.020,y,f"\u2212{n}  {why}",fontsize=7.2,color=C["bad"],
                ha="left",va="center")
    arrow(x,0.687,x,0.652,col); drop(0.6695,d['excluded_by_size'],"set size")
    box(x,0.622,0.26,0.050,f"{d['size_eligible']} of size 15 to 1000",fc="white",ec=col)
    arrow(x,0.597,x,0.562,col); drop(0.5795,d['not_on_platform'],"fewer than 10 symbols")
    box(x,0.532,0.26,0.050,
        f"{d['size_eligible']-d['not_on_platform']} evaluable on the platform",fc="white",ec=col)
    arrow(x,0.507,x,0.472,col); drop(0.4895,d['null_effect'],"no unadjusted effect")
    box(x,0.440,0.26,0.056,f"{d['evaluable']} with a tumor-adjacent\neffect to retain",
        fc="white",ec=col,bold=True,tc=col)
    arrow(x,0.412,x,0.378,col)
    box(x,0.348,0.26,0.048,f"{d['evaluable']} analyzed",fc=col,ec=col,bold=True,tc="white")

ax.text(0.5,0.302,"Rules applied, in this order and fixed before any retention value existed",
        ha="center",fontsize=7.6,fontweight="bold")
rules=[
 "1. Keyword. The set's MSigDB standard name contains the tissue keyword. No set is added or removed by hand.",
 "2. Size. Sets of fewer than 15 or more than 1000 symbols are excluded, because the score is a mean of gene-wise z scores.",
 "3. Platform. A set needs at least 10 of its symbols present in the cohort's matrix to be scored.",
 "4. Null effect. A set whose unadjusted paired tumor-adjacent shift does not reach two-sided Wilcoxon P < 0.05 has",
 "     no effect to retain and is excluded from the retention distribution. Its number is reported.",
 "",
 "No set was excluded on the basis of its retention value, and none was added after any result was seen. The source",
 "collection is deposited in full with the analysis code, so the whole enumeration can be reproduced from it.",
 "Every enumerated set and its inclusion status is listed in Supplementary Table S8.",
 "",
 "A pre-submission audit found that six of the 36 sets analyzed in lung are not lung sets: five matched on the",
 "surname McClung and one on the word alveolar. They were kept because the rule was fixed in advance. \u00a76ac",
 "reports the sensitivity with them removed, which leaves the lung distribution unchanged and reverses one comparison.",
]
for i,r in enumerate(rules):
    ax.text(0.055,0.272-i*0.0205,r,ha="left",fontsize=7.2,
            color=C["muted"] if i>=6 else C["ink"])

stem="/home/claude/FigureS1_enumeration"
fig.savefig(stem+".png",dpi=300,facecolor="white"); fig.savefig(stem+".svg",facecolor="white")
plt.close(fig)
im=Image.open(stem+".png").convert("RGB")
im.save(stem+".png",dpi=(300,300)); im.save(stem+".tiff",compression="tiff_lzw",dpi=(300,300))
print("written",im.size)
