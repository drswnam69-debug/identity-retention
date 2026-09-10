"""Shared style and drawing helpers for the concept figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib import font_manager as fm

# Liberation Sans is metric-compatible with Arial; substitute true Arial at
# production if the journal requires the exact face.
for f in ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf"):
    fm.fontManager.addfont(f)

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 7.2,
    "svg.fonttype": "none",      # keep text editable in the SVG
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})

MM = 1 / 25.4
WIDTH_MM = 174            # journal single-figure full width

# Validated with the dataviz palette validator (light surface, categorical):
#   node scripts/validate_palette.js "#4a3aa7,#199e70,#d95926" --mode light
#   -> PASS on all six checks (worst adjacent CVD dE 9.4 deutan / 21.8 tritan,
#      normal-vision dE 26.5, all >= 3:1 contrast).
# The three categorical roles are supply / protective / pathogenic. Outcome and
# hub are neutral inks, deliberately NOT a fourth series color.
C = {
    "supply":   "#4a3aa7",   # mevalonate supply      (violet)
    "supply_f": "#e3e0f4",
    "protect":  "#199e70",   # protective arm         (green-teal)
    "protect_f": "#d7f0e6",
    "patho":    "#d95926",   # pathogenic arm         (orange)
    "patho_f":  "#fbe4d8",
    "hub":      "#3A3A3A",
    "hub_f":    "#EDEDEB",
    "block":    "#d03b3b",   # inhibition             (status critical)
    "out":      "#0b0b0b",   # convergence / outcome  (terminal node, not a series)
    "out_f":    "#ffffff",
    "ink":      "#0b0b0b",
    "muted":    "#898781",
    "rule":     "#c3c2b7",
}

# Sequential ramp for the ordered MASLD spectrum (one hue, light -> dark)
STAGE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281"]


def box(ax, x, y, w, h, text, ec, fc, *, fs=7, bold=False, lw=0.9,
        radius=0.012, ha="center", tc=None, zorder=3):
    """Rounded box centred on (x, y) with wrapped label text."""
    ax.add_patch(FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=zorder))
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs,
            color=tc or C["ink"], zorder=zorder + 1,
            fontweight="bold" if bold else "normal", linespacing=1.35)


def arrow(ax, p0, p1, color, *, lw=1.0, style="-|>", rad=0.0, ls="-",
          ms=6, zorder=2, alpha=1.0):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle=style, mutation_scale=ms, linewidth=lw,
        color=color, linestyle=ls, zorder=zorder, alpha=alpha,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=1.5, shrinkB=1.5))


def blocked(ax, x, y, color, *, w=0.022, lw=1.6):
    """Perpendicular bar denoting pharmacological inhibition."""
    ax.plot([x, x], [y - w, y + w], color=color, lw=lw,
            solid_capstyle="butt", zorder=6)


def panel_label(ax, x, y, letter):
    ax.text(x, y, letter, fontsize=10, fontweight="bold",
            color=C["ink"], ha="left", va="top")


def blank_axes(fig, rect):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    return ax


def save_all(fig, stem):
    """Write PNG 300 dpi, TIFF LZW 300 dpi, and SVG, all flattened to RGB."""
    from PIL import Image
    fig.savefig(f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(f"{stem}.svg", facecolor="white")
    plt.close(fig)
    # journals expect RGB, not RGBA; flatten onto white and re-emit
    im = Image.open(f"{stem}.png")
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1] if im.mode == "RGBA" else None)
        im = bg
    im.save(f"{stem}.png", dpi=(300, 300))
    im.save(f"{stem}.tiff", compression="tiff_lzw", dpi=(300, 300))
    return [f"{stem}.png", f"{stem}.tiff", f"{stem}.svg"]
