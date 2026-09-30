#!/usr/bin/env python3
"""82_build_MCP.py -- assemble submission_MCP/ deterministically.

Every file the editor receives is produced here from a named source, and the
built DOCX is re-read and searched for markers unique to the current text, so a
stale artifact cannot pass silently. Written after a submission once went out
carrying a stale build.
"""
import hashlib, os, re, shutil, subprocess, sys, zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = f"{ROOT}/manuscript/manuscript_MCP.md"
FIGS = f"{ROOT}/figures/mcp"
PKG = "/home/claude/submission_MCP"
FIG_PNG = {1: "Figure1_panel_comes_apart.png",
           2: "Figure2_eight_datasets.png",
           3: "Figure3_cholangiocarcinoma.png"}
MARKERS = [
    "Albumin is not a hepatocyte identity marker",
    "764 patient-matched pairs", "cholangiocytes do not synthesize albumin".capitalize(),
    "a surviving share of at least 0.95", "the 14 proteins under test",
    "+8.71 in adjacent tissue", "1.2th and 1.6th percentile",
    "Spearman rho = \u22120.084", "abundance accounts for 17.5%",
    "gene-disjoint subfamily contains only 2 sets",
    "Experimental Design and Statistical Rationale",
    "Werner, T. et al.", "MassIVE repository under MSV000095336",
    "80 tumor and 69 adjacent samples from 85 patients",
    "not a third-party registration",
]


def run(*a):
    r = subprocess.run(a, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"FAILED: {' '.join(a)}\n{r.stderr[:2000]}")
    return r


def docx_text(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    t = re.sub(r"<[^>]+>", "", re.sub(r"</w:p>", "\n", xml))
    import html
    return html.unescape(t)


def superscript_runs(path):
    with zipfile.ZipFile(path) as z:
        return z.read("word/document.xml").decode("utf-8").count("w:vertAlign")


def tighten(path, twips=1080):
    """narrower margins so the front matter and legends sit comfortably"""
    import shutil as sh
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for it in zin.infolist():
            data = zin.read(it.filename)
            if it.filename == "word/document.xml":
                x = data.decode("utf-8")
                x = re.sub(r'<w:pgMar [^/]*/>',
                           f'<w:pgMar w:top="{twips}" w:right="{twips}" '
                           f'w:bottom="{twips}" w:left="{twips}" w:header="720" '
                           f'w:footer="720" w:gutter="0"/>', x)
                data = x.encode("utf-8")
            zout.writestr(it, data)
    sh.move(tmp, path)


def main() -> int:
    os.makedirs(PKG, exist_ok=True)
    md = open(SRC, encoding="utf-8").read().rstrip() + "\n\n"
    for n, png in FIG_PNG.items():
        marker = f"**Figure {n}."
        i = md.index(marker)
        end = md.index("\n\n", i)
        md = md[:end] + f"\n\n![]({FIGS}/{png}){{width=16cm}}" + md[end:]
    tmp = "/tmp/manuscript_MCP_build.md"
    open(tmp, "w", encoding="utf-8").write(md)

    out = f"{PKG}/01_Manuscript.docx"
    run("pandoc", tmp, "-o", out, "--standalone")
    tighten(out)
    run("pandoc", tmp, "-o", f"{PKG}/02_Manuscript.pdf",
        "--pdf-engine=xelatex", "-V", "geometry:margin=2cm",
        "-V", "mainfont=DejaVu Sans")

    txt = docx_text(out)
    missing = [m for m in MARKERS if m not in txt and m.lower() not in txt.lower()]
    if missing:
        sys.exit(f"STALE BUILD: markers absent from the DOCX: {missing}")
    runs = superscript_runs(out)
    if runs < 3:
        sys.exit(f"the affiliation markers are not superscripted ({runs} runs)")
    with zipfile.ZipFile(out) as z:
        n_media = len([x for x in z.namelist() if x.startswith("word/media/")])
    if n_media != 3:
        sys.exit(f"expected 3 embedded figures, found {n_media}")
    print(f"  {runs} superscript runs, {n_media} figures embedded, "
          f"every current-text marker present")

    for n, png in FIG_PNG.items():
        shutil.copy(f"{FIGS}/{png.replace('.png', '.tiff')}",
                    f"{PKG}/Figure{n}.tiff")
    for f in ("MANUSCRIPT_NUMBERS.json", "MANUSCRIPT_EXTRAS.json",
              "HOLDOUT_Yi2023_EncyclopeDIA.json", "ICCA_VALIDATION.json"):
        shutil.copy(f"{ROOT}/results/{f}", f"{PKG}/{f}")
    for f in ("70_panel_coherence.py", "77_holdout_validation.py",
              "78_icca_validation.py", "79_manuscript_numbers.py",
              "79b_manuscript_extras.py", "80_verify_MCP.py",
              "81_figures_MCP.py"):
        shutil.copy(f"{ROOT}/code/{f}", f"{PKG}/{f}")
    shutil.copy(SRC, f"{PKG}/manuscript_MCP.md")

    print(f"\n  {PKG}")
    for f in sorted(os.listdir(PKG)):
        p = os.path.join(PKG, f)
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
        print(f"  {os.path.getsize(p):>10,}  {h}  {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
