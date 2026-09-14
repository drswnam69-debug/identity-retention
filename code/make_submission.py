#!/usr/bin/env python3
"""make_submission.py -- assemble submission_GigaScience/ deterministically.

The package used to be assembled by hand, one shell command at a time, which is
how a stale supplementary PDF and a set of file names the journal does not allow
both survived several audits. Everything the editor receives is built here, from
the manuscript and the repository, and nothing else is left in the folder.

GigaScience requires that "every supplementary material file contains the phrase
'supplementary material' as part of the actual file name", so the supplementary
names below carry it and 50_consistency_check.py verifies that they do.
"""
import hashlib
import os
import re
import shutil
import subprocess
import zipfile

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
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, CODE as IR_CODE, DOCS as IR_DOCS

HOME = IR_DOCS
PKG = f"{HOME}/submission_GigaScience"
SF2 = f"{HOME}/SF2build/SupplementaryFile2"
BUILD = f"{HOME}/build"

SF1_PDF = "03_Supplementary Material 1_Preregistered protocol.pdf"
SF2_XLSX = "04_Supplementary Material 2_Supplementary tables.xlsx"
SF3_PDF = "05_Supplementary Material 3_Supplementary Note.pdf"
FIGS1 = "Figure S1_Supplementary Material.tif"
# The code deposit is no longer a supplementary file. Editorial Manager unpacks
# an uploaded archive into its members, so the 259-file code and result deposit
# arrived as 259 separate submission items. It now reaches the reader through
# the repository under its Zenodo digital object identifier, which is what the
# journal asks for anyway, and this build writes the archive for that upload
# next to the package rather than inside it.
ZENODO_ZIP = "identity-retention_archive.zip"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"{cmd[0]} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")


def pandoc(src, out, extra=()):
    run(["pandoc", src, "--from", "markdown+pipe_tables", *extra, "-o", out])


def zip_dir(root, out, arc_root):
    if os.path.exists(out):
        os.remove(out)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for base, dirs, files in os.walk(root):
            dirs[:] = [d for d in sorted(dirs) if d != "__pycache__"]
            for f in sorted(files):
                if f.endswith(".pyc"):
                    continue
                full = os.path.join(base, f)
                z.write(full, os.path.join(arc_root, os.path.relpath(full, root)))


def refresh_sf2():
    """Bring the deposited copy of the code back in line with the working tree."""
    shutil.rmtree(f"{SF2}/code", ignore_errors=True)
    shutil.copytree(IR_CODE, f"{SF2}/code",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    # The results were a stale snapshot: code/ was re-synced on every build and
    # results/ never was, so the deposit shipped a checker that raised
    # FileNotFoundError on files the manuscript cites. Both are re-synced now.
    shutil.rmtree(f"{SF2}/results", ignore_errors=True)
    shutil.copytree(IR_RESULTS, f"{SF2}/results",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy(f"{HOME}/SupplementaryFile1_v2.md", f"{SF2}/SupplementaryFile1_v2.md")
    # This used to be copied in after the archive had been zipped and
    # checksummed, so the release carried the previous build's note and
    # SHA256SUMS.txt did not cover the file that shipped.
    shutil.copy(f"{HOME}/SupplementaryNote_v1.md", f"{SF2}/SupplementaryNote_v1.md")
    # every top-level analysis script the manuscript points at
    for f in ("build_preprint.py", "sync_figures.py", "make_metadata.py",
              "renumber_refs.py", "make_submission.py",
              "make_supp_tables_xlsx.py", "insert_alt_text.py",
              "make_release.py"):
        src = os.path.join(HOME, f)
        if os.path.exists(src):
            shutil.copy(src, os.path.join(SF2, "code", f))
    for f in sorted(os.listdir(HOME)):
        if f.startswith("make_figure") and f.endswith(".py"):
            shutil.copy(os.path.join(HOME, f), os.path.join(SF2, "code", f))
    # checksums last, over everything else
    sums = os.path.join(SF2, "SHA256SUMS.txt")
    if os.path.exists(sums):
        os.remove(sums)
    lines = []
    for base, dirs, files in os.walk(SF2):
        dirs[:] = [d for d in sorted(dirs) if d != "__pycache__"]
        for f in sorted(files):
            if f.endswith(".pyc"):
                continue
            full = os.path.join(base, f)
            h = hashlib.sha256(open(full, "rb").read()).hexdigest()
            lines.append(f"{h}  {os.path.relpath(full, SF2)}")
    open(sums, "w").write("\n".join(sorted(lines, key=lambda l: l.split("  ", 1)[1])) + "\n")
    print(f"  SF2 refreshed, {len(lines)} files checksummed")


def figure_legends_docx(ms, out):
    block = ms[ms.index("\n## Figure legends\n"):ms.index("\n## Additional files\n")]
    tmp = f"{BUILD}/_legends.md"
    open(tmp, "w").write(block.strip() + "\n")
    pandoc(tmp, out)


def main():
    ms = open(f"{HOME}/manuscript_v4.md", encoding="utf-8").read()
    os.makedirs(PKG, exist_ok=True)
    os.makedirs(BUILD, exist_ok=True)

    # -- Supplementary Material 1: the protocol, as one PDF ------------------
    # It was a two-file ZIP holding the same document twice, in Markdown and in
    # PDF. The Markdown source stays in the repository archive, where a reader
    # who wants to diff it can find it, and the editor receives one file.
    pandoc(f"{HOME}/SupplementaryFile1_v2.md", f"{PKG}/{SF1_PDF}",
           extra=["--pdf-engine", "xelatex", "-V", "geometry:margin=2.4cm",
                  "-V", "mainfont=DejaVu Serif", "-V", "monofont=DejaVu Sans Mono",
                  "-V", "fontsize=10pt"])
    print(f"  {SF1_PDF}")

    # -- the repository archive, for Zenodo rather than for the journal ------
    # refresh_sf2 keeps the build tree current; make_release turns it into the
    # tree that is actually deposited, which is not the same set of files.
    refresh_sf2()
    run([_sys.executable, f"{HOME}/make_release.py"])
    zip_dir(f"{HOME}/release/identity-retention", f"{HOME}/{ZENODO_ZIP}",
            "identity-retention")
    print(f"  {ZENODO_ZIP} (upload target: GitHub and Zenodo, not the journal)")

    # -- Supplementary Material 2: the result tables, as one workbook --------
    run([_sys.executable, f"{HOME}/make_supp_tables_xlsx.py"])
    shutil.copy(f"{BUILD}/{SF2_XLSX}", f"{PKG}/{SF2_XLSX}")
    print(f"  {SF2_XLSX}")

    # -- manuscript, cover letter, metadata, legends -------------------------
    shutil.copy(f"{BUILD}/Nam_2026_manuscript_v4.docx", f"{PKG}/01_Manuscript.docx")
    shutil.copy(f"{BUILD}/Nam_2026_preprint_v4.pdf", f"{PKG}/02_Manuscript_typeset.pdf")
    pandoc(f"{HOME}/CoverLetter_GigaScience.md", f"{PKG}/00_CoverLetter.docx")
    pandoc(f"{HOME}/CoverLetter_GigaScience.md", f"{PKG}/00_CoverLetter.pdf",
           extra=["--pdf-engine", "xelatex", "-V", "geometry:margin=2.4cm",
                  "-V", "mainfont=DejaVu Serif"])
    # Supplementary Material 3: the analyses the main text reports in summary
    pandoc(f"{HOME}/SupplementaryNote_v1.md", f"{PKG}/{SF3_PDF}",
           extra=["--pdf-engine", "xelatex", "-V", "geometry:margin=2.4cm",
                  "-V", "mainfont=DejaVu Serif", "-V", "monofont=DejaVu Sans Mono",
                  "-V", "fontsize=10pt", "-V", "linestretch=1.15"])
    shutil.copy(f"{HOME}/Figure_alt_text.md", f"{PKG}/06_Figure_alt_text.md")
    figure_legends_docx(ms, f"{PKG}/07_Figure_legends.docx")

    # -- supplementary figure, under a compliant name ------------------------
    # sync_figures owns the figure copies, including the resolution that puts them
    # at the journal's 170 mm page width; do not overwrite its output with a
    # straight copy, which is how the supplementary figure went back to 174 mm.
    if not os.path.exists(f"{PKG}/{FIGS1}"):
        raise SystemExit("run sync_figures.py first: no packaged Supplementary Figure S1")

    # -- remove anything the package should not carry ------------------------
    n_fig = len(set(re.findall(r"\*\*Figure (\d+)\.", ms)))
    keep = {"00_CoverLetter.docx", "00_CoverLetter.pdf", "01_Manuscript.docx",
            "02_Manuscript_typeset.pdf", SF1_PDF, SF2_XLSX,
            SF3_PDF, "08_title_abstract_keywords.txt", "06_Figure_alt_text.md",
            "07_Figure_legends.docx", FIGS1}
    keep |= {f"Figure{i}.tif" for i in range(1, n_fig + 1)}
    for f in sorted(os.listdir(PKG)):
        if f not in keep:
            os.remove(os.path.join(PKG, f))
            print(f"  removed stale {f}")
    total = sum(os.path.getsize(os.path.join(PKG, f)) for f in os.listdir(PKG))
    print(f"\npackage: {len(os.listdir(PKG))} files, {total/1e6:.1f} MB")


if __name__ == "__main__":
    main()
