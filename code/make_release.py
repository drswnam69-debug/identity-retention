#!/usr/bin/env python3
"""make_release.py -- build the repository archive the article cites.

Until this submission the code reached the reader twice: once as the GitHub
repository under its Zenodo digital object identifier, and once as a
supplementary file. The two were assembled by different code and had drifted:
the pushed repository was missing six analysis scripts and six result files that
the manuscript names, `63_what_the_spread_is_not.py` among them, because only
the supplementary copy was ever re-synced. The supplementary copy is gone now,
so the repository is the only deposit and this script is the only thing that
builds it.

It is the supplementary archive with three differences, each deliberate:

  * `README.txt`, which addressed the reader of a supplementary file and still
    carried an earlier title of the article, is dropped;
  * `README_repository.md` becomes `README.md`;
  * `results/figures/` is dropped. It is a working directory of every figure
    ever drawn, including superseded ones and three files left under a British
    spelling, and the figures that belong to the article are in `figures/`.

`docs/` is carried over from the previous release, which is where it has always
lived. `figures/` is rebuilt from `sync_figures.py`'s map of article figure
number to file name, because the hand-written table in `figures/README.md` had
drifted from that map and named the wrong file for five of the thirteen figures.
"""
import hashlib
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "rsi", "code"))
from paths import DOCS as IR_DOCS

HOME = IR_DOCS
SRC = f"{HOME}/SF2build/SupplementaryFile2"
OUT = f"{HOME}/release/identity-retention"
KEEP_FROM_PREVIOUS = ("docs",)
DROP_FILES = {"README.txt", "SHA256SUMS.txt"}
DROP_DIRS = {os.path.join("results", "figures")}


def write_figures() -> None:
    """Rebuild figures/ and its index from the one map that the build uses.

    The index used to be written by hand. It named the wrong file for Figures 5,
    6, 8, 11 and 12, which is what a second copy of a mapping does."""
    sys.path.insert(0, HOME)
    from sync_figures import MAP, SEARCH
    d = os.path.join(OUT, "figures")
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    rows, missing = [], []
    # the graphical abstract is an article item with no figure number, so it is
    # not in sync_figures' map and has to be named here
    ga = "GraphicalAbstract_identity_retention"
    for n, stem in list(MAP.items()) + [("Graphical abstract", ga)]:
        for src in SEARCH:
            p = os.path.join(src, stem + ".png")
            if os.path.exists(p):
                shutil.copy2(p, os.path.join(d, stem + ".png"))
                rows.append((n if n == "Graphical abstract" else f"Figure {n}",
                             stem + ".png"))
                break
        else:
            missing.append(stem)
    if missing:
        raise SystemExit("no PNG for: " + ", ".join(missing))
    with open(os.path.join(d, "README.md"), "w", encoding="utf-8") as fh:
        fh.write("# Figures in this directory\n\n"
                 "The file names are the names the plotting scripts write. The "
                 "numbers are the\nnumbers the article uses. This table is "
                 "generated from the map in\n`code/sync_figures.py`, which is "
                 "what puts each figure into the built document,\nso the two "
                 "cannot disagree.\n\n| Article | File |\n|---|---|\n")
        for a, f in rows:
            fh.write(f"| {a} | `{f}` |\n")
    print(f"  figures/: {len(rows)} figures and a generated index")


def main() -> int:
    carried = {}
    for d in KEEP_FROM_PREVIOUS:
        p = os.path.join(OUT, d)
        if not os.path.isdir(p):
            raise SystemExit(f"the previous release has no {d}/ to carry over")
        tmp = f"{HOME}/build/_carry_{d}"
        shutil.rmtree(tmp, ignore_errors=True)
        shutil.copytree(p, tmp)
        carried[d] = tmp

    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)

    n = 0
    for base, dirs, files in os.walk(SRC):
        rel = os.path.relpath(base, SRC)
        rel = "" if rel == "." else rel
        dirs[:] = [d for d in sorted(dirs)
                   if d != "__pycache__" and os.path.join(rel, d) not in DROP_DIRS]
        for f in sorted(files):
            if f.endswith(".pyc") or (rel == "" and f in DROP_FILES):
                continue
            name = "README.md" if (rel == "" and f == "README_repository.md") else f
            dst = os.path.join(OUT, rel, name)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(os.path.join(base, f), dst)
            n += 1

    for d, tmp in carried.items():
        shutil.copytree(tmp, os.path.join(OUT, d))
        shutil.rmtree(tmp, ignore_errors=True)
        n += sum(len(fs) for _b, _d, fs in os.walk(os.path.join(OUT, d)))

    write_figures()

    lines = []
    for base, dirs, files in os.walk(OUT):
        dirs[:] = [x for x in sorted(dirs) if x != "__pycache__"]
        for f in sorted(files):
            if f.endswith(".pyc"):
                continue
            full = os.path.join(base, f)
            h = hashlib.sha256(open(full, "rb").read()).hexdigest()
            lines.append(f"{h}  {os.path.relpath(full, OUT)}")
    open(os.path.join(OUT, "SHA256SUMS.txt"), "w").write(
        "\n".join(sorted(lines, key=lambda l: l.split("  ", 1)[1])) + "\n")
    print(f"release/identity-retention: {len(lines)} files, checksummed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
