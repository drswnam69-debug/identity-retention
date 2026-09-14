"""paths.py -- where this archive's files are, resolved from the archive itself.

Every script here used to name absolute paths on the machine the analysis was
run on. That is invisible to the author and fatal to everyone else: a reader
who downloads this repository got FileNotFoundError from a third of the code,
including every figure script, while the file it wanted sat in the download.

Resolution walks up from this file, so it is correct wherever the archive is
unpacked and whichever directory a script is run from. IR_ROOT overrides it.

Author-side documents (the manuscript, the cover letter, the submission
package) are deliberately NOT part of this archive before publication. The
scripts that read them take a path argument; IR_DOCS names their directory.
"""
from __future__ import annotations

import os

__all__ = ["ROOT", "CODE", "RESULTS", "GENESETS", "DATA", "FIGURES", "TOOL",
           "DOCS", "doc", "res", "describe"]


def _find_root() -> str:
    env = os.environ.get("IR_ROOT")
    if env:
        return os.path.abspath(env)
    here = os.path.dirname(os.path.abspath(__file__))
    cur = here
    for _ in range(6):
        if all(os.path.isdir(os.path.join(cur, d)) for d in ("results", "genesets")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    # last resort: the directory above code/, whatever it holds
    return os.path.dirname(here)


ROOT = _find_root()
CODE = os.path.join(ROOT, "code")
RESULTS = os.path.join(ROOT, "results")
GENESETS = os.path.join(ROOT, "genesets")
DATA = os.path.join(ROOT, "data")
FIGURES = os.path.join(ROOT, "figures")
TOOL = os.path.join(ROOT, "tool")

# Documents that live outside the archive until the article is published.
DOCS = os.path.abspath(os.environ.get("IR_DOCS", os.path.dirname(ROOT)))


def res(*parts: str) -> str:
    """A path inside results/."""
    return os.path.join(RESULTS, *parts)


def doc(name: str) -> str:
    """A path to an author-side document, under IR_DOCS."""
    return os.path.join(DOCS, name)


def describe() -> str:
    return (f"archive root : {ROOT}\n"
            f"results      : {RESULTS}\n"
            f"genesets     : {GENESETS}\n"
            f"data         : {DATA}\n"
            f"documents    : {DOCS}  (set IR_DOCS to move)")


if __name__ == "__main__":
    print(describe())
