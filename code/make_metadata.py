#!/usr/bin/env python3
"""make_metadata.py -- rebuild 05_title_abstract_keywords.txt from the manuscript.

This file is what the submission system's title and abstract boxes are filled
from. It was written by hand once and then went stale when the abstract was
revised, which no check would have caught, so it is generated from the
manuscript instead and verified against it by 50_consistency_check.py.
"""
import re


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
MS = f"{IR_DOCS}/manuscript_v4.md"
OUT = f"{IR_DOCS}/submission_GigaScience/05_title_abstract_keywords.txt"

ms = open(MS, encoding="utf-8").read()
title = re.search(r"^# (.+)$", ms, re.M).group(1).strip()
run = re.search(r"\*\*Running title\.?\*\*\s*(.+)", ms)
running = run.group(1).strip() if run else "Identity retention in tumor expression signatures"
ab = re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", ms, re.S).group(1).strip()
kw = re.search(r"\*\*Keywords\.\*\*\s*(.+)", ms).group(1).strip()

plain = re.sub(r"\*\*(.+?)\.\*\*", r"\1.", ab)
plain = plain.replace("**", "").replace("*", "")
n = len(ab.split())

open(OUT, "w", encoding="utf-8").write(f"""ARTICLE TYPE
------------
Research

TITLE
-----
{title}

RUNNING TITLE
-------------
{running}

ABSTRACT (structured, {n} words)
{'-' * (len(f'ABSTRACT (structured, {n} words)'))}

{plain}

KEYWORDS
--------
{kw}
""")
print(f"wrote {OUT}  ({n} words)")
