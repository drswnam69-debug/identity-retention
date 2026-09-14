#!/usr/bin/env python3
"""Renumber the reference list into first-appearance order and rewrite every
citation to match. Fails loudly rather than guessing."""
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
P = f"{IR_DOCS}/manuscript_v4.md"
s = open(P, encoding="utf-8").read()
head = "\n## References\n"
tail_head = "\n## Tables\n"
i, j = s.index(head), s.index(tail_head)
refs, tail = s[i:j], s[j:]
body = s[:i]  # citations before the list
tail_cites = tail    # tables, figure legends and additional files cite too

entries = {int(m.group(1)): m.group(2).strip()
           for m in re.finditer(r"^(\d+)\. (.*)$", refs, re.M)}
assert len(entries) == len(re.findall(r"^\d+\. ", refs, re.M)), "duplicate reference numbers"

def expand(tok):
    out = []
    for part in tok.split(","):
        part = part.strip().replace("–", "-")
        if "-" in part:
            a, b = part.split("-"); out.extend(range(int(a), int(b) + 1))
        elif part.isdigit(): out.append(int(part))
    return out

order = []
for m in re.finditer(r"\[([0-9][0-9,\s–-]*)\]", body + tail_cites):
    for n in expand(m.group(1)):
        if n not in order: order.append(n)
missing = sorted(set(entries) - set(order))
extra = sorted(set(order) - set(entries))
assert not extra, f"cited but not in the list: {extra}"
assert not missing, f"in the list but never cited: {missing}"

mapping = {old: new for new, old in enumerate(order, start=1)}

def rewrite(m):
    nums = expand(m.group(1))
    new = sorted(mapping[n] for n in nums)
    # collapse runs of three or more
    out, i2 = [], 0
    while i2 < len(new):
        j = i2
        while j + 1 < len(new) and new[j + 1] == new[j] + 1: j += 1
        out.append(f"{new[i2]}-{new[j]}" if j - i2 >= 2 else ", ".join(str(x) for x in new[i2:j + 1]))
        i2 = j + 1
    return "[" + ", ".join(out) + "]"

body = re.sub(r"\[([0-9][0-9,\s–-]*)\]", rewrite, body)
tail = re.sub(r"\[([0-9][0-9,\s–-]*)\]", rewrite, tail)

# The supplementary note cites the manuscript's list and lives in its own file,
# so it has to be renumbered in the same pass or its numbers silently drift.
NOTE = f"{IR_DOCS}/SupplementaryNote_v1.md"
if _os.path.exists(NOTE):
    nt = open(NOTE, encoding="utf-8").read()
    nt = re.sub(r"\[([0-9][0-9,\s–-]*)\]", rewrite, nt)
    open(NOTE, "w", encoding="utf-8").write(nt)
    print("renumbered the supplementary note's citations too")
lines = [f"{mapping[old]}. {entries[old]}" for old in sorted(entries, key=lambda o: mapping[o])]
open(P, "w", encoding="utf-8").write(body + head + "\n" + "\n\n".join(lines) + "\n" + tail)
print(f"renumbered {len(entries)} references into first-appearance order")
print("map:", {k: v for k, v in sorted(mapping.items())})
