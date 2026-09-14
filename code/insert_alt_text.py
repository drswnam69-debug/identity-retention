#!/usr/bin/env python3
"""insert_alt_text.py -- put each figure's alt text into the manuscript.

GigaScience requires that each alt text description sit "in your main manuscript
file, directly under the figure legend for the relevant figure, preceded by
'Alt text:'". The descriptions were maintained only as a standalone companion
file, which satisfies the submission form and not the manuscript rule, so this
copies them under their legends and rewrites the block on every run rather than
leaving two copies to drift apart.
"""
import re
import sys

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
from paths import DOCS as IR_DOCS

MS = f"{IR_DOCS}/manuscript_v4.md"
ALT = f"{IR_DOCS}/Figure_alt_text.md"


def main():
    ms = open(MS, encoding="utf-8").read()
    alt = open(ALT, encoding="utf-8").read()
    alts = {m.group(1): m.group(2).strip() for m in
            re.finditer(r"\*\*(?:Supplementary )?Figure (\d+|S\d+)\.\*\*(.*?)(?=\n\n|\Z)",
                        alt, re.S)}
    if not alts:
        raise SystemExit("no alt text entries found")

    # drop any previously inserted block so the manuscript never carries two
    ms = re.sub(r"\n\n\*\*Alt text:\*\*[^\n]*", "", ms)

    added = 0
    for n in sorted(alts, key=lambda x: (x.startswith("S"), int(x.lstrip("S")))):
        # the numbered legends read "**Figure 3. Title.** text"; the supplementary
        # entry in Additional files reads "**Supplementary Figure S1** (file name ...)"
        pat = re.compile(r"(\*\*(?:Supplementary )?Figure " + re.escape(n) +
                         r"[.*].*?)(?=\n\n)", re.S)
        m = pat.search(ms)
        if not m:
            print(f"  no legend for Figure {n}; skipped")
            continue
        block = "\n\n**Alt text:** " + " ".join(alts[n].split())
        ms = ms[:m.end()] + block + ms[m.end():]
        added += 1
    open(MS, "w", encoding="utf-8").write(ms)
    print(f"inserted {added} alt text blocks under their legends")
    return 0


if __name__ == "__main__":
    sys.exit(main())
