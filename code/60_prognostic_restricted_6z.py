#!/usr/bin/env python3
"""60_prognostic_restricted_6z.py -- the §6z correlation among signatures that
are actually prognostic.

The registered test correlates retention with the surviving fraction of the log
hazard ratio across all 103 evaluable signatures. For the 21 that are not
prognostic before adjustment, that surviving fraction is a ratio of two
quantities indistinguishable from zero, which is the instability the §6v guard
exists to catch. This recomputes the same correlation on the prognostic subsets
so the manuscript can say how much of the registered estimate rests on them.
Nothing here replaces the registered verdict; it qualifies a sentence about it.
"""
import csv
import json
import math
import os
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
from paths import RESULTS as IR_RESULTS, DOCS as IR_DOCS

SRC = f"{IR_DOCS}/SF2build/SupplementaryFile2/SupplementaryTable_S7_prognostic.csv"


def spearman(a, b):
    n = len(a)

    def rank(x):
        idx = sorted(range(n), key=lambda i: x[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and x[idx[j + 1]] == x[idx[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[idx[k]] = avg
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((x - mb) ** 2 for x in rb))
    r = num / den
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(n - 3)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return round(r, 4), p, n


def main():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    out = {"plan": "PREREGISTRATION 6z, sensitivity added before submission",
           "note": "The registered estimate is the all-signature one. The two "
                   "restricted estimates show how much of it rests on signatures "
                   "whose unadjusted hazard ratio is indistinguishable from zero, "
                   "for which the surviving fraction is a ratio of two near-zero "
                   "quantities.",
           "subsets": {}}
    for label, keep in (("all", lambda r: True),
                        ("prognostic_p05", lambda r: float(r["p_unadj"]) < 0.05),
                        ("prognostic_bh05", lambda r: float(r["q_unadj"]) < 0.05)):
        sub = [r for r in rows if keep(r)]
        rho, p, n = spearman([float(r["retention_joint"]) for r in sub],
                             [float(r["surviving_fraction"]) for r in sub])
        out["subsets"][label] = {"n": n, "spearman_rho": rho, "spearman_p": p}
        print(f"  {label:16s} n={n:3d}  rho={rho:+.4f}  P={p:.4f}")
    out["registered_threshold"] = 0.30
    path = os.path.join(IR_RESULTS, "PROGNOSTIC_RESTRICTED_6z.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
