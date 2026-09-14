#!/usr/bin/env python3
"""61_signature_independence.py -- how independent the benchmark signatures are.

Every across-signature correlation in this paper treats the panel as that many
independent observations. The C2:CGP collection is not built that way: sets come
in families from the same source study, and members of a family share genes. This
counts the source studies behind each cohort's panel and measures pairwise gene
overlap, so the manuscript can state the dependence rather than assume it away.
Written after a pre-submission review; not pre-registered.
"""
import collections
import csv
import json
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
from paths import RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, DOCS as IR_DOCS

S8 = f"{IR_DOCS}/SF2build/SupplementaryFile2/SupplementaryTable_S8_enumeration.csv"


def load_gmt():
    """Every gene set in the deposited C2:CGP collection, by name."""
    sets = {}
    for base, _dirs, files in os.walk(IR_GENESETS):
        for f in files:
            if not f.endswith(".gmt"):
                continue
            for line in open(os.path.join(base, f), encoding="utf-8"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 2:
                    sets[parts[0]] = set(parts[2:])
    return sets


def main():
    rows = [r for r in csv.DictReader(open(S8, encoding="utf-8"))
            if r["status"] == "included"]
    gmt = load_gmt()
    out = {"note": "Source-study composition and pairwise gene overlap of each "
                   "cohort's comparator panel. Across-signature statistics in the "
                   "manuscript treat these as independent; they are not.",
           "cohorts": {}}
    for cohort in sorted({r["cohort"] for r in rows}):
        sub = [r for r in rows if r["cohort"] == cohort]
        names = [r["set_name"] for r in sub]
        fam = collections.Counter(n.split("_")[0] for n in names)
        rec = {"n_sets": len(names), "n_source_prefixes": len(fam),
               "largest_family": fam.most_common(1)[0][1],
               "largest_family_name": fam.most_common(1)[0][0],
               "top_families": fam.most_common(5)}
        have = [n for n in names if n in gmt]
        if len(have) > 1:
            j_over_10 = j_over_25 = 0
            mx = 0.0
            for i in range(len(have)):
                for k in range(i + 1, len(have)):
                    a, b = gmt[have[i]], gmt[have[k]]
                    u = len(a | b)
                    j = len(a & b) / u if u else 0.0
                    mx = max(mx, j)
                    if j > 0.10:
                        j_over_10 += 1
                    if j > 0.25:
                        j_over_25 += 1
            rec.update(n_sets_with_membership=len(have),
                       n_pairs_jaccard_over_0_10=j_over_10,
                       n_pairs_jaccard_over_0_25=j_over_25,
                       max_jaccard=round(mx, 4))
        out["cohorts"][cohort] = rec
        print(f"  {cohort:10s} {rec['n_sets']:3d} sets from {rec['n_source_prefixes']:2d} "
              f"source prefixes, largest {rec['largest_family']} "
              f"({rec['largest_family_name']})")
    path = os.path.join(IR_RESULTS, "SIGNATURE_INDEPENDENCE.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
