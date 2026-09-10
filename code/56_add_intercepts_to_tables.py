#!/usr/bin/env python3
"""56_add_intercepts_to_tables.py -- carry the signed intercept into S2, S5, S6.

The retention fraction is |intercept| over |unadjusted shift|, so it is
positive whether adjustment shrinks a signature's shift or reverses its sign.
The manuscript says so, and tells the reader that Supplementary Tables S2, S5
and S6 carry the signed intercept beside every fraction so the two readings can
be told apart. They did not: no supplementary table carried an intercept
column at all, so the one sentence in the paper that resolves the ambiguity
pointed at data that was not there. This adds the column, from the archived
benchmarks, and a flag for the rows where the sign reverses.

Run after the benchmarks; it reads results/ and rewrites the CSVs in place.
"""
from __future__ import annotations

import csv
import json
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_cands = [_here, os.path.join(_here, "code")]
if os.environ.get("IR_ROOT"):
    _cands.insert(0, os.path.join(os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if os.path.exists(os.path.join(_c, "paths.py")):
        if _c not in sys.path:
            sys.path.insert(0, _c)
        break
from paths import RESULTS, DOCS  # noqa: E402

SF2 = os.environ.get("IR_SF2", f"{DOCS}/SF2build/SupplementaryFile2")
JOBS = [
    ("SupplementaryTable_S2_covariate_comparison_GSE14520.csv", "GSE14520"),
    ("SupplementaryTable_S2_covariate_comparison_GSE76427.csv", "GSE76427"),
    ("SupplementaryTable_S5_TCGA_LIHC_retention.csv", "TCGA_LIHC"),
    ("SupplementaryTable_S6_lung_retention.csv", "TCGA_LUAD"),
    ("SupplementaryTable_S6_kidney_retention.csv", "TCGA_KIRC"),
]
NEW = ["joint_intercept", "sign_reverses", "ratio_unstable"]


def index(cohort: str):
    p = f"{RESULTS}/SIGNATURE_BENCHMARK_{cohort}.json"
    b = {s["name"]: s for s in json.load(open(p, encoding="utf-8"))["signatures"]}
    unstable: set[str] = set()
    ip = f"{RESULTS}/RETENTION_INTERVALS_6v_{cohort}.json"
    if os.path.exists(ip):
        d = json.load(open(ip, encoding="utf-8"))
        unstable = set(d.get("unstable_names", []))
        for s in d.get("signatures", []):
            if s.get("ratio_unstable"):
                unstable.add(s["name"])
    return b, unstable


def main() -> int:
    for fname, cohort in JOBS:
        path = os.path.join(SF2, fname)
        if not os.path.exists(path):
            print(f"  skip (absent): {fname}")
            continue
        bench, unstable = index(cohort)
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        head, body = rows[0], rows[1:]
        key = head.index("signature")
        for c in NEW:
            if c in head:
                col = head.index(c)
                head.pop(col)
                for r in body:
                    r.pop(col)
        head = head + NEW
        n_rev = 0
        for r in body:
            s = bench.get(r[key], {})
            ic = s.get("joint_intercept")
            un = s.get("unadjusted_mean_delta")
            beta = ic["beta"] if isinstance(ic, dict) else None
            rev = ""
            if beta is not None and un not in (None, 0):
                rev = "1" if beta * un < 0 else "0"
                n_rev += rev == "1"
            r += ["" if beta is None else f"{beta:+.4f}", rev,
                  "1" if r[key] in unstable else "0"]
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(head)
            w.writerows(body)
        print(f"  {fname}: {len(body)} rows, {n_rev} sign reversals, "
              f"{sum(1 for r in body if r[-1] == '1')} unstable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
