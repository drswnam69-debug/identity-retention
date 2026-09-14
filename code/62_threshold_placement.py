#!/usr/bin/env python3
"""62_threshold_placement.py -- how many signatures each cohort actually places
against the 50% threshold.

The manuscript reported "placed against the threshold" as the count of intervals
lying entirely BELOW 0.5, and drew a conclusion about cohort size from it. That
is half the question: an interval lying entirely ABOVE 0.5 places a signature
just as definitely, on the other side. Counting only one side made the 50-pair
cohort look far less informative than it is, and put a false number, 109 of 119,
into the recommendation that rests on it. This computes all three counts.

Not pre-registered; written after a pre-submission review found the conflation.
"""
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
from paths import RESULTS as IR_RESULTS

COHORTS = [("GSE14520", 213), ("GSE76427", 52), ("TCGA_LIHC", 50),
           ("TCGA_LUAD", 58), ("TCGA_KIRC", 72)]
THRESHOLD = 0.50


def main():
    out = {"plan": "PREREGISTRATION 6v, reporting step",
           "threshold": THRESHOLD,
           "note": "A signature is placed against the threshold when its 95% "
                   "interval lies entirely on one side of it. The manuscript "
                   "previously counted only the below side.",
           "cohorts": {}}
    for name, pairs in COHORTS:
        d = json.load(open(os.path.join(IR_RESULTS,
                      f"RETENTION_INTERVALS_6v_{name}.json"), encoding="utf-8"))
        sigs = [s for s in d["signatures"] if s.get("ci95_retention_joint")]
        below = above = span = 0
        for s in sigs:
            lo, hi = s["ci95_retention_joint"]
            if hi < THRESHOLD:
                below += 1
            elif lo > THRESHOLD:
                above += 1
            else:
                span += 1
        rec = {"n_pairs": pairs, "n_with_interval": len(sigs),
               "entirely_below": below, "entirely_above": above,
               "spanning": span, "placed": below + above,
               "placed_fraction": round((below + above) / len(sigs), 4),
               "median_ci_width": d["median_ci_width"]}
        out["cohorts"][name] = rec
        print(f"  {name:10s} {pairs:4d} pairs  n={len(sigs):4d}  "
              f"below {below:3d}  above {above:3d}  spanning {span:3d}  "
              f"placed {rec['placed']:3d} ({rec['placed_fraction']:.0%})")
    path = os.path.join(IR_RESULTS, "THRESHOLD_PLACEMENT_6v.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
