#!/usr/bin/env python3
"""59_direction_matched_percentile.py -- where this study's own two modules sit
inside the direction-matched benchmark distribution.

§6ad established that in liver a signature's retention tracks the direction of
its unadjusted shift, and Potential implications tells readers to compare
retention only between signatures that move the same way. The manuscript was
still quoting its own modules against the pooled distribution, which mixes both
directions. This computes and archives both versions so the text can report the
direction-matched one.
"""
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

import csv, json, statistics as st
rows=[r for r in csv.DictReader(open(f"{IR_DOCS}/SF2build/SupplementaryFile2/SupplementaryTable_S2_covariate_comparison_GSE14520.csv"))
      if r["retention_joint"] not in ("","NA")]
vals=[(float(r["retention_joint"]), float(r["unadjusted_mean_delta"])) for r in rows]
def pct(x, pool): return round(100*sum(1 for v in pool if v < x)/len(pool),1)
fall=[v for v,d in vals if d<0]; rise=[v for v,d in vals if d>0]
out={"plan":"PREREGISTRATION 6ad, reporting step",
     "note":"Percentile position of this study's own two modules within the "
            "direction-matched part of the GSE14520 benchmark distribution. §6ad "
            "established that retention tracks the direction of a signature's "
            "unadjusted shift in liver, so a percentile taken across both "
            "directions is a percentile against a mixture.",
     "n_all":len(vals),"n_falling":len(fall),"n_rising":len(rise),
     "median_falling":round(st.median(fall),4),"median_rising":round(st.median(rise),4),
     "DRAIN_retention_joint":0.3103,"REDUCTION_retention_joint":0.8905,
     "DRAIN_percentile_all":pct(0.3103,[v for v,_ in vals]),
     "DRAIN_percentile_falling":pct(0.3103,fall),
     "REDUCTION_percentile_all":pct(0.8905,[v for v,_ in vals]),
     "REDUCTION_percentile_rising":pct(0.8905,rise)}
json.dump(out, open(f"{IR_RESULTS}/DIRECTION_MATCHED_PERCENTILE_6ad.json","w"), indent=2)
print(json.dumps(out, indent=2))
