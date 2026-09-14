#!/usr/bin/env python3
"""58_covariate_movement.py -- how far each cohort's own covariates move, beside
the retention median that cohort produced.

Written after a pre-submission review pointed out that the three tissue medians
order themselves by the summed movement of that cohort's identity and
composition covariates rather than by tissue, which §6ae already showed governs
how much a covariate removes within one cohort. The manuscript reports these
sums, so they are computed and archived here rather than done by hand in the
text. Nothing here is a new analysis: every input is an already-archived paired
shift or median.
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

R = IR_RESULTS


def j(path):
    with open(os.path.join(R, path), encoding="utf-8") as f:
        return json.load(f)


# identity shift, composition shift, retention median: each from the file of record
COHORTS = [
    ("TCGA_LIHC", "TCGA_LIHC/PREMISE_6w.json", "D1_paired_delta", "C1_paired_delta", 0.8752),
    ("GSE76427", None, None, None, 0.6080),
    ("GSE14520", None, None, None, 0.7579),
    ("TCGA_LUAD", "TCGA_LUAD/PREMISE_6y.json", "identity_paired_delta", "C1_paired_delta", 0.5761),
    ("TCGA_KIRC", "TCGA_KIRC/PREMISE_6x.json", "identity_paired_delta", "C1_paired_delta", 0.4665),
]


def main():
    out = {"plan": "PREREGISTRATION 6ae, reporting step",
           "note": "Summed absolute paired shift of each cohort's identity and "
                   "composition covariates, beside that cohort's median joint "
                   "retention. Inputs are archived elsewhere; nothing is refitted.",
           "cohorts": {}}
    for name, premise, kd, kc, median in COHORTS:
        if name == "GSE14520":
            # identity shift from the simulation record, C1 from the composition run
            d = j("ESTIMATOR_SIMULATION_6v.json")["mean_dD1"]
            c1 = -0.2164
        elif name == "GSE76427":
            # this cohort archives D1 as the two tissue means rather than a delta
            pc = j("GSE76427_differentiation.json")["positive_control"]
            d = pc["d1_tumor"] - pc["d1_adjacent"]
            c1 = j("GSE76427_composition.json")["positive_control"]["delta"]
        else:
            p = j(premise)
            d, c1 = p[kd], p[kc]
        out["cohorts"][name] = {
            "identity_paired_delta": round(d, 4),
            "C1_paired_delta": round(c1, 4),
            "summed_absolute_movement": round(abs(d) + abs(c1), 4),
            "median_retention_joint": median,
        }
    order = sorted(out["cohorts"], key=lambda k: out["cohorts"][k]["summed_absolute_movement"])
    mv = [out["cohorts"][k]["summed_absolute_movement"] for k in order]
    med = [out["cohorts"][k]["median_retention_joint"] for k in order]
    out["ordered_by_movement"] = order
    out["movement_ascending"] = mv
    out["median_retention_in_that_order"] = med
    out["median_is_monotone_decreasing"] = all(med[i] > med[i + 1] for i in range(len(med) - 1))
    n = len(med)
    rm = sorted(range(n), key=lambda i: med[i])
    rank = [0] * n
    for r, i in enumerate(rm):
        rank[i] = r + 1
    dsq = sum((i + 1 - rank[i]) ** 2 for i in range(n))
    out["spearman_movement_vs_median"] = round(1 - 6 * dsq / (n * (n * n - 1)), 4)
    # GSE76427's median comes from the archived gene-by-gene assembly that the
    # benchmark code does not reproduce, so it is the one cohort whose position
    # here should be read with that caveat.
    out["caveat_GSE76427"] = ("median from the archived assembly this study reports "
                              "as irreproducible; excluded from per-signature claims "
                              "elsewhere")
    path = os.path.join(R, "COVARIATE_MOVEMENT_6ae.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    for k in order:
        v = out["cohorts"][k]
        print(f"  {k:11s} movement {v['summed_absolute_movement']:.3f}   "
              f"median retention {v['median_retention_joint']:.4f}")
    print("  medians decrease monotonically as movement grows:",
          out["median_is_monotone_decreasing"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
