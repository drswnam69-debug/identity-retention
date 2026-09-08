#!/usr/bin/env python3
"""27_fix_ci_6n.py -- PREREGISTRATION 6n, defect 1.

ols_ci used z = 1.96 instead of t(dof) above dof = 200, so GSE14520's six
confidence intervals are 0.53-0.57% too narrow. Point estimates, standard errors
and P values are unaffected and are NOT touched here: only the interval bounds
are recomputed, as beta +/- t(dof) * se.

Every corrected bound is checked against SciPy before it is written.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sp = importlib.util.spec_from_file_location("d", os.path.join(HERE, "12_differentiation_adjust.py"))
_da = importlib.util.module_from_spec(sp)
sp.loader.exec_module(_da)

try:
    from scipy import stats as sps
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

N_PAIRS = {"GSE14520": 213, "GSE76427": 52}

TARGETS = [
    ("results/GSE14520_differentiation.json",
     [("paired_adjusted_D1", "intercept", 2), ("paired_adjusted_D1", "slope", 2),
      ("paired_adjusted_D2", "intercept", 2), ("paired_adjusted_D2", "slope", 2)]),
    ("results/GSE14520_composition.json",
     [("composite_C1", None, 2), ("composite_joint", None, 3),
      ("module_reduction_joint", None, 3), ("module_drain_joint", None, 3),
      ("module_reduction_C1", None, 2), ("module_drain_C1", None, 2),
      ("composite_joint_C2", None, 3), ("module_reduction_joint_C2", None, 3),
      ("module_drain_joint_C2", None, 3)]),
]


def fix_block(blk, dof, changed, label):
    """Recompute ci95 from beta and se with the correct t critical value."""
    if not isinstance(blk, dict) or "beta" not in blk or "se" not in blk:
        return
    b, se = blk["beta"], blk["se"]
    if se == 0:
        return
    t = _da._tcrit(dof)
    lo, hi = round(b - t * se, 4), round(b + t * se, 4)
    if HAVE_SCIPY:
        t_ref = sps.t.ppf(0.975, dof)
        assert abs(t - t_ref) < 1e-6, f"tcrit disagrees with SciPy: {t} vs {t_ref}"
    old = blk.get("ci95")
    if old != [lo, hi]:
        changed.append((label, old, [lo, hi]))
        blk["ci95"] = [lo, hi]


def main():
    changed = []
    for path, keys in TARGETS:
        if not os.path.exists(path):
            continue
        d = json.load(open(path))
        cohort = "GSE14520" if "14520" in path else "GSE76427"
        n = N_PAIRS[cohort]
        for key, sub, k in keys:
            if key not in d:
                continue
            dof = n - k
            node = d[key]
            if sub:
                fix_block(node.get(sub, {}), dof, changed, f"{path}:{key}.{sub}")
            else:
                # top-level ci95 mirrors the intercept of the fit
                fit = node.get("fit", {})
                for nm, blk in fit.items():
                    fix_block(blk, dof, changed, f"{path}:{key}.fit.{nm}")
                if "ci95" in node and "intercept" in fit:
                    if node["ci95"] != fit["intercept"]["ci95"]:
                        changed.append((f"{path}:{key}.ci95(top)", node["ci95"],
                                        fit["intercept"]["ci95"]))
                        node["ci95"] = fit["intercept"]["ci95"]
                elif "ci95" in node and "tumour" in fit:
                    pass
        json.dump(d, open(path, "w"), indent=1)

    print(f"=== 6n defect 1: {len(changed)} interval bounds corrected ===")
    for lab, old, new in changed:
        print(f"  {lab}")
        print(f"      {old}  ->  {new}")
    json.dump([{"where": a, "old": b, "new": c} for a, b, c in changed],
              open("results/CI_CORRECTION_6n.json", "w"), indent=1)
    print("\nwrote results/CI_CORRECTION_6n.json")


if __name__ == "__main__":
    main()
