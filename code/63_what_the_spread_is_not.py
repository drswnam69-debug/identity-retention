#!/usr/bin/env python3
"""63_what_the_spread_is_not.py -- two objections to the per-signature spread.

The paper's surviving claim is that retention differs between signatures with
the covariate held fixed. A referee has two cheap ways to dismiss it: the spread
is sampling noise, or it is the unadjusted effect size wearing a new name.

The first version of this script decomposed the variance on the RATIO scale,
using each bootstrap interval's width over 3.92 as a within-signature standard
deviation. A pre-submission referee showed that this is not defensible: the
quantity is a ratio whose denominator can approach zero, so its bootstrap
distribution is heavy-tailed and its true variance is dominated by a few sets.
Re-running the §6v bootstrap and taking the actual standard deviation of the
resampled ratios gives a mean within-signature variance of about 1.31 against a
between-signature variance of 0.14, i.e. a NEGATIVE variance component. The
interval-width proxy only appeared to work because the 2.5/97.5 percentile cut
truncates exactly those tails.

On the log scale the resampled ratio is close to symmetric and the decomposition
is stable: it gives the same answer whether the within-signature term comes from
the interval width or from the true bootstrap standard deviation. That is the
version reported. A cluster bootstrap over the 33 source studies is included
because the sets are not independent, and the source-study effect is tested on
both scales because it is detectable on one and not the other.

Not pre-registered. Written after a pre-submission review asked what a reader is
meant to do with the number, and rewritten after a second review showed the
first calculation was wrong.
"""
import collections
import importlib.util
import json
import math
import os
import random
import statistics as st
import sys

import numpy as np
import pandas as pd

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
from paths import RESULTS as IR_RESULTS, DATA as IR_DATA, CODE as IR_CODE, \
    GENESETS as IR_GENESETS

B, SEED = 2000, 20260907          # the §6v settings, unchanged
BOOT, BSEED = 2000, 20260913      # for the interval on the share


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(IR_CODE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def batched_intercept(y, X, idx):
    Xb, yb = X[idx], y[idx]
    XtX = np.einsum("bij,bik->bjk", Xb, Xb)
    Xty = np.einsum("bij,bi->bj", Xb, yb)
    return np.linalg.solve(XtX, Xty[..., None])[..., 0][:, 0]


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
            for k in range(i, j + 1):
                r[idx[k]] = (i + j) / 2 + 1
            i = j + 1
        return r

    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((x - mb) ** 2 for x in rb))
    r = num / den
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(n - 3)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    lo = math.tanh(math.atanh(r) - 1.96 / math.sqrt(n - 3))
    hi = math.tanh(math.atanh(r) + 1.96 / math.sqrt(n - 3))
    return round(r, 4), p, [round(lo, 4), round(hi, 4)]


def main():
    _da = _load("diffadj", "12_differentiation_adjust.py")
    _ca = _load("compadj", "18_composition_adjust.py")
    _bm = _load("bench", "30_signature_benchmark.py")
    D1, C1, zmean = _da.D1, _ca.C1, _bm.zmean

    base = os.path.join(IR_RESULTS, "GSE14520")
    rsi = pd.read_csv(os.path.join(base, "rsi.tsv"), sep="\t", index_col=0)
    ph = pd.read_csv(os.path.join(base, "phenotype.tsv"), sep="\t",
                     index_col=0).loc[rsi.index]
    expr = pd.read_csv(os.path.join(IR_DATA, "GSE14520_symbols.tsv.gz"),
                       sep="\t", index_col=0)[rsi.index]
    is_t = (~ph["tissue"].astype(str).str.lower().str.contains("adjacent")).to_numpy()
    pid = (ph["patient_id"].astype(str) if "patient_id" in ph.columns
           else pd.Series(rsi.index, index=rsi.index).astype(str))
    fr = pd.DataFrame({"pid": pid.values, "t": is_t}, index=rsi.index)

    def paired(s):
        f = fr.assign(v=s.to_numpy())
        tt = f[f.t].set_index("pid")["v"]; tt = tt[~tt.index.duplicated()]
        nn = f[~f.t].set_index("pid")["v"]; nn = nn[~nn.index.duplicated()]
        k = sorted(set(tt.index) & set(nn.index))
        return (tt.loc[k] - nn.loc[k]).to_numpy(float)

    d_d1, d_c1 = paired(zmean(expr, D1)[0]), paired(zmean(expr, C1)[0])
    n = len(d_d1)
    X_J = np.column_stack([np.ones(n), d_d1, d_c1])
    IDX = np.random.default_rng(SEED).integers(0, n, size=(B, n))

    gs = json.load(open(os.path.join(IR_GENESETS, "eligible.json"), encoding="utf-8"))
    V = json.load(open(os.path.join(IR_RESULTS,
                  "RETENTION_INTERVALS_6v_GSE14520.json"), encoding="utf-8"))
    ok = [s for s in V["signatures"] if s.get("status") == "ok"]
    stable = [s for s in ok if s.get("ci95_retention_joint")]

    rows, reproduced = [], 0
    for s in stable:
        genes = gs.get(s["name"])
        if genes is None:
            continue
        dv = paired(zmean(expr, genes)[0])
        den = dv[IDX].mean(axis=1)
        ratio = np.abs(batched_intercept(dv, X_J, IDX)) / np.abs(den)
        lo, hi = (round(float(x), 4) for x in np.percentile(ratio, [2.5, 97.5]))
        if [lo, hi] == s["ci95_retention_joint"]:
            reproduced += 1
        rows.append({"name": s["name"], "pt": s["retention_joint"],
                     "ci": s["ci95_retention_joint"],
                     "sd_ratio": float(ratio.std(ddof=1)),
                     "sd_log": float(np.log(ratio).std(ddof=1)),
                     "delta": s["unadjusted_mean_delta"]})
    print(f"  reproduced {reproduced} of {len(rows)} archived intervals")

    lp = [math.log(r["pt"]) for r in rows]
    lw_sd = [r["sd_log"] ** 2 for r in rows]
    lw_ci = [((math.log(r["ci"][1]) - math.log(r["ci"][0])) / 3.9199) ** 2 for r in rows]
    tot = st.variance(lp)
    out = {"plan": "PREREGISTRATION 6v, reporting step; NOT pre-registered",
           "cohort": "GSE14520", "scale": "natural log of the retention fraction",
           "why_log": "on the ratio scale the within-signature variance is "
                      "dominated by a few sets whose denominator approaches zero "
                      "(mean 1.31 against a total of 0.14, giving a negative "
                      "variance component); on the log scale the decomposition is "
                      "stable under both estimators of the within term",
           "n": len(rows),
           "variance_of_log_point_estimates": round(tot, 4),
           "mean_within_from_bootstrap_sd": round(st.mean(lw_sd), 4),
           "mean_within_from_interval_width": round(st.mean(lw_ci), 4),
           "between_signature_variance": round(tot - st.mean(lw_sd), 4),
           "between_signature_share": round((tot - st.mean(lw_sd)) / tot, 4),
           "between_signature_share_interval_estimator":
               round((tot - st.mean(lw_ci)) / tot, 4),
           "between_signature_sd_log": round(math.sqrt(tot - st.mean(lw_sd)), 4),
           "between_signature_fold":
               round(math.exp(math.sqrt(tot - st.mean(lw_sd))), 2),
           "ratio_scale_mean_within_variance": round(st.mean(
               [r["sd_ratio"] ** 2 for r in rows]), 3)}
    print(f"  log-scale between share {out['between_signature_share']:.3f} "
          f"(interval estimator {out['between_signature_share_interval_estimator']:.3f}), "
          f"fold {out['between_signature_fold']}")

    def share(idx):
        a = [lp[i] for i in idx]
        b = [lw_sd[i] for i in idx]
        t = st.variance(a)
        return (t - st.mean(b)) / t

    rng = random.Random(BSEED)
    bs = sorted(share([rng.randrange(len(rows)) for _ in range(len(rows))])
                for _ in range(BOOT))
    out["share_ci95_signature_bootstrap"] = [round(bs[int(.025 * BOOT)], 3),
                                             round(bs[int(.975 * BOOT)], 3)]
    fam = collections.defaultdict(list)
    for i, r in enumerate(rows):
        fam[r["name"].split("_")[0]].append(i)
    keys = list(fam)
    cb = []
    for _ in range(BOOT):
        pick = []
        for _ in range(len(keys)):
            pick += fam[keys[rng.randrange(len(keys))]]
        cb.append(share(pick))
    cb.sort()
    out["n_source_studies"] = len(keys)
    out["share_ci95_cluster_bootstrap_by_source_study"] = [
        round(cb[int(.025 * BOOT)], 3), round(cb[int(.975 * BOOT)], 3)]
    print(f"  95% CI {out['share_ci95_signature_bootstrap']} by signature, "
          f"{out['share_ci95_cluster_bootstrap_by_source_study']} clustered by study")

    # is the between-signature part really between source studies?
    def anova(vals):
        g = collections.defaultdict(list)
        for r, v in zip(rows, vals):
            g[r["name"].split("_")[0]].append(v)
        allv = [x for v in g.values() for x in v]
        k, m = len(g), len(allv)
        gm = st.mean(allv)
        ssb = sum(len(v) * (st.mean(v) - gm) ** 2 for v in g.values())
        sst = sum((x - gm) ** 2 for x in allv)
        F = (ssb / (k - 1)) / ((sst - ssb) / (m - k))
        return {"k_groups": k, "n": m, "eta_squared": round(ssb / sst, 4),
                "F": round(F, 3), "df": [k - 1, m - k],
                "eta_squared_expected_under_null": round((k - 1) / (m - 1), 4)}

    out["source_study_effect_log"] = anova(lp)
    out["source_study_effect_ratio"] = anova([r["pt"] for r in rows])
    try:
        from scipy.stats import f as fdist
        for key in ("source_study_effect_log", "source_study_effect_ratio"):
            e = out[key]
            e["p"] = float(1 - fdist.cdf(e["F"], *e["df"]))
    except Exception:
        pass
    print(f"  source study: eta2 {out['source_study_effect_ratio']['eta_squared']} on the "
          f"ratio scale, {out['source_study_effect_log']['eta_squared']} on the log scale, "
          f"null expectation {out['source_study_effect_log']['eta_squared_expected_under_null']}")

    # objection two: is retention the unadjusted effect size renamed?
    for label, sel in (("all", lambda s: True),
                       ("rising", lambda s: s["unadjusted_mean_delta"] > 0),
                       ("falling", lambda s: s["unadjusted_mean_delta"] < 0)):
        sub = [s for s in ok if sel(s)]
        r, p, ci = spearman([abs(s["unadjusted_mean_delta"]) for s in sub],
                            [s["retention_joint"] for s in sub])
        out.setdefault("retention_vs_absolute_shift", {})[label] = {
            "n": len(sub), "spearman_rho": r, "spearman_p": p, "ci95": ci}
        print(f"  {label:8s} n={len(sub):3d}  rho={r:+.4f}  P={p:.3f}  CI {ci}")

    # Does the spread survive the negative control? §6ae reported only medians
    # and counts for the random pools, so the manuscript's claim that no control
    # removes the between-signature spread had no number behind it. S9 carries,
    # for each signature, the median retention across the 200 matched random
    # draws, which is exactly the comparison needed.
    import csv as _csv
    s9 = os.path.join(os.path.dirname(IR_RESULTS), "..", "SF2build",
                      "SupplementaryFile2", "SupplementaryTable_S9_negative_control.csv")
    if not os.path.exists(s9):
        s9 = os.path.join(IR_RESULTS, "..", "..", "SF2build", "SupplementaryFile2",
                          "SupplementaryTable_S9_negative_control.csv")
    if os.path.exists(s9):
        def f(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return None
        pr = [(f(r["retention_joint_D1"]), f(r["retention_random_matched_median"]))
              for r in _csv.DictReader(open(s9, encoding="utf-8"))]
        pr = [(a, b) for a, b in pr if a and b and a > 0 and b > 0]

        def qq(v, x):
            v = sorted(v)
            k = (len(v) - 1) * x
            i = int(k)
            j = min(i + 1, len(v) - 1)
            return v[i] + (k - i) * (v[j] - v[i])

        rec = {"n": len(pr)}
        for lab, vals in (("D1", [a for a, _ in pr]),
                          ("matched_random", [b for _, b in pr])):
            lv = [math.log(v) for v in vals]
            rec[lab] = {"p05": round(qq(vals, .05), 4), "p95": round(qq(vals, .95), 4),
                        "fold_p05_p95": round(qq(vals, .95) / qq(vals, .05), 1),
                        "variance_of_log": round(st.variance(lv), 4),
                        "sd_fold": round(math.exp(math.sqrt(st.variance(lv))), 2)}
        out["spread_under_the_negative_control"] = rec
        print(f"  spread under D1: {rec['D1']['fold_p05_p95']}-fold; under a matched "
              f"random covariate: {rec['matched_random']['fold_p05_p95']}-fold "
              f"(log variance {rec['D1']['variance_of_log']} vs "
              f"{rec['matched_random']['variance_of_log']})")

    # robust spread, since the 558-fold range rests on one point estimate
    pts = sorted(s["retention_joint"] for s in ok)

    def q(v, x):
        k = (len(v) - 1) * x
        f = int(k)
        c = min(f + 1, len(v) - 1)
        return v[f] + (k - f) * (v[c] - v[f])

    out["spread"] = {"min": pts[0], "second_smallest": pts[1], "max": pts[-1],
                     "fold_full": round(pts[-1] / pts[0], 1),
                     "fold_excluding_smallest": round(pts[-1] / pts[1], 1),
                     "p05": round(q(pts, .05), 3), "p95": round(q(pts, .95), 3),
                     "fold_p05_p95": round(q(pts, .95) / q(pts, .05), 1)}
    print(f"  spread: {out['spread']['fold_full']}-fold overall, "
          f"{out['spread']['fold_excluding_smallest']}-fold without the smallest, "
          f"{out['spread']['fold_p05_p95']}-fold between the 5th and 95th percentiles")

    path = os.path.join(IR_RESULTS, "SPREAD_IS_NOT_6v.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
