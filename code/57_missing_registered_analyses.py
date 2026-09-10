#!/usr/bin/env python3
"""57_missing_registered_analyses.py -- run two analyses §6d fixed and I skipped.

A pre-submission audit found that the manuscript explained two unmet
pre-registered commitments with reasons that were not true.

  1. §6d fixed recurrence-free survival as a reported secondary endpoint. The
     manuscript said GSE76427 carried overall survival only. The tidy
     phenotype table this study built carries only that, but the deposited
     series matrix carries event_rfs and duryears_rfs for the cohort, so the
     endpoint was skipped and not unavailable.
  2. §6d fixed H6 in every cohort already prepared, not only the prognostic
     one. GSE76427 was left off the command line, and the manuscript later
     explained its absence by a panel-restricted probe map. Its matrix carries
     20,757 symbols.

Both are run here, on the data that was there all along.
"""
from __future__ import annotations

import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

_here = os.path.dirname(os.path.abspath(__file__))
_c = [_here] + ([os.path.join(os.environ["IR_ROOT"], "code")]
                if os.environ.get("IR_ROOT") else [])
for _d in _c:
    if os.path.exists(os.path.join(_d, "paths.py")):
        sys.path.insert(0, _d)
        break
from paths import RESULTS, DATA, ROOT  # noqa: E402

sys.path.insert(0, _here)
import stats_lite as sl  # noqa: E402
import survival as sv  # noqa: E402

SERIES = f"{DATA}/GSE76427_series_matrix.txt.gz"
PANEL_OVERLAP = {"AIFM2", "GPX4", "NQO1", "SLC7A11", "ACSL4"}


def rfs_from_series() -> pd.DataFrame:
    """event_rfs and duryears_rfs, keyed by GSM, from the deposited matrix."""
    gsm, rows = None, {}
    with gzip.open(SERIES, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsm = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                if vals and ":" in vals[0]:
                    key = vals[0].split(":", 1)[0].strip()
                    if key in ("event_rfs", "duryears_rfs"):
                        rows[key] = [v.split(":", 1)[1].strip() if ":" in v else ""
                                     for v in vals]
    df = pd.DataFrame(rows, index=gsm)
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna()


def main() -> int:
    out: dict = {"plan": "§6d secondary endpoint and §6d H6, both run late",
                 "note": "the manuscript previously reported these as unmet"}

    # ---- 1. recurrence-free survival, the registered secondary endpoint -----
    rfs = rfs_from_series()
    rsi = pd.read_csv(f"{RESULTS}/GSE76427/rsi.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(f"{RESULTS}/GSE76427/phenotype.tsv", sep="\t", index_col=0)
    tum = ph[~ph["tissue"].astype(str).str.contains("adjacent", case=False)]
    ids = [i for i in tum.index if i in rfs.index and i in rsi.index]
    z = rsi.loc[ids, "RSI"].to_numpy(float)
    z = (z - z.mean()) / z.std(ddof=1)
    ev = rfs.loc[ids, "event_rfs"].to_numpy(int)
    tt = rfs.loc[ids, "duryears_rfs"].to_numpy(float)
    fit = sv.cox_ph(tt, ev, np.column_stack([z]), ["RSI_z"])
    b = float(fit["coef"][0]); se = float(fit["se"][0])
    out["rfs_secondary"] = {
        "endpoint": "recurrence-free survival, §6d secondary",
        "source": "GSE76427 series matrix characteristics, event_rfs / duryears_rfs",
        "n": len(ids), "events": int(ev.sum()),
        "hr_per_sd": round(float(np.exp(b)), 4),
        "ci95": [round(float(np.exp(b - 1.96 * se)), 4),
                 round(float(np.exp(b + 1.96 * se)), 4)],
        "p": float(fit["p"][0]),
    }
    r = out["rfs_secondary"]
    print(f"  RFS: n={r['n']}, {r['events']} events, HR {r['hr_per_sd']:.3f} "
          f"({r['ci95'][0]:.3f} to {r['ci95'][1]:.3f}), P = {r['p']:.3f}")

    # ---- 2. H6 in the cohort §6d made primary --------------------------------
    expr = pd.read_csv(f"{DATA}/GSE76427_symbols.tsv.gz", sep="\t", index_col=0)
    gmt = {}
    with open(f"{ROOT}/genesets/liver_cgp.gmt", encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gmt[p[0]] = [g for g in p[2:] if g]
    h6 = json.load(open(f"{RESULTS}/phase_e_h6.json", encoding="utf-8"))["H6"]
    sets = {k: v["dropped_overlapping_with_RSI_panel"] for k, v in h6["sets"].items()}
    out["h6_GSE76427"] = {"n_symbols_on_platform": int(expr.shape[0]), "sets": {}}
    for name in ("GOBP_NEGATIVE_REGULATION_OF_FERROPTOSIS",
                 "GOBP_POSITIVE_REGULATION_OF_FERROPTOSIS"):
        members = [g for g in gmt.get(name, []) if g not in PANEL_OVERLAP]
        if not members:
            src = f"{ROOT}/genesets/ferroptosis_{name}.txt"
            if os.path.exists(src):
                members = [l.strip() for l in open(src) if l.strip()
                           and l.strip() not in PANEL_OVERLAP]
        present = [g for g in members if g in expr.index]
        if len(present) < 5:
            out["h6_GSE76427"]["sets"][name] = {
                "n_members_after_overlap_removal": len(members),
                "n_on_platform": len(present),
                "status": "fewer than 5 members on this platform"}
            print(f"  H6 {name[:44]:46s} {len(present):3d} on platform -- skipped")
            continue
        sub = expr.loc[present, ids]
        zz = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1), axis=0).mean(axis=0)
        rho, p = sl.spearman(rsi.loc[ids, "RSI"].to_numpy(float), zz.to_numpy(float))
        out["h6_GSE76427"]["sets"][name] = {
            "n_members_after_overlap_removal": len(members),
            "n_on_platform": len(present), "n": len(ids),
            "spearman_rho": round(float(rho), 4), "p": float(p)}
        print(f"  H6 {name[:44]:46s} {len(present):3d} genes, rho {rho:+.3f}, "
              f"P = {p:.3g}")

    p = f"{RESULTS}/MISSING_REGISTERED_6d.json"
    json.dump(out, open(p, "w", encoding="utf-8"), indent=2)
    print("\nwritten", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
