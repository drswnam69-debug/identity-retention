#!/usr/bin/env python3
"""test_retentionfrac.py -- exercise the deposited package against this study.

The manuscript makes four behavioral claims about `retentionfrac`: the premise
check and the ratio-stability guard live inside the function, it agrees with
this study's own archived module value, and it refuses a reversed pairing.
Until this file existed, nothing in the archive tested the package at all, and
the "agrees exactly" claim went three revisions without being checked. It does
not agree exactly, for a reason worth recording: the pipeline divides an
intercept that `ols_ci` has already rounded to four decimals, so the archived
fraction and a fraction computed from the raw intercept can differ in the
fourth place. That tolerance is asserted here rather than hidden.

Run:  python3 code/test_retentionfrac.py
"""
from __future__ import annotations

import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))
_cands = [_here, _os.path.join(_here, "code")]
if _os.environ.get("IR_ROOT"):
    _cands.insert(0, _os.path.join(_os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        if _c not in _sys.path:
            _sys.path.insert(0, _c)
        break
from paths import ROOT, RESULTS, DATA  # noqa: E402

import json  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_sys.path.insert(0, _os.path.join(ROOT, "tool"))
from retentionfrac import identity_retention  # noqa: E402

import rsi_config as cfg  # noqa: E402

D1 = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2",
      "CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1", "ASGR1", "HNF4A", "HNF1A",
      "FOXA1", "FOXA2", "NR1H4"]
TUMOR = "HCC tumor"
OK, BAD = [], []


def check(cond: bool, label: str) -> None:
    (OK if cond else BAD).append(label)
    print(("  PASS  " if cond else "  FAIL  ") + label)


def load():
    expr = pd.read_csv(f"{DATA}/GSE14520_symbols.tsv.gz", sep="\t", index_col=0)
    ph = pd.read_csv(f"{RESULTS}/GSE14520/phenotype.tsv", sep="\t", index_col=0)
    t = ph[ph.tissue == TUMOR]
    a = ph[ph.tissue != TUMOR]
    pt = dict(zip(t.patient_id, t.index))
    pn = dict(zip(a.patient_id, a.index))
    ids = sorted(set(pt) & set(pn))
    return expr, [pt[i] for i in ids], [pn[i] for i in ids]


def main() -> int:
    print("=== retentionfrac against this study's own archive ===")
    expr, tum, adj = load()
    print(f"  GSE14520, {len(tum)} pairs, {expr.shape[0]} symbols")

    arch = json.load(open(f"{RESULTS}/SIGNATURE_BENCHMARK_GSE14520.json",
                          encoding="utf-8"))
    own = arch.get("own_modules", {})

    # 1. the study's own modules, against the archived identity-only fractions
    for mod in ("DRAIN", "REDUCTION"):
        genes = getattr(cfg, f"MODULE_{mod}")
        r = identity_retention(expr, tum, adj, genes, D1, min_genes=2, n_boot=0)
        want = own.get(mod, {}).get("retention_D1")
        if want is None:
            check(False, f"{mod}: archive carries an identity-only fraction")
            continue
        check(r.retention is not None and abs(r.retention - want) <= 2e-4,
              f"{mod}: package {r.retention:.6f} agrees with archived {want} "
              f"to four decimals (difference {abs(r.retention - want):.2e})")

    # 2. the premise check fires on a reversed pairing
    r = identity_retention(expr, adj, tum, cfg.MODULE_DRAIN, D1,
                           min_genes=2, n_boot=0)
    check(r.retention is None and not r.premise_ok,
          f"a reversed pairing is refused (premise delta {r.premise_delta:+.4f})")

    # 3. the premise check fires on a matrix where identity rises in tumor
    rng = np.random.default_rng(7)
    idx = list(expr.index[:400]) + D1 + list(cfg.MODULE_DRAIN)
    idx = list(dict.fromkeys(idx))
    m = pd.DataFrame(rng.normal(size=(len(idx), 40)), index=idx,
                     columns=[f"s{i}" for i in range(40)])
    tt, aa = list(m.columns[:20]), list(m.columns[20:])
    m.loc[[g for g in D1 if g in m.index], tt] += 2.0      # identity UP in tumor
    r = identity_retention(m, tt, aa, [g for g in cfg.MODULE_DRAIN if g in m.index],
                           [g for g in D1 if g in m.index], min_genes=2, n_boot=0)
    check(r.retention is None and not r.premise_ok,
          "a matrix where identity rises in tumor is refused")

    # 4. the ratio-stability guard fires when the denominator spans zero
    m2 = pd.DataFrame(rng.normal(size=(len(idx), 60)), index=idx,
                      columns=[f"s{i}" for i in range(60)])
    tt2, aa2 = list(m2.columns[:30]), list(m2.columns[30:])
    m2.loc[[g for g in D1 if g in m2.index], tt2] -= 2.0   # premise holds
    sig = [g for g in m2.index[:30]]                        # no real shift
    r = identity_retention(m2, tt2, aa2, sig, [g for g in D1 if g in m2.index],
                           min_genes=2, n_boot=200)
    check(r.retention is None or not r.stable,
          "a signature with no shift to retain returns no fraction")

    # 5. the guard is not silently bypassable: the reason is always given
    check(bool(r.notes), "a refusal carries a stated reason")

    print(f"\n{len(OK)} of {len(OK) + len(BAD)} checks passed")
    if BAD:
        print("FAILED:")
        for b in BAD:
            print("  - " + b)
    return 1 if BAD else 0


if __name__ == "__main__":
    raise SystemExit(main())
