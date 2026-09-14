#!/usr/bin/env python3
"""48_gene_panel_archive.py -- archive the per-gene contrasts behind Table 4.

Table 4 quotes a paired tumor-minus-adjacent difference and a Wilcoxon P for
each of the seven module genes in three datasets, but only the two array
cohorts had those numbers written to a results file. The proteome columns and
the TCGA-LIHC per-gene values quoted in the Discussion existed only inside the
table, which means no automatic check could reach them. This script recomputes
them from the same matrices the original analyses read, with the same z scoring
and the same pairing, and writes them to results/GENE_PANEL_2026-09-09.json so
that the consistency checker can bind every cell of Table 4 to an archived
quantity.

This changes no reported result. It fills an archive gap.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

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
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, \
    DATA as IR_DATA, CODE as IR_CODE, FIGURES as IR_FIGURES, DOCS as IR_DOCS
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

RES = f"{IR_RESULTS}"

PANEL = ["CYB5R3", "CYB5R1", "AIFM2", "NQO1", "MTARC1", "MTARC2", "POR"]
ALIAS = {"MTARC1": "MARC1", "MTARC2": "MARC2", "G6PC1": "G6PC"}


def zrows(mat: pd.DataFrame) -> pd.DataFrame:
    sd = mat.std(axis=1)
    keep = mat[(sd > 0).values]
    return keep.sub(keep.mean(axis=1), axis=0).div(keep.std(axis=1), axis=0)


def contrast(z: pd.DataFrame, gene: str, pt, pn, pids):
    key = gene if gene in z.index else ALIAS.get(gene)
    if key is None or key not in z.index:
        return {"status": "not measured"}
    s = z.loc[key]
    if isinstance(s, pd.DataFrame):
        s = s.mean(axis=0)
    d = np.array([s[pt[p]] - s[pn[p]] for p in pids], float)
    d = d[~np.isnan(d)]
    if len(d) < 10:
        return {"status": f"only {len(d)} usable pairs"}
    return {"symbol_used": key,
            "n_pairs": int(len(d)),
            "paired_mean_delta": round(float(d.mean()), 4),
            "wilcoxon_p": float(wilcoxon(d, alternative="two-sided").pvalue)}


def proteome(path: str, label: str) -> dict:
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    df.index = [str(i).strip() for i in df.index]
    df = df[[c for c in df.columns if str(c).strip()]].apply(
        pd.to_numeric, errors="coerce")
    if float(np.nanmax(df.values)) > 100:
        df = np.log2(df.clip(lower=0) + 1.0)
    tum = [c for c in df.columns if str(c).startswith("T")]
    non = [c for c in df.columns if str(c).startswith("N")]
    pt = {c[1:]: c for c in tum}
    pn = {c[1:]: c for c in non}
    pids = sorted(set(pt) & set(pn))
    z = zrows(df)
    return {"source": label, "n_pairs_total": len(pids),
            "per_gene": {g: contrast(z, g, pt, pn, pids) for g in PANEL}}


def tcga_lihc() -> dict:
    """TCGA-LIHC paired tumor and adjacent normal, barcode-matched by patient."""
    df = pd.read_csv(f"{IR_DATA}/TCGA_LIHC_symbols.tsv.gz",
                     sep="\t", index_col=0, low_memory=False)
    df = df.apply(pd.to_numeric, errors="coerce")
    pt, pn = {}, {}
    for c in df.columns:
        parts = str(c).split("-")
        if len(parts) < 4:
            continue
        pid, typ = "-".join(parts[:3]), parts[3][:2]
        if typ == "01":
            pt.setdefault(pid, c)
        elif typ == "11":
            pn.setdefault(pid, c)
    pids = sorted(set(pt) & set(pn))
    z = zrows(df)
    return {"source": "TCGA-LIHC", "n_pairs_total": len(pids),
            "per_gene": {g: contrast(z, g, pt, pn, pids) for g in PANEL}}


def main() -> None:
    out = {"note": "archive gap fill; recomputed, no reported result changes",
           "Gao2019": proteome(f"{IR_DOCS}/gao_proteins.tsv", "Gao2019"),
           "Jiang2019": proteome(f"{IR_DOCS}/jiang_proteins_renamed.tsv",
                                 "Jiang2019"),
           "TCGA_LIHC": tcga_lihc()}
    for k in ("Gao2019", "Jiang2019", "TCGA_LIHC"):
        print(f"=== {k}  ({out[k]['n_pairs_total']} pairs) ===")
        for g, r in out[k]["per_gene"].items():
            if "paired_mean_delta" in r:
                print(f"  {g:<8} {r['paired_mean_delta']:+.4f}  "
                      f"P = {r['wilcoxon_p']:.3g}  (n = {r['n_pairs']})")
            else:
                print(f"  {g:<8} {r['status']}")
    with open(os.path.join(RES, "GENE_PANEL_2026-09-09.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote results/GENE_PANEL_2026-09-09.json")


if __name__ == "__main__":
    main()
