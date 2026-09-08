#!/usr/bin/env python3
"""33_jiang_prepare.py -- build the Jiang 2019 protein matrix for PREREG 6s.

Source: PRIDE PXD006512, MaxQuant proteinGroups.txt (CC0).

The processing follows the normalization the source paper describes for its own
analysis: iBAQ intensities are extracted per sample, quantile normalized, then
log2-transformed. Reverse hits, potential contaminants, and groups identified
only by site are removed first, as MaxQuant output requires.

The paper analyzed 98 non-tumor and 101 tumor samples; the deposited MaxQuant
run carries more, so the pairing here is derived from the sample identifiers and
the number of usable pairs is reported rather than assumed.
"""
from __future__ import annotations

import csv
import sys

csv.field_size_limit(10_000_000)

import numpy as np
import pandas as pd

SRC = ("/home/claude/jiang/MaxQuant results and its relative supplementary "
       "materials/proteinGroups.txt")
OUT = "/home/claude/jiang_proteins.tsv"


def main() -> None:
    src = sys.argv[1] if len(sys.argv) > 1 else SRC
    with open(src, newline="", encoding="latin-1") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        idx = {c: i for i, c in enumerate(hdr)}
        ibaq = [c for c in hdr
                if c.startswith("iBAQ ") and c != "iBAQ peptides"]
        need = {k: idx[k] for k in ("Gene names", "Only identified by site",
                                    "Reverse", "Potential contaminant")
                if k in idx}
        print(f"  {len(hdr)} columns, {len(ibaq)} iBAQ sample columns")
        rows, genes = [], []
        n_all = n_drop = 0
        for r in rdr:
            n_all += 1
            if any(r[need[k]].strip() == "+" for k in need if k != "Gene names"):
                n_drop += 1
                continue
            g = r[need["Gene names"]].split(";")[0].strip()
            if not g:
                continue
            genes.append(g)
            rows.append([r[idx[c]] for c in ibaq])
        print(f"  {n_all} protein groups, {n_drop} removed as reverse, "
              f"contaminant or site-only, {len(rows)} carry a gene name")

    mat = pd.DataFrame(rows, index=genes, columns=[c[5:] for c in ibaq])
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    mat = mat.groupby(level=0).max()
    print(f"  {mat.shape[0]} gene symbols x {mat.shape[1]} samples")

    # drop samples and proteins that are entirely zero before normalizing
    mat = mat.loc[:, mat.sum(axis=0) > 0]
    keep = (mat > 0).sum(axis=1)
    mat = mat[keep > 0]
    print(f"  after dropping all-zero rows and columns: {mat.shape}")

    # Quantile normalization across samples, then log2, as the source paper did.
    # Zeros are not measurements: MaxQuant reports a protein group as 0 in a
    # sample where it was not quantified. Ranking them with method="first" would
    # spread tied zeros across distinct low values and invent structure, so ties
    # take the average rank and every zero in a sample maps to the same value.
    zero_frac = float((mat == 0).values.mean())
    print(f"  zeros (not quantified) make up {zero_frac*100:.1f}% of the matrix")
    ranks = mat.rank(method="average", axis=0)
    mean_sorted = np.sort(mat.values, axis=0).mean(axis=1)
    n = mat.shape[0]
    grid = np.arange(1, n + 1)
    qn = pd.DataFrame(
        np.column_stack([np.interp(ranks[c].values, grid, mean_sorted)
                         for c in mat.columns]),
        index=mat.index, columns=mat.columns)
    qn = np.log2(qn + 1.0)
    print("  quantile normalized (average ranks for ties) then log2 transformed")

    tum = [c for c in qn.columns if c.endswith("T")]
    non = [c for c in qn.columns if c.endswith("P")]
    pids = sorted(set(c[:-1] for c in tum) & set(c[:-1] for c in non))
    print(f"  {len(tum)} tumor, {len(non)} non-tumor, {len(pids)} pairs")

    qn.index.name = "symbol"
    qn.to_csv(OUT, sep="\t")
    print(f"  wrote {OUT}")


if __name__ == "__main__":
    main()
