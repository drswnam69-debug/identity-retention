#!/usr/bin/env python3
"""31_prepare_full_matrix.py -- full symbol-level matrix for PREREG 6t.

The archived pipeline kept only the panel genes, so the signature benchmark
needs the full matrix rebuilt from the GEO series matrix. The preprocessing
here reproduces Methods 2.2 of the manuscript exactly: array data are
log2-transformed and quantile-normalized across the full probe matrix BEFORE
mapping to symbols, and the maximum is taken where several probes map to one
symbol.

Whether this reproduction is faithful is not asserted, it is tested: the
benchmark script runs this study's own REDUCTION and DRAIN modules through the
identical code path and aborts unless it recovers the archived retention values.

Usage:
  python3 code/31_prepare_full_matrix.py \
      --series data/GSE14520-GPL3921_series_matrix.txt.gz \
      --annot  data/GPL3921.annot.gz \
      --units array-log \
      --out    data/GSE14520_symbols.tsv.gz
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import io_utils  # noqa: E402


def read_matrix(path: str) -> pd.DataFrame:
    """The value table of a GEO series matrix file."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        lines = fh.readlines()
    try:
        start = next(i for i, l in enumerate(lines)
                     if l.startswith("!series_matrix_table_begin")) + 1
    except StopIteration:
        raise SystemExit(f"{path}: no series_matrix_table_begin marker")
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("!series_matrix_table_end"):
            end = i
            break
    body = "".join(lines[start:end])
    mat = pd.read_csv(io.StringIO(body), sep="\t", index_col=0,
                      na_values=["", "null", "NA"])
    mat.index = mat.index.astype(str).str.strip('"')
    mat.columns = [c.strip('"') for c in mat.columns]
    return mat.apply(pd.to_numeric, errors="coerce")


def read_annot(path: str) -> pd.Series:
    """probe id -> gene symbol, from a GEO platform .annot or a SOFT table."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        lines = fh.readlines()
    start = 0
    for i, l in enumerate(lines):
        if l.startswith("!platform_table_begin"):
            start = i + 1
            break
    end = len(lines)
    for i in range(start, len(lines)):
        if lines[i].startswith("!platform_table_end"):
            end = i
            break
    tab = pd.read_csv(io.StringIO("".join(lines[start:end])), sep="\t",
                      dtype=str, low_memory=False)
    tab.columns = [c.strip().strip('"') for c in tab.columns]
    idcol = tab.columns[0]
    cand = [c for c in tab.columns
            if c.lower() in ("gene symbol", "symbol", "ilmn_gene",
                             "gene_symbol", "genesymbol")]
    if not cand:
        raise SystemExit(f"{path}: no gene symbol column; saw {list(tab.columns)[:12]}")
    col = cand[0]
    print(f"  annotation: {len(tab)} rows, id column '{idcol}', "
          f"symbol column '{col}'")
    s = tab.set_index(idcol)[col].dropna().astype(str).str.strip()
    # GEO writes multi-symbol entries as "A///B"; the first is the primary.
    s = s.str.split("///").str[0].str.strip()
    s = s[(s != "") & (s.str.lower() != "nan")]
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--units", choices=["array-log", "array-linear"],
                    required=True,
                    help="array-log: values are already log2, quantile "
                         "normalize only. array-linear: log2(x+1) first.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-quantile", action="store_true",
                    help="skip quantile normalization. GSE14520 was analyzed "
                         "from the RMA log2 values as served by GEO, without "
                         "an added quantile step; using one changes every "
                         "value and breaks reproduction of the archived run.")
    ap.add_argument("--alias-annotation", default=None,
                    help="NCBI gene annotation table used to rename legacy "
                         "symbols to current HGNC symbols, as in the archived "
                         "pipeline (io_utils.resolve_symbol_aliases)")
    a = ap.parse_args()

    mat = read_matrix(a.series)
    print(f"  {a.series}: {mat.shape[0]} probes x {mat.shape[1]} samples")
    before = mat.shape[0]
    mat = mat.dropna(how="all")
    mat = mat.loc[mat.notna().all(axis=1)]
    if mat.shape[0] != before:
        print(f"  dropped {before - mat.shape[0]} probes with missing values")

    logged = io_utils.already_logged(mat)
    print(f"  already_logged() says {logged}; --units {a.units}")
    if a.units == "array-linear":
        if logged:
            print("  !! values look logged already but --units array-linear "
                  "was given; proceeding as instructed")
        mat = np.log2(mat.clip(lower=0) + 1.0)
    else:
        if not logged:
            print("  !! values do not look logged but --units array-log "
                  "was given; proceeding as instructed")
    if a.no_quantile:
        print("  quantile normalization SKIPPED (--no-quantile)")
    else:
        mat = io_utils.quantile_normalize(mat)
        print("  quantile normalized across the full probe matrix")

    ann = read_annot(a.annot)
    common = mat.index.intersection(ann.index)
    print(f"  {len(common)} of {mat.shape[0]} probes carry a symbol")
    if len(common) < 0.5 * mat.shape[0]:
        raise SystemExit("fewer than half the probes mapped; wrong annotation "
                         "file for this platform?")
    mat = mat.loc[common]
    mat.index = ann.loc[common].values
    n_probes = mat.shape[0]
    mat = mat.groupby(level=0).max()
    print(f"  {n_probes} probes collapsed to {mat.shape[0]} symbols "
          f"(maximum across probes, as in Methods 2.2)")

    if a.alias_annotation:
        before = set(mat.index)
        mat = io_utils.resolve_symbol_aliases(mat, a.alias_annotation)
        renamed = len(before - set(mat.index))
        print(f"  legacy symbol aliases resolved; {renamed} symbols renamed")

    mat.index.name = "symbol"
    mat.to_csv(a.out, sep="\t", compression="gzip" if a.out.endswith(".gz") else None)
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
