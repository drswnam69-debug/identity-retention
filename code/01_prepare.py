#!/usr/bin/env python3
"""01_prepare.py -- load one cohort, QC it, and write a tidy expression matrix.

Example
-------
python3 code/01_prepare.py \
    --counts data/GSE135251_RAW.tar \
    --series data/GSE135251_series_matrix.txt.gz \
    --cohort GSE135251 \
    --outdir results
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io_utils          # noqa: E402
import phenotype as ph   # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--series", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--annotation", default=None,
                    help="two-column Ensembl gene_id -> symbol map")
    ap.add_argument("--array", action="store_true",
                    help="input is a log-scale microarray matrix, not counts")
    ap.add_argument("--units", default="counts",
                    choices=["counts", "fpkm", "logged", "array-linear"],
                    help="counts: log2(CPM+1). fpkm: the same arithmetic, "
                         "which converts FPKM to TPM before the log. "
                         "array-linear: log2(x+1) then quantile "
                         "normalization. logged: already log-scaled.")
    ap.add_argument("--stage-source", default="auto",
                    choices=["auto", "group", "fibrosis", "subtype", "tissue", "tissue-hcc"],
                    help="which ordered severity axis this cohort supports "
                         "(see PREREGISTRATION.md, Replication axes)")
    args = ap.parse_args()

    out = os.path.join(args.outdir, args.cohort)
    os.makedirs(out, exist_ok=True)

    print(f"[{args.cohort}] reading phenotype ...")
    pheno = io_utils.read_series_matrix(args.series)

    print(f"[{args.cohort}] reading expression ...")
    mat = io_utils.read_counts(args.counts)
    mat = io_utils.drop_summary_rows(mat)

    # align samples: GEO titles sometimes replace accessions as column names
    shared = [s for s in pheno.index if s in mat.columns]
    if len(shared) < 0.5 * len(pheno):
        title_map = {}
        if "title" in pheno.columns:
            title_map = {str(t): g for g, t in pheno["title"].items()}
        renamed = {c: title_map[c] for c in mat.columns if c in title_map}
        if renamed:
            mat = mat.rename(columns=renamed)
            shared = [s for s in pheno.index if s in mat.columns]
    if not shared:
        raise SystemExit(
            "No overlap between series-matrix accessions and expression "
            f"columns.\n  pheno[:3] = {list(pheno.index[:3])}"
            f"\n  expr[:3]  = {list(mat.columns[:3])}")
    if len(shared) < len(pheno):
        print(f"  ! {len(pheno) - len(shared)} phenotype rows have no "
              f"expression column and were dropped")
    pheno, mat = pheno.loc[shared], mat[shared]

    # Normalize BEFORE mapping identifiers: library size must be the total
    # over all measured genes. Mapping first would compute CPM over only the
    # mapped subset and silently rescale every sample by a different factor.
    if args.units == "array-linear":
        expr = np.log2(mat.clip(lower=0) + 1.0)
        expr = io_utils.quantile_normalize(expr)
        print("  log2(x+1) then quantile normalization across samples")
    elif args.array or args.units == "logged" or io_utils.already_logged(mat):
        expr = mat
        print("  values taken as already log-scaled")
    else:
        expr = io_utils.log2_cpm(mat)
        print("  log2(CPM+1)" if args.units == "counts" else
              "  rescaled each sample to a fixed total, i.e. FPKM -> TPM, "
              "then log2(TPM+1)")
    expr = io_utils.to_symbols(expr, annotation=args.annotation)
    print(f"  expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

    pheno = pheno.copy()
    pheno["stage"], stage_order, axis = ph.build_stage_column(
        pheno, axis=args.stage_source)
    print(f"  severity axis: {axis}  ({' < '.join(stage_order)})")
    for name, cands in (("nas", ("nas_score", "nafld_activity_score")),
                        ("fibrosis", ("fibrosis_stage", "fibrosis")),
                        ("age", ("age", "age_at_biopsy")),
                        ("bmi", ("bmi",))):
        v = ph.numeric(pheno, *cands)
        if v is not None:
            pheno[name] = v

    counts_by_stage = (pd.Series(pheno["stage"]).value_counts()
                       .reindex(stage_order).fillna(0).astype(int))
    print("\n  stage counts (verify against the publication):")
    for k, v in counts_by_stage.items():
        print(f"    {k:<12} {v}")
    if pd.isna(pheno["stage"]).any():
        print(f"    {'UNASSIGNED':<12} {int(pd.isna(pheno['stage']).sum())}")

    # --- QC -----------------------------------------------------------------
    lib = mat.sum(axis=0)
    centred = expr.sub(expr.mean(axis=1), axis=0)
    sv = np.linalg.svd(centred.values, full_matrices=False)
    pcs = pd.DataFrame(sv[2][:4].T * sv[1][:4], index=expr.columns,
                       columns=[f"PC{i}" for i in range(1, 5)])
    var_expl = (sv[1] ** 2 / (sv[1] ** 2).sum())[:4]
    z_pc1 = (pcs["PC1"] - pcs["PC1"].mean()) / pcs["PC1"].std(ddof=1)
    outliers = list(pcs.index[np.abs(z_pc1) > 4])

    sex_genes = [g for g in ("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D")
                 if g in expr.index]
    sex = expr.loc[sex_genes].T if sex_genes else pd.DataFrame(index=expr.columns)

    expr.to_csv(os.path.join(out, "expr_log.tsv.gz"), sep="\t")
    pheno.to_csv(os.path.join(out, "phenotype.tsv"), sep="\t")
    pcs.join(sex).assign(library_size=lib).to_csv(
        os.path.join(out, "qc.tsv"), sep="\t")

    qc = {"cohort": args.cohort, "axis": axis, "stage_order": stage_order,
          "n_samples": int(expr.shape[1]),
          "n_genes": int(expr.shape[0]),
          "stage_counts": {k: int(v) for k, v in counts_by_stage.items()},
          "n_stage_unassigned": int(pd.isna(pheno["stage"]).sum()),
          "pc_variance_explained": [round(float(v), 4) for v in var_expl],
          "library_size_median": float(lib.median()),
          "library_size_min": float(lib.min()),
          "pc1_outliers_z_gt_4": outliers,
          "sex_marker_genes_found": sex_genes}
    with open(os.path.join(out, "qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)

    print(f"\n  PC1-4 variance explained: "
          f"{', '.join(f'{v:.1%}' for v in var_expl)}")
    print(f"  PC1 outliers (|z|>4): {outliers or 'none'}")
    print(f"  wrote -> {out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
