#!/usr/bin/env python3
"""37_tcga_prepare.py -- PREREGISTRATION 6w.

Build a symbol-level matrix and the tumor/adjacent pairing for TCGA-LIHC from the
UCSC Xena GDC hub STAR - TPM file, under the rules fixed in 6w before this file
was opened:

  * source is TCGA-LIHC.star_tpm.tsv, unit log2(TPM + 1), aliquots already
    averaged by the hub. No alternative quantification is substituted.
  * symbols come from gencode.v36.annotation.gtf.gene.probemap.
  * where one symbol carries several Ensembl identifiers, the MAXIMUM across
    identifiers on the log2 scale is taken, the same collapse rule this study
    applied to multiple probes on the arrays.
  * the single alias G6PC1 -> G6PC is applied, since GENCODE v36 carries the
    retired symbol and G6PC1 is a member of the locked D1 list.
  * pairing is from the TCGA barcode: patient = first 12 characters, sample type
    = characters 14 and 15, 01 = primary tumor, 11 = solid tissue normal. Any
    other sample type is discarded before pairing. Where a patient still has more
    than one sample of a type, the first in sorted order is taken and the number
    of such patients is reported.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

ALIAS = {"G6PC1": "G6PC"}          # applied in reverse when indexing, see below
TUMOR_CODE, NORMAL_CODE = "01", "11"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--probemap", required=True)
    ap.add_argument("--out-expr", required=True)
    ap.add_argument("--out-pheno", required=True)
    ap.add_argument("--report", required=True)
    a = ap.parse_args()

    print(f"=== TCGA preparation (PREREG 6w for LIHC, 6x for KIRC): "
          f"{os.path.basename(a.matrix)} ===")
    pm = pd.read_csv(a.probemap, sep="\t")
    id2sym = dict(zip(pm["id"], pm["gene"]))
    print(f"  probemap: {len(pm)} identifiers, {pm['gene'].nunique()} symbols")

    df = pd.read_csv(a.matrix, sep="\t", index_col=0)
    print(f"  matrix:   {df.shape[0]} identifiers x {df.shape[1]} samples")
    lo, hi = float(np.nanmin(df.values)), float(np.nanmax(df.values))
    print(f"  value range {lo:.4g} to {hi:.4g}  (declared unit log2(TPM + 1))")
    if lo < -1e-9:
        raise SystemExit("negative values: this is not a log2(TPM + 1) matrix")
    n_na = int(df.isna().sum().sum())
    print(f"  missing values: {n_na}")

    sym = pd.Series([id2sym.get(i) for i in df.index], index=df.index)
    n_unmapped = int(sym.isna().sum())
    df = df[sym.notna().values]
    sym = sym[sym.notna()]
    print(f"  {n_unmapped} identifiers had no symbol and were dropped")

    dup_before = int(sym.duplicated().sum())
    df.index = sym.values
    expr = df.groupby(level=0).max()          # 6w collapse rule
    print(f"  {dup_before} duplicate symbol rows collapsed by max on the log2 "
          f"scale -> {expr.shape[0]} symbols")

    # 6w alias: the locked list says G6PC1, GENCODE v36 says G6PC
    renamed = []
    for current, legacy in ALIAS.items():
        if current not in expr.index and legacy in expr.index:
            expr = expr.rename(index={legacy: current})
            renamed.append(f"{legacy} -> {current}")
    print(f"  alias applied: {renamed if renamed else 'none needed'}")

    # --- barcode parsing, 6w rule -------------------------------------------
    recs = []
    other = {}
    for c in expr.columns:
        parts = str(c).split("-")
        if len(parts) < 4:
            other[c] = "malformed"
            continue
        patient = "-".join(parts[:3])
        code = parts[3][:2]
        if code == TUMOR_CODE:
            recs.append((c, patient, "tumor"))
        elif code == NORMAL_CODE:
            recs.append((c, patient, "adjacent"))
        else:
            other[code] = other.get(code, 0) + 1
    ph = pd.DataFrame(recs, columns=["sample", "patient_id", "tissue"]).set_index("sample")
    print(f"\n  sample types kept: "
          f"{int((ph.tissue=='tumor').sum())} tumor (01), "
          f"{int((ph.tissue=='adjacent').sum())} adjacent (11)")
    print(f"  sample types discarded: {other if other else 'none'}")

    multi = 0
    keep = []
    for (pid, tis), grp in ph.groupby(["patient_id", "tissue"]):
        s = sorted(grp.index)
        if len(s) > 1:
            multi += 1
        keep.append(s[0])
    ph = ph.loc[sorted(keep)]
    print(f"  patients with more than one sample of a type: {multi} "
          f"(first in sorted order kept)")

    t = set(ph[ph.tissue == "tumor"]["patient_id"])
    n = set(ph[ph.tissue == "adjacent"]["patient_id"])
    paired = sorted(t & n)
    print(f"  patients with tumor: {len(t)}; with adjacent: {len(n)}; "
          f"PAIRED: {len(paired)}")

    ph = ph[ph["patient_id"].isin(paired)]
    expr = expr[ph.index]
    print(f"  final matrix: {expr.shape[0]} symbols x {expr.shape[1]} samples")

    os.makedirs(os.path.dirname(a.out_pheno), exist_ok=True)
    expr.to_csv(a.out_expr, sep="\t")
    ph.to_csv(a.out_pheno, sep="\t")
    # the benchmark scripts read results/<cohort>/rsi.tsv only for its index
    pd.DataFrame(index=ph.index).assign(placeholder=0.0).to_csv(
        os.path.join(os.path.dirname(a.out_pheno), "rsi.tsv"), sep="\t")

    rep = {"plan": "PREREGISTRATION 6w",
           "amendment_sha256_of_text_as_written":
               "f854c263bcffbdd80c249bbcbf59558480cd0b047f2073910e45be20a262498c",
           "source_file": os.path.basename(a.matrix),
           "declared_unit": "log2(TPM + 1)",
           "value_range": [round(lo, 4), round(hi, 4)],
           "n_missing_values": n_na,
           "n_identifiers_in": int(df.shape[0]) + n_unmapped,
           "n_unmapped_dropped": n_unmapped,
           "n_duplicate_rows_collapsed": dup_before,
           "collapse_rule": "max across Ensembl identifiers on the log2 scale",
           "alias_applied": renamed,
           "n_symbols": int(expr.shape[0]),
           "n_tumor_01": int((ph.tissue == "tumor").sum()),
           "n_adjacent_11": int((ph.tissue == "adjacent").sum()),
           "sample_types_discarded": {str(k): v for k, v in other.items()},
           "n_patients_multi_sample": multi,
           "n_pairs": len(paired)}
    with open(a.report, "w") as fh:
        json.dump(rep, fh, indent=1)
    print(f"\n  wrote {a.out_expr}\n  wrote {a.out_pheno}\n  wrote {a.report}")


if __name__ == "__main__":
    main()
