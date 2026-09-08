#!/usr/bin/env python3
"""25_hpa_celltype.py -- PREREGISTRATION 6m.

6m-A: cell-type attribution of the RSI panel genes in liver, from the Human
Protein Atlas single-cell data. Discharges the 6j-S obligation.
6m-B: protein-level (IHC) corroboration in normal liver, same genes.

The decision rule was fixed in 6m before any value here was read:
  hepatocyte fraction = hepatocyte nCPM / sum of nCPM over liver cell types
  >= 0.50 for the majority of genes in a module -> premise supported for it.

This is SUPPORTING EVIDENCE ABOUT THE GENES. It cannot change any verdict in
6g-6l, and nothing here is fed back into the index.
"""
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rsi_config as cfg

REDUCTION = cfg.MODULE_REDUCTION          # CYB5R3 CYB5R1 AIFM2 NQO1
DRAIN = cfg.MODULE_DRAIN                  # MTARC1 MTARC2 POR
MODULE = {g: "REDUCTION" for g in REDUCTION}
MODULE.update({g: "DRAIN" for g in DRAIN})

THRESHOLD = 0.50                          # fixed in 6m


def main():
    sc = pd.read_csv("data/hpa/hpa_liver_singlecell.csv")
    ihc = pd.read_csv("data/hpa/hpa_liver_protein_ihc.csv")

    assert set(sc["gene"]) == set(MODULE), "6m covers exactly the 7 panel genes"
    n_types = sc.groupby("gene")["cell_type"].nunique()
    assert n_types.nunique() == 1, "cell-type set must be identical across genes"

    rows = []
    for g, sub in sc.groupby("gene"):
        tot = float(sub["nCPM"].sum())
        hep = float(sub.loc[sub["cell_type"] == "Hepatocytes", "nCPM"].iloc[0])
        # as registered: unweighted share of summed per-cell-type nCPM
        frac = hep / tot if tot > 0 else float("nan")
        # secondary, descriptive only: weighted by the atlas's own cell counts
        sub = sub.assign(w=sub["nCPM"] * sub["n_cells"])
        wfrac = float(sub.loc[sub["cell_type"] == "Hepatocytes", "w"].iloc[0]) / \
            float(sub["w"].sum()) if sub["w"].sum() > 0 else float("nan")
        top = sub.sort_values("nCPM", ascending=False).iloc[0]
        rows.append({"gene": g, "module": MODULE[g], "hep_nCPM": hep,
                     "sum_nCPM": round(tot, 1),
                     "hep_fraction": round(frac, 3),
                     "hep_fraction_weighted": round(wfrac, 3),
                     "top_cell_type": top["cell_type"],
                     "top_nCPM": float(top["nCPM"]),
                     "hep_is_top": bool(top["cell_type"] == "Hepatocytes")})
    df = pd.DataFrame(rows).sort_values(["module", "gene"])

    print("=== 6m-A: hepatocyte share of liver single-cell expression ===")
    print(df.to_string(index=False))

    verdict = {}
    for mod in ("REDUCTION", "DRAIN"):
        sub = df[df["module"] == mod]
        n_pass = int((sub["hep_fraction"] >= THRESHOLD).sum())
        n = len(sub)
        supported = n_pass > n / 2
        verdict[mod] = {"genes": n, "n_at_or_above_threshold": n_pass,
                        "majority_supported": bool(supported)}
        print(f"\n  {mod}: {n_pass}/{n} genes with hepatocyte fraction >= "
              f"{THRESHOLD:.2f}  ->  premise "
              f"{'SUPPORTED' if supported else 'NOT SUPPORTED'} for this module")

    # --- 6m-B and the dropout control -------------------------------------
    print("\n=== 6m-B: protein (IHC) in normal liver, and the dropout control ===")
    m = df.merge(ihc, on=["gene", "module"])
    m["rna_says"] = m["hep_fraction"].map(lambda f: "low" if f < 0.15 else
                                          ("mid" if f < 0.5 else "high"))
    DETECTED = {"Low", "Medium", "High"}
    m["protein_detected_hep"] = m["hepatocytes_IHC"].isin(DETECTED)
    m["discordant"] = m["protein_detected_hep"] & (m["hep_fraction"] < 0.15)
    print(m[["gene", "module", "hep_nCPM", "hep_fraction", "hepatocytes_IHC",
             "cholangiocytes_IHC", "discordant"]].to_string(index=False))
    n_disc = int(m["discordant"].sum())
    print(f"\n  genes with hepatocyte protein detected but hepatocyte RNA "
          f"fraction < 0.15: {n_disc}")

    # DRAIN hepatocyte-identity check at protein level
    order = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
    dr = m[m["module"] == "DRAIN"]
    hep_gt_chol = int(sum(order[a] > order[b] for a, b in
                          zip(dr["hepatocytes_IHC"], dr["cholangiocytes_IHC"])))
    print(f"  DRAIN genes with hepatocyte staining > cholangiocyte staining: "
          f"{hep_gt_chol}/{len(dr)}")
    red = m[m["module"] == "REDUCTION"]
    n_red_nd = int((red["hepatocytes_IHC"] == "Not detected").sum())
    print(f"  REDUCTION genes not detected in normal hepatocytes: "
          f"{n_red_nd}/{len(red)}")

    out = {"plan": "PREREGISTRATION 6m",
           "source": "Human Protein Atlas, liver single-cell and normal-liver IHC",
           "threshold": THRESHOLD,
           "per_gene": df.to_dict("records"),
           "module_verdict_6mA": verdict,
           "protein_ihc": ihc.to_dict("records"),
           "rna_protein_discordant_genes": n_disc,
           "drain_hep_gt_cholangiocyte": f"{hep_gt_chol}/{len(dr)}",
           "reduction_not_detected_normal_hepatocyte": f"{n_red_nd}/{len(red)}",
           "6mB_tumour_vs_normal_IHC": "NOT PERFORMED - the tumour staining "
           "table is not extractable from HPA's current page structure; 6m "
           "forbids substituting another source"}
    json.dump(out, open("results/HPA_celltype_6m.json", "w"), indent=1)
    print("\nwrote results/HPA_celltype_6m.json")


if __name__ == "__main__":
    main()
