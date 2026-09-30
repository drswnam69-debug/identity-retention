#!/usr/bin/env python3
"""79_manuscript_numbers.py

Merge every result file the new manuscript draws on into one JSON, so that the
text and the verifier read the same source. Adds nothing: every value here is
copied or arithmetically derived from a file written by 70, 75, 76, 77 or 78.
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = f"{ROOT}/results"


def load(name):
    with open(f"{R}/{name}") as fh:
        return json.load(fh)


def main() -> int:
    pc = load("PANEL_COHERENCE.json")
    st = load("SECRETOME_TRANSFER.json")
    sv = load("SECRETOME_VERIFY.json")
    yi = {s: load(f"HOLDOUT_Yi2023_{s}.json")
          for s in ("EncyclopeDIA", "Spectronaut", "DIANN")}
    ic = load("ICCA_VALIDATION.json")

    out = {"what": "every number the manuscript prints, assembled from the "
                   "result files; nothing is entered by hand"}

    # ---- the eight datasets, one row each
    rows = []
    for name, v in pc["datasets"].items():
        bc, cc = v["panel_by_class"], v["class_coherence"]
        met, sec = bc.get("metabolic enzyme", {}), bc.get("secreted plasma protein", {})
        rows.append({
            "dataset": name, "level": v["level"], "tumor_type": "HCC",
            "n_pairs": v["n_pairs"],
            "metabolic_n": met.get("n"), "metabolic_mean_z": met.get("mean_z"),
            "metabolic_share": (cc.get("metabolic enzyme") or {}).get("cancellation"),
            "secreted_n": sec.get("n"), "secreted_mean_z": sec.get("mean_z"),
            "secreted_share": (cc.get("secreted plasma protein") or {}).get("cancellation"),
            "source": "70_panel_coherence.py"})
    e = yi["EncyclopeDIA"]
    rows.append({
        "dataset": "Yi 2023 proteome (held out)", "level": "proteome",
        "tumor_type": "HCC", "n_pairs": e["premises"]["n_pairs"],
        "metabolic_n": e["H1_metabolic_module"]["n"],
        "metabolic_mean_z": e["H1_metabolic_module"]["mean_z"],
        "metabolic_share": e["H1_metabolic_module"]["coherence"]["cancellation"],
        "secreted_n": e["H2_plasma_proteins"]["n"],
        "secreted_mean_z": e["H2_plasma_proteins"]["mean_z"],
        "secreted_share": e["H2_plasma_proteins"]["coherence"]["cancellation"],
        "source": "77_holdout_validation.py"})
    for label, c in ic["cohorts"].items():
        rows.append({
            "dataset": f"{label} proteome (exploratory)", "level": "proteome",
            "tumor_type": "intrahepatic cholangiocarcinoma",
            "n_pairs": len([k for k in c["provenance"]]) and c["provenance"]["patients_with_both_arms"],
            "metabolic_n": c["G1_module"]["n"],
            "metabolic_mean_z": c["G1_module"]["mean_z"],
            "metabolic_share": c["G1_module"]["coherence"]["cancellation"],
            "secreted_n": c["G3_plasma"]["n"],
            "secreted_mean_z": c["G3_plasma"]["mean_z"],
            "secreted_share": c["G3_plasma"]["coherence"]["cancellation"],
            "source": "78_icca_validation.py"})
    out["datasets"] = rows
    out["n_datasets"] = len(rows)
    out["gap_secreted_minus_metabolic"] = {
        r["dataset"]: round(r["secreted_mean_z"] - r["metabolic_mean_z"], 4)
        for r in rows if None not in (r["secreted_mean_z"], r["metabolic_mean_z"])}

    # ---- the comparison sets and the original panel
    cbl = pc["coherence_by_level"]
    out["comparison_sets"] = {
        "n_sets": cbl["n_sets_scored_at_both_levels"],
        "median_share_transcriptome": cbl["cancellation_median_transcriptome"],
        "median_share_proteome": cbl["cancellation_median_proteome"],
        "mean_difference": cbl["cancellation_mean_difference"],
        "difference_ci95": cbl["cancellation_difference_ci95"],
        "wilcoxon_p": cbl["cancellation_wilcoxon_p"],
        "sign_agreement_median_transcriptome": cbl["sign_agreement_median_transcriptome"],
        "sign_agreement_median_proteome": cbl["sign_agreement_median_proteome"],
        "panel_drop": pc["panel_within_shared_sets"],
        "source": "70_panel_coherence.py"}

    # ---- the secretome question, which came out null
    out["secretome_question"] = {
        "h1_drop_vs_secreted_fraction":
            st["h1_drop_vs_secreted_fraction"],
        "gene_level_secreted_minus_other": {
            k: {"level": v["level"], "n_genes_scored": v["n_genes_scored"],
                "mean_secreted": v.get("mean_secreted"),
                "mean_other": v.get("mean_other"),
                "difference_of_means": v.get("difference_of_means"),
                "mannwhitney_p": v.get("mannwhitney_p")}
            for k, v in sv["A_gene_level"].items()},
        "abundance_matched": sv["B_abundance_matched"],
        "set_overlap": {k: v for k, v in sv["D_set_overlap"].items()
                        if k != "gene_disjoint_subfamily"},
        "n_sets_in_gene_disjoint_subfamily":
            len(sv["D_set_overlap"]["gene_disjoint_subfamily"]),
        "annotation_check": st["annotation_check"],
        "source": "75_secretome_transfer.py and 76_secretome_verify.py"}

    # ---- the held-out test, all three software outputs
    out["held_out_test"] = {}
    for s, d in yi.items():
        row = {"gate": d["premises"]["sample_identity_gate"],
               "n_pairs": d["premises"]["n_pairs"],
               "n_genes": d["premises"]["n_genes"],
               "missingness": d["premises"]["missingness"],
               "matrix_diagnostic": d["premises"]["matrix_diagnostic"]}
        for k in ("H1_metabolic_module", "H2_plasma_proteins", "H3_extremity",
                  "H4_composite", "H5_abundance"):
            if k in d:
                row[k] = d[k]
        row["verdict"] = d.get("verdict")
        out["held_out_test"][s] = row
    out["held_out_primary_sheet"] = "EncyclopeDIA"

    # ---- the cholangiocarcinoma cohorts
    out["cholangiocarcinoma"] = ic

    # ---- blocks quoted from the first analysis that are not in the tables above
    first = next(iter(pc["datasets"]))
    out["controls"] = {
        name: {"erythrocyte_mean_z": v.get("erythrocyte_mean_z"),
               "immunoglobulin_mean_z": v.get("immunoglobulin_mean_z"),
               "off_panel_liver_plasma_mean_z":
                   v.get("off_panel_plasma_mean_z")}
        for name, v in pc["datasets"].items() if v["level"] == "proteome"}
    out["jiang_normalization_sensitivity"] = {
        norm: {ms: {"paired_delta": blk.get("paired_delta"),
                    "wilcoxon_p": blk.get("wilcoxon_p")}
               for ms, blk in v.items() if isinstance(blk, dict)}
        for norm, v in pc["jiang_normalization_sensitivity"].items()}
    out["matrix_scale"] = {name: v.get("matrix_scale")
                           for name, v in pc["datasets"].items()}
    out["abundance"] = {name: v.get("abundance")
                        for name, v in pc["datasets"].items()}
    out["class_difference_secreted_minus_metabolic"] = {
        name: v.get("class_difference_secreted_minus_metabolic")
        for name, v in pc["datasets"].items()}
    with open(f"{R}/MANUSCRIPT_EXTRAS.json") as fh:
        out["provenance_and_robustness"] = json.load(fh)

    with open(f"{R}/MANUSCRIPT_NUMBERS.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print(f"{out['n_datasets']} datasets assembled\n")
    hdr = f"{'dataset':34s} {'level':13s} {'pairs':>5s} {'metab z':>8s} {'share':>6s} {'plasma z':>9s} {'share':>6s}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['dataset'][:34]:34s} {r['level']:13s} {r['n_pairs']:>5d} "
              f"{r['metabolic_mean_z']:>+8.3f} {r['metabolic_share']:>6.3f} "
              f"{r['secreted_mean_z']:>+9.3f} {r['secreted_share']:>6.3f}")
    print(f"\nwrote {R}/MANUSCRIPT_NUMBERS.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
