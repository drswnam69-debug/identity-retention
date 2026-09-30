#!/usr/bin/env python3
"""78_icca_validation.py

The second pre-registered test, on two paired intrahepatic cholangiocarcinoma
proteomes (Nat Commun 2026, doi 10.1038/s41467-026-70817-1). Predictions G1-G5
were recorded on 2026-09-30 before any of these files were downloaded.

The point of this cohort is that the tumor is biliary, not hepatocellular, so
the plasma-protein prediction runs the OTHER WAY from the HCC cohorts. If the
rise of plasma proteins in HCC proteomes reflects tumor cells that still
secrete, it should not happen here.

Inputs (Supplementary Data of that article):
  MOESM3   MSKCC sample table   (Sample, Patient_No, Group in {Tumor, normal})
  MOESM4   MSKCC protein matrix (log2, sample-median normalized)
  MOESM11  UKF sample table     (Sample, Patient_No, Group in {Tumor, TANM})
  MOESM7   UKF protein matrix

Writes results/ICCA_VALIDATION.json.
"""
import csv, gzip, importlib.util, json, os, re, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location(
    "panel_coherence", os.path.join(HERE, "70_panel_coherence.py"))
p70 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p70)

MODULE = ["CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1"]
PLASMA = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH",
          "FGA", "FGB", "FGG", "F2"]
PROLIFERATION = ["MKI67", "PCNA", "TOP2A", "RRM2", "TYMS", "MCM2",
                 "CCNB1", "AURKA"]
# DEVIATION, declared 2026-09-30 before G1-G5 were computed. The gate above
# failed on coverage in both cohorts: only PCNA and MCM2 are quantified in 10+
# pairs, because the rest are low-abundance nuclear proteins that tissue DIA
# does not reach. The replacement adds biliary epithelial markers, whose
# direction in an intrahepatic cholangiocarcinoma against adjacent liver is
# fixed by the biology without looking at the data, and none of which belongs
# to any class of the panel under test. Results obtained under this gate are
# EXPLORATORY, not confirmatory; the pre-registered verdict is reported too.
GATE_AMENDED = ["KRT19", "KRT7", "EPCAM", "S100P", "PCNA", "MCM2"]
SECRETED_TERM = "GO:0005615"
GATE_MIN = 3
G1_MIN_SHARE = 0.9
G5_MAX_PCT = 10.0
# for G2 and G4, the values already in hand from the HCC proteomes
HCC_MODULE_MEAN = {"Gao 2019": -1.419, "Jiang 2019": -0.899, "Yi 2023": -0.715}
HCC_GAP = {"Gao 2019": 1.827, "Jiang 2019": 1.413, "Yi 2023": 0.857}


def secreted_set():
    path = os.path.join(ROOT, "data", "Human.GRCh38.p13.annot.tsv.gz")
    out = set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if SECRETED_TERM in (row.get("GOComponentID") or "").split("///"):
                out.add(row["Symbol"])
    return out


def acc_map(path, paper_tables=()):
    """UniProt accession -> gene symbol.

    Two sources, the paper's own first. `paper_tables` are (xlsx, sheet) pairs
    from this article's Supplementary Data that carry a `Protein` column and a
    `Gene.names` column; ONLY those two columns are read, so no result column of
    theirs enters this analysis. The fallback is the map built from the Jiang
    MaxQuant proteinGroups.txt already deposited in this archive.
    """
    m = {}
    for xlsx, sheet in paper_tables:
        try:
            df = pd.read_excel(xlsx, sheet_name=sheet,
                               usecols=["Protein", "Gene.names"])
        except Exception:
            continue
        for acc, gene in zip(df["Protein"], df["Gene.names"]):
            if not isinstance(gene, str) or not isinstance(acc, str):
                continue
            g = gene.split(";")[0].strip().split()[0] if gene.strip() else ""
            if not g:
                continue
            g = {"G6PC": "G6PC1"}.get(g, g)
            for a in re.split(r"[;,]", acc):
                a = re.sub(r"-\d+$", "", a.strip())
                if a:
                    m.setdefault(a, g)
    n_paper = len(m)
    with gzip.open(path, "rt") as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for a, g in r:
            m.setdefault(a, {"G6PC": "G6PC1"}.get(g, g))
    print(f"accession map: {n_paper} from the article's own tables, "
          f"{len(m) - n_paper} added from the deposited Jiang proteinGroups, "
          f"{len(m)} total")
    return m


def to_symbol(protein, amap):
    for a in re.split(r"[;,]", str(protein)):
        a = re.sub(r"^SWISS-PROT:", "", a.strip())
        a = re.sub(r"-\d+$", "", a)
        if a in amap:
            return amap[a]
    return None


def load_cohort(matrix_xlsx, sheet, sample_xlsx, tumor_label, normal_label, amap):
    import openpyxl
    samples = pd.read_excel(sample_xlsx)
    grp = {str(r["Sample"]): (str(r["Group"]), str(r["Patient_No"]))
           for _, r in samples.iterrows()}
    wb = openpyxl.load_workbook(matrix_xlsx, read_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    cols = [str(c) for c in hdr[1:]]
    data, unmapped = {}, 0
    for row in it:
        if row[0] is None:
            continue
        g = to_symbol(row[0], amap)
        if g is None:
            unmapped += 1
            continue
        v = np.array([x if isinstance(x, (int, float)) else np.nan
                      for x in row[1:]], float)
        data[g] = np.fmax(data[g], v) if g in data else v
    wb.close()
    mat = pd.DataFrame.from_dict(data, orient="index", columns=cols)
    # one column per patient arm; no patient has more than one sample per arm
    tmap = {p: s for s, (gr, p) in grp.items() if gr == tumor_label and s in mat.columns}
    nmap = {p: s for s, (gr, p) in grp.items() if gr == normal_label and s in mat.columns}
    paired = sorted(set(tmap) & set(nmap))
    tc = [f"T_{p}" for p in paired]
    pc = [f"N_{p}" for p in paired]
    out = pd.DataFrame(index=mat.index)
    for p in paired:
        out[f"T_{p}"] = mat[tmap[p]]
        out[f"N_{p}"] = mat[nmap[p]]
    return out, tc, pc, {"proteins_in_file": len(data) + unmapped,
                         "mapped_to_symbol": len(data),
                         "map_rate": round(len(data) / (len(data) + unmapped), 4),
                         "samples_in_matrix": len(cols),
                         "patients_with_both_arms": len(paired)}


def score(mat, tc, pc, secreted, label):
    res = {"cohort": label}
    res["matrix_diagnostic"] = p70.matrix_diagnostic(mat, tc + pc)
    res["missingness"] = p70.missingness(mat, tc, pc)

    prol = p70.member_shifts(mat, tc, pc, PROLIFERATION)
    zs = [v["z"] for v in prol.values()]
    res["gate_preregistered"] = {
        "members": {g: round(v["z"], 4) for g, v in prol.items()},
        "n": len(zs), "mean_z": float(np.mean(zs)) if zs else None,
        "rule": f">= {GATE_MIN} quantified members and mean z > 0",
        "passes": bool(len(zs) >= GATE_MIN and zs and np.mean(zs) > 0),
        "why_it_failed": "only PCNA and MCM2 reach 10 pairs; the other six are "
                         "low-abundance nuclear proteins this assay does not "
                         "quantify. The two that are quantified are strongly "
                         "positive, so the count failed, not the direction."}
    am = p70.member_shifts(mat, tc, pc, GATE_AMENDED)
    az = [v["z"] for v in am.values()]
    res["gate_amended"] = {
        "members": {g: round(v["z"], 4) for g, v in am.items()},
        "n": len(az), "mean_z": float(np.mean(az)) if az else None,
        "passes": bool(len(az) >= GATE_MIN and az and np.mean(az) > 0)}
    res["status"] = ("EXPLORATORY: the pre-registered gate failed on coverage "
                     "and was amended after the fact. The confirmatory "
                     "validation is Yi 2023 alone.")
    if not res["gate_amended"]["passes"]:
        res["verdict"] = "NOT TESTABLE: the amended gate failed as well"
        return res

    allsh = p70.member_shifts(mat, tc, pc, list(mat.index))
    nonsec = np.array([v["z"] for g, v in allsh.items() if g not in secreted])

    mod = p70.member_shifts(mat, tc, pc, MODULE)
    mz = [v["z"] for v in mod.values()]
    pla = p70.member_shifts(mat, tc, pc, PLASMA)
    pz = [v["z"] for v in pla.values()]

    res["G1_module"] = {
        "members": {g: round(v["z"], 4) for g, v in mod.items()},
        "n": len(mz), "mean_z": float(np.mean(mz)) if mz else None,
        "coherence": p70.coherence(mz) if mz else None,
        "passes": bool(mz and np.mean(mz) < 0
                       and p70.coherence(mz)["cancellation"] >= G1_MIN_SHARE)}
    res["G2_deeper_than_hcc"] = {
        "module_mean_z": res["G1_module"]["mean_z"],
        "hcc_values": HCC_MODULE_MEAN,
        "deeper_than_all_three": bool(res["G1_module"]["mean_z"] is not None and
                                      res["G1_module"]["mean_z"]
                                      < min(HCC_MODULE_MEAN.values())),
        "deeper_than_the_median_hcc": bool(res["G1_module"]["mean_z"] is not None and
                                           res["G1_module"]["mean_z"] < -0.899)}
    res["G3_plasma"] = {
        "members": {g: round(v["z"], 4) for g, v in pla.items()},
        "n": len(pz), "mean_z": float(np.mean(pz)) if pz else None,
        "coherence": p70.coherence(pz) if pz else None,
        "prediction": "mean z < 0, the opposite of the HCC cohorts",
        "passes": bool(pz and np.mean(pz) < 0)}
    gap = (res["G3_plasma"]["mean_z"] - res["G1_module"]["mean_z"]
           if None not in (res["G3_plasma"]["mean_z"], res["G1_module"]["mean_z"])
           else None)
    res["G4_gap"] = {"secreted_minus_metabolic": gap, "hcc_values": HCC_GAP,
                     "smaller_than_all_hcc": bool(gap is not None
                                                  and gap < min(HCC_GAP.values()))}
    res["G5_extremity"] = {
        "n_non_secreted": int(nonsec.size),
        "percentile_of_module_mean": float(100.0 * np.mean(nonsec <= res["G1_module"]["mean_z"]))
        if nonsec.size and res["G1_module"]["mean_z"] is not None else None,
        "median_percentile_of_members": float(np.median(
            [100.0 * np.mean(nonsec <= v["z"]) for v in mod.values()]))
        if nonsec.size and mod else None}
    res["G5_extremity"]["passes"] = bool(
        res["G5_extremity"]["percentile_of_module_mean"] is not None
        and res["G5_extremity"]["percentile_of_module_mean"] < G5_MAX_PCT)
    return res


def main() -> int:
    U = "/root/.claude/uploads/f26f35bd-7d55-5294-804e-c3a8249fe52a"
    amap = acc_map("/home/claude/holdout/uniprot_to_symbol.tsv.gz",
                   paper_tables=[
                       (f"{U}/f0fd543f-41467_2026_70817_MOESM5_ESM.xlsx",
                        "Limmaresult_Tumor_TANM"),
                       (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx",
                        "Coxph_TTR"),
                       (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx",
                        "Coxph_OS")])
    secreted = secreted_set()
    cohorts = [
        ("MSKCC-ICC", f"{U}/95db7fac-41467_2026_70817_MOESM4_ESM.xlsx",
         "TimsTOF_ICC_ExprMat_log2_median",
         f"{U}/a7ddcca0-41467_2026_70817_MOESM3_ESM.xlsx", "Tumor", "normal"),
        ("UKF-ICC", f"{U}/8ae8e955-41467_2026_70817_MOESM7_ESM.xlsx",
         "UKF_ICC_ExprMat_log2_median",
         f"{U}/f53c7329-41467_2026_70817_MOESM11_ESM.xlsx", "Tumor", "TANM"),
    ]
    out = {"what": "second pre-registered test; predictions recorded 2026-09-30 "
                   "before the files were downloaded",
           "tumor_type": "intrahepatic cholangiocarcinoma (biliary, not hepatocellular)",
           "cohorts": {}}
    for label, mx, sheet, sx, tl, nl in cohorts:
        mat, tc, pc, prov = load_cohort(mx, sheet, sx, tl, nl, amap)
        r = score(mat, tc, pc, secreted, label)
        r["provenance"] = prov
        out["cohorts"][label] = r
        print(f"\n=== {label}  {prov}")
        gp, ga = r["gate_preregistered"], r["gate_amended"]
        print(f"  gate (pre-registered) n={gp['n']} mean_z={gp['mean_z']:+.3f} "
              f"-> {'PASS' if gp['passes'] else 'FAIL (count, not direction)'}")
        print(f"  gate (amended, exploratory) n={ga['n']} mean_z={ga['mean_z']:+.3f} "
              f"-> {'PASS' if ga['passes'] else 'FAIL'}")
        print(f"     {ga['members']}")
        if not ga["passes"]:
            continue
        print(f"  G1 module  n={r['G1_module']['n']} mean_z={r['G1_module']['mean_z']:+.3f} "
              f"share={r['G1_module']['coherence']['cancellation']:.3f} "
              f"-> {'PASS' if r['G1_module']['passes'] else 'FAIL'}")
        print(f"     {r['G1_module']['members']}")
        print(f"  G2 deeper than all three HCC: {r['G2_deeper_than_hcc']['deeper_than_all_three']}"
              f"  (than the median HCC: {r['G2_deeper_than_hcc']['deeper_than_the_median_hcc']})")
        print(f"  G3 plasma  n={r['G3_plasma']['n']} mean_z={r['G3_plasma']['mean_z']:+.3f} "
              f"share={r['G3_plasma']['coherence']['cancellation']:.3f} "
              f"-> {'PASS (falls, as predicted)' if r['G3_plasma']['passes'] else 'FAIL (rises)'}")
        print(f"     {r['G3_plasma']['members']}")
        print(f"  G4 gap={r['G4_gap']['secreted_minus_metabolic']:+.3f} "
              f"smaller than all HCC: {r['G4_gap']['smaller_than_all_hcc']}")
        print(f"  G5 module pct={r['G5_extremity']['percentile_of_module_mean']:.2f} "
              f"median member pct={r['G5_extremity']['median_percentile_of_members']:.2f} "
              f"-> {'PASS' if r['G5_extremity']['passes'] else 'FAIL'}")
    with open(f"{ROOT}/results/ICCA_VALIDATION.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"\nwrote {ROOT}/results/ICCA_VALIDATION.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
