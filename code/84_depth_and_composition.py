#!/usr/bin/env python3
"""84_depth_and_composition.py

Two alternative explanations that have not been addressed, and a decision rule
written before the numbers are seen.

THE TWO WORRIES
1. DEPTH. The tumor arm is quantified more deeply than the adjacent arm in every
   proteome with real missingness. With sample-median normalization and
   left-censoring, that asymmetry alone can raise anything near the detection
   limit and can bias any paired contrast computed on complete cases only.
2. COMPOSITION. The module members are among the most abundant proteins in liver
   parenchyma and they collapse in tumor. In a median-normalized matrix, losing
   that much protein mass mechanically raises everything that does not collapse.
   The plasma class rising could be the arithmetic shadow of the enzymes falling.

THE TESTS
A. Detection asymmetry per panel member: detections in each arm.
B. Complete-case analysis: restrict to proteins detected in BOTH arms of EVERY
   pair, then recompute the class means, the module and the module's percentile
   against a background restricted the same way. No member is then scored on a
   different set of patients from any other, and no protein enters through
   differential dropout.
C. Compositional accounting: the share of linear signal carried by the module in
   each arm, and the global log2 shift that renormalization must produce once
   that share is lost. Compare that shift with the observed rise of the plasma
   class.

DECISION RULE, fixed now.
  D1  If under B the module's mean shift loses its sign, or its percentile rises
      above the 25th, in the majority of the proteomes where B can be run, the
      module claim is driven by differential dropout and the work is closed.
  D2  If under C the implied renormalization shift accounts for more than half
      of the plasma class's observed rise in a cohort, no claim about plasma
      proteins rising may be made from that cohort.
  D3  The module claim survives only if it passes D1 and its fall is larger than
      the implied compositional shift by at least a factor of three.

Writes results/DEPTH_AND_COMPOSITION.json.
"""
import csv, gzip, importlib.util, json, os, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def _load(name, fn):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


p70 = _load("p70", "70_panel_coherence.py")
m78 = _load("m78", "78_icca_validation.py")

MODULE = ["CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1"]
PLASMA = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2"]
SECRETED_TERM = "GO:0005615"
U = "/root/.claude/uploads/f26f35bd-7d55-5294-804e-c3a8249fe52a"
HOLD = "/home/claude/holdout"


def secreted_set():
    out = set()
    with gzip.open(f"{ROOT}/data/Human.GRCh38.p13.annot.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if SECRETED_TERM in (row.get("GOComponentID") or "").split("///"):
                out.add(row["Symbol"])
    return out


def cohorts():
    out = {}
    mat, tc, pc = p70.load_jiang("raw")
    out["Jiang 2019"] = (mat, list(tc), list(pc), True)

    hcc = set(json.load(open(f"{HOLD}/hcc_paired_patients.json")))
    m = pd.read_csv(f"{HOLD}/yi2023_EncyclopeDIA_proteins_A.tsv.gz",
                    sep="\t", index_col=0).groupby(level=0).max()
    ids = sorted(hcc, key=int)
    out["Yi 2023"] = (m, [f"T{p}" for p in ids], [f"N{p}" for p in ids], True)

    amap = m78.acc_map(f"{HOLD}/uniprot_to_symbol.tsv.gz", paper_tables=[
        (f"{U}/f0fd543f-41467_2026_70817_MOESM5_ESM.xlsx", "Limmaresult_Tumor_TANM"),
        (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx", "Coxph_TTR"),
        (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx", "Coxph_OS")])
    for label, mx, sheet, sx, tl, nl in [
        ("MSKCC-ICC", f"{U}/95db7fac-41467_2026_70817_MOESM4_ESM.xlsx",
         "TimsTOF_ICC_ExprMat_log2_median",
         f"{U}/a7ddcca0-41467_2026_70817_MOESM3_ESM.xlsx", "Tumor", "normal"),
        ("UKF-ICC", f"{U}/8ae8e955-41467_2026_70817_MOESM7_ESM.xlsx",
         "UKF_ICC_ExprMat_log2_median",
         f"{U}/f53c7329-41467_2026_70817_MOESM11_ESM.xlsx", "Tumor", "TANM")]:
        mm, t, p, _ = m78.load_cohort(mx, sheet, sx, tl, nl, amap)
        out[label] = (mm, list(t), list(p), False)
    return out


def main() -> int:
    sec = secreted_set()
    res = {"what": "depth and composition checks; the decision rule was written "
                   "into this script before the numbers were seen",
           "cohorts": {}}

    for label, (mat, tc, pc, ref_median) in cohorts().items():
        T = mat[tc].to_numpy(float); A = mat[pc].to_numpy(float)
        okT = np.isfinite(T); okA = np.isfinite(A)
        rec = {"n_pairs": len(tc), "n_proteins": int(mat.shape[0]),
               "detections_per_tumor_sample": float(okT.sum(axis=0).mean()),
               "detections_per_adjacent_sample": float(okA.sum(axis=0).mean())}
        rec["tumor_deeper_by"] = float(rec["detections_per_tumor_sample"] /
                                       rec["detections_per_adjacent_sample"] - 1)

        # ---- A. detection asymmetry for the genes under test
        det = {}
        idx = list(mat.index)
        for g in MODULE + PLASMA:
            if g not in mat.index:
                continue
            i = idx.index(g)
            det[g] = {"detected_tumor": int(okT[i].sum()),
                      "detected_adjacent": int(okA[i].sum()),
                      "detected_both": int((okT[i] & okA[i]).sum())}
        rec["A_detection"] = det

        # ---- B. complete cases only
        complete = okT.all(axis=1) & okA.all(axis=1)
        genes_c = [g for g, k in zip(mat.index, complete) if k]
        rec["B_complete_case"] = {"n_proteins_complete": int(complete.sum()),
                                  "share_complete": float(complete.mean())}
        if complete.sum() >= 200:
            sub = mat.loc[genes_c]
            d = p70.member_shifts(sub, tc, pc, list(sub.index))
            allz = np.array([v["z"] for v in d.values()])
            nonsec = np.array([v["z"] for g, v in d.items() if g not in sec])
            mod = [d[g]["z"] for g in MODULE if g in d]
            pla = [d[g]["z"] for g in PLASMA if g in d]
            rec["B_complete_case"].update({
                "module_n": len(mod),
                "module_members": {g: round(d[g]["z"], 4) for g in MODULE if g in d},
                "module_mean_z": float(np.mean(mod)) if mod else None,
                "module_surviving_share": p70.coherence(mod)["cancellation"] if mod else None,
                "module_percentile_among_non_secreted":
                    float(100 * np.mean(nonsec <= np.mean(mod))) if mod and nonsec.size else None,
                "plasma_n": len(pla),
                "plasma_mean_z": float(np.mean(pla)) if pla else None,
                "plasma_members": {g: round(d[g]["z"], 4) for g in PLASMA if g in d},
                "background_mean_z": float(allz.mean())})

        # ---- C. compositional accounting
        if ref_median:
            # Jiang and Yi are log2 intensities, not centered; center on the
            # per-sample median so that the mass share is comparable
            med = np.nanmedian(np.where(okT, T, np.nan), axis=0)
            Tc = T - med
            med2 = np.nanmedian(np.where(okA, A, np.nan), axis=0)
            Ac = A - med2
        else:
            Tc, Ac = T, A           # already sample-median centered
        linT = np.where(np.isfinite(Tc), 2.0 ** Tc, 0.0)
        linA = np.where(np.isfinite(Ac), 2.0 ** Ac, 0.0)
        totT, totA = linT.sum(axis=0), linA.sum(axis=0)

        def share(genes):
            rows = [idx.index(g) for g in genes if g in mat.index]
            if not rows:
                return None, None
            return (float(np.mean(linT[rows].sum(axis=0) / totT)),
                    float(np.mean(linA[rows].sum(axis=0) / totA)))

        mT, mA = share(MODULE)
        pT, pA = share(PLASMA)
        c = {"module_signal_share_tumor": mT, "module_signal_share_adjacent": mA,
             "plasma_signal_share_tumor": pT, "plasma_signal_share_adjacent": pA}
        if mT is not None:
            lost = mA - mT                       # mass share the module gives up
            # renormalizing a matrix that loses `lost` of its mass lifts every
            # surviving protein by this much on a log2 scale
            c["module_mass_share_lost_in_tumor"] = float(lost)
            c["implied_global_log2_lift"] = float(-np.log2(1.0 - lost)) if lost < 1 else None
        res_plasma = rec.get("B_complete_case", {}).get("plasma_mean_z")
        c["observed_plasma_rise_complete_case"] = res_plasma
        rec["C_composition"] = c

        res["cohorts"][label] = rec
        print(f"\n=== {label}  ({len(tc)} pairs, {mat.shape[0]} proteins)")
        print(f"  depth: tumor {rec['detections_per_tumor_sample']:.0f} vs adjacent "
              f"{rec['detections_per_adjacent_sample']:.0f} "
              f"({100*rec['tumor_deeper_by']:+.1f}%)")
        b = rec["B_complete_case"]
        print(f"  complete-case proteins: {b['n_proteins_complete']} "
              f"({100*b['share_complete']:.1f}%)")
        if "module_mean_z" in b and b["module_mean_z"] is not None:
            print(f"    module  n={b['module_n']} mean_z={b['module_mean_z']:+.3f} "
                  f"share={b['module_surviving_share']:.3f} "
                  f"percentile={b['module_percentile_among_non_secreted']:.2f}")
            print(f"      {b['module_members']}")
            print(f"    plasma  n={b['plasma_n']} mean_z={b['plasma_mean_z']:+.3f}")
            print(f"    background mean over complete proteins: {b['background_mean_z']:+.3f}")
        if c.get("implied_global_log2_lift") is not None:
            print(f"  module signal share: adjacent {mA:.4f} -> tumor {mT:.4f}; "
                  f"mass lost {c['module_mass_share_lost_in_tumor']:+.4f}")
            print(f"  implied global lift from renormalization: "
                  f"{c['implied_global_log2_lift']:+.4f} log2")

    with open(f"{ROOT}/results/DEPTH_AND_COMPOSITION.json", "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    print(f"\nwrote {ROOT}/results/DEPTH_AND_COMPOSITION.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
