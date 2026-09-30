#!/usr/bin/env python3
"""77_holdout_validation.py

The pre-registered test of the metabolic module on a held-out paired liver
proteome. Written on 2026-09-29, BEFORE the validation matrix was obtained,
against the predictions recorded the same day. The only part added after the
file arrives is `load_holdout`, the adapter that turns that particular
supplementary table into a gene-by-sample matrix and a pairing; every decision
rule below is fixed now.

Usage:
    python code/77_holdout_validation.py <matrix.xlsx> <samples.xlsx>

Definitions are imported from 70_panel_coherence.py so that member_shifts,
coherence and matrix_diagnostic mean exactly what they mean in the main
analysis. Writes results/HOLDOUT_VALIDATION.json.
"""
import importlib.util, json, os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location(
    "panel_coherence", os.path.join(HERE, "70_panel_coherence.py"))
p70 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p70)

MODULE = ["CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1"]          # H1, H3
PLASMA = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH",
          "FGA", "FGB", "FGG", "F2"]                              # H2
PANEL = list(p70.D1)                                              # H4
# orthogonal sample-identity gate: nothing here is in the panel
PROLIFERATION = ["MKI67", "PCNA", "TOP2A", "RRM2", "TYMS", "MCM2",
                 "CCNB1", "AURKA"]
SECRETED_TERM = "GO:0005615"

GATE_MIN_MEMBERS = 3        # the gate needs at least this many proliferation proteins
H1_MIN_SHARE = 0.9          # surviving share the module must keep
H3_MAX_PERCENTILE = 10.0    # the module must sit below this percentile


def secreted_set():
    import csv, gzip
    path = os.path.join(ROOT, "data", "Human.GRCh38.p13.annot.tsv.gz")
    out = set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ids = (row.get("GOComponentID") or "").split("///")
            if SECRETED_TERM in ids:
                out.add(row["Symbol"])
    return out


def load_holdout(matrix_path, samples_path):
    """ADAPTER for Yi 2023 (Mol Cell Proteomics), written 2026-09-30 after the
    supplementary tables arrived. Everything below this function was fixed the
    day before, without the data.

    matrix_path: a gene-symbol-indexed TSV built from one protein sheet of
                 Supplemental Table S2 (mmc4.xlsx), columns named T<id>/N<id>
                 for the primary (Tech_rep "A") run of each patient arm.
    samples_path: Supplemental Table S1 (mmc3.xlsx), which carries
                 Patient_DIA (patient id -> tumor type) and Sample_DIA
                 (file -> tissue type, patient, technical replicate).

    Returns the matrix restricted to the HCC patients that have both arms.
    """
    import openpyxl
    wb = openpyxl.load_workbook(samples_path, read_only=True)
    it = wb["Patient_DIA"].iter_rows(values_only=True); next(it)
    ptype = {str(r[0]).strip(): str(r[1]).strip() for r in it if r[0] is not None}
    it = wb["Sample_DIA"].iter_rows(values_only=True); next(it)
    arms = {}
    for r in it:
        if r[0] is None or str(r[4]) != "A":
            continue
        pid, tis = str(r[2]).strip(), str(r[1]).strip()
        if ptype.get(pid) != "HCC":          # 41 HCC of 56; CCA and MCA excluded
            continue
        if tis not in ("N", "T"):
            continue
        arms.setdefault(pid, set()).add(tis)
    wb.close()
    paired = sorted([p for p, a in arms.items() if a == {"N", "T"}], key=int)

    mat = pd.read_csv(matrix_path, sep="\t", index_col=0)
    mat = mat.groupby(level=0).max()
    tc = [f"T{p}" for p in paired if f"T{p}" in mat.columns and f"N{p}" in mat.columns]
    pcs = [f"N{p}" for p in [c[1:] for c in tc]]
    print(f"  {len(ptype)} patients in the table "
          f"({sum(1 for v in ptype.values() if v == 'HCC')} HCC); "
          f"{len(paired)} HCC patients with both arms; {len(tc)} pairs usable")
    return mat[tc + pcs], tc, pcs


def report_premises(mat, tc, pc):
    """Established before any contrast is computed, per the standing rule."""
    diag = p70.matrix_diagnostic(mat, list(tc) + list(pc))
    miss = p70.missingness(mat, tc, pc)
    prol = p70.member_shifts(mat, tc, pc, PROLIFERATION)
    zs = [v["z"] for v in prol.values()]
    gate = {
        "proliferation_members_quantified": {g: round(v["z"], 4)
                                             for g, v in prol.items()},
        "n_quantified": len(zs),
        "mean_z": float(np.mean(zs)) if zs else None,
        "required": f">= {GATE_MIN_MEMBERS} members and mean z > 0",
        "passes": bool(len(zs) >= GATE_MIN_MEMBERS and zs and np.mean(zs) > 0),
    }
    return {"n_pairs": len(tc), "n_genes": int(mat.shape[0]),
            "matrix_diagnostic": diag, "missingness": miss,
            "sample_identity_gate": gate}


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    mat, tc, pc = load_holdout(sys.argv[1], sys.argv[2])

    pre = report_premises(mat, tc, pc)
    print(json.dumps(pre["sample_identity_gate"], indent=2))
    out = {"what": "pre-registered held-out validation, predictions recorded "
                   "2026-09-29 before the data were obtained",
           "premises": pre}
    if not pre["sample_identity_gate"]["passes"]:
        out["verdict"] = ("NOT TESTABLE: the orthogonal sample-identity gate "
                          "failed, so the tumor and adjacent labels or the "
                          "pairing cannot be relied on. No contrast computed.")
        print("\n" + out["verdict"])
        with open(f"{ROOT}/results/HOLDOUT_VALIDATION.json", "w") as fh:
            json.dump(out, fh, indent=1, sort_keys=True)
        return 1

    secreted = secreted_set()
    allshift = p70.member_shifts(mat, tc, pc, list(mat.index))

    mod = p70.member_shifts(mat, tc, pc, MODULE)
    mz = [v["z"] for v in mod.values()]
    h1 = {"members_quantified": {g: round(v["z"], 4) for g, v in mod.items()},
          "n": len(mz), "mean_z": float(np.mean(mz)) if mz else None,
          "coherence": p70.coherence(mz) if mz else None,
          "prediction": f"mean z < 0 and surviving share >= {H1_MIN_SHARE}",
          "passes": bool(mz and np.mean(mz) < 0
                         and p70.coherence(mz)["cancellation"] >= H1_MIN_SHARE)}

    pla = p70.member_shifts(mat, tc, pc, PLASMA)
    pz = [v["z"] for v in pla.values()]
    h2 = {"members_quantified": {g: round(v["z"], 4) for g, v in pla.items()},
          "n": len(pz), "mean_z": float(np.mean(pz)) if pz else None,
          "coherence": p70.coherence(pz) if pz else None,
          "prediction": "mean z > 0",
          "passes": bool(pz and np.mean(pz) > 0)}

    nonsec = np.array([v["z"] for g, v in allshift.items() if g not in secreted])
    h3 = {"n_non_secreted_proteins": int(nonsec.size),
          "module_mean_z": h1["mean_z"],
          "percentile_of_module_mean_among_non_secreted":
              float(100.0 * np.mean(nonsec <= h1["mean_z"]))
              if nonsec.size and h1["mean_z"] is not None else None,
          "median_percentile_of_module_members": float(np.median(
              [100.0 * np.mean(nonsec <= v["z"]) for v in mod.values()]))
              if nonsec.size and mod else None,
          "prediction": f"module mean below the {H3_MAX_PERCENTILE:.0f}th percentile",
          }
    h3["passes"] = bool(h3["percentile_of_module_mean_among_non_secreted"]
                        is not None
                        and h3["percentile_of_module_mean_among_non_secreted"]
                        < H3_MAX_PERCENTILE)

    scored = p70.member_shifts(mat, tc, pc, PANEL)
    h4 = {"composite": p70.composite_sample_test(mat, tc, pc, PANEL,
                                                 restrict=list(scored)),
          "prediction": "Wilcoxon P not significant, or sign opposite to the "
                        "transcriptome sign (which is negative)"}

    h5 = {"applicable": bool(pre["matrix_diagnostic"].get("on_an_abundance_scale"))}
    if h5["applicable"]:
        ab = mat[list(tc) + list(pc)].mean(axis=1).reindex(list(allshift)).dropna()
        pctile = ab.rank(pct=True) * 100.0
        zser = pd.Series({g: allshift[g]["z"] for g in ab.index})
        others = [g for g in ab.index if g not in secreted and g not in MODULE]
        op = np.array([pctile[g] for g in others]); oz = np.array([zser[g] for g in others])
        order = np.argsort(op); op, oz = op[order], oz[order]
        cum = np.concatenate([[0.0], np.cumsum(oz)])
        rows = {}
        for g in [x for x in MODULE if x in ab.index]:
            lo = np.searchsorted(op, pctile[g] - 1.5, "left")
            hi = np.searchsorted(op, pctile[g] + 1.5, "right")
            if hi - lo >= 10:
                exp = (cum[hi] - cum[lo]) / (hi - lo)
                rows[g] = {"abundance_percentile": float(pctile[g]),
                           "observed_z": float(zser[g]),
                           "expected_from_abundance_matched_neighbors": float(exp),
                           "n_neighbors": int(hi - lo),
                           "share_explained": float(exp / zser[g])
                           if zser[g] else None}
        h5["members"] = rows
        if rows:
            # CORRECTED 2026-09-30. The first version averaged the per-member
            # ratios, which is wrong: a member whose observed shift is near
            # zero (TAT) returns a ratio of several hundred percent and drags
            # the average with it. The quantity that answers the question is
            # the ratio of the means: how much of the module's mean fall the
            # abundance-matched neighborhoods reproduce.
            obs = float(np.mean([r["observed_z"] for r in rows.values()]))
            exp = float(np.mean([r["expected_from_abundance_matched_neighbors"]
                                 for r in rows.values()]))
            h5["module_mean_observed"] = obs
            h5["module_mean_expected_from_abundance"] = exp
            h5["share_of_module_fall_explained_by_abundance"] = (
                float(exp / obs) if obs else None)
            h5["superseded_mean_of_per_member_ratios"] = float(
                np.mean([r["share_explained"] for r in rows.values()
                         if r["share_explained"] is not None]))
    else:
        h5["skipped"] = "matrix is not on an abundance scale"

    out.update({"H1_metabolic_module": h1, "H2_plasma_proteins": h2,
                "H3_extremity": h3, "H4_composite": h4, "H5_abundance": h5})
    out["verdict"] = {
        "module_claim_survives": bool(h1["passes"] and h3["passes"]),
        "two_opposing_axes_framing_survives": bool(h2["passes"]),
        "decision_rule": "H1 and H3 must both pass for the module claim; "
                         "H2 decides the framing; if H1 or H3 fails, 방향2B "
                         "is closed (recorded 2026-09-29)"}
    print("\n" + json.dumps(out["verdict"], indent=2))
    with open(f"{ROOT}/results/HOLDOUT_VALIDATION.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"wrote {ROOT}/results/HOLDOUT_VALIDATION.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
