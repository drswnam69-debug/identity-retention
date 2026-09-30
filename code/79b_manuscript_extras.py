#!/usr/bin/env python3
"""79b_manuscript_extras.py

The manuscript quotes a few numbers that were computed while checking the
analysis rather than by 70, 77 or 78: matrix completeness, accession map rates,
the leave-one-out range, the technical-replicate agreement, the raw-scale
separation in the cholangiocarcinoma cohorts and the per-member pair counts.
Rather than leave them hand-entered, this script recomputes them and writes
results/MANUSCRIPT_EXTRAS.json so that 80_verify_MCP.py can trace them too.
"""
import csv, gzip, importlib.util, json, os, re
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
U = "/root/.claude/uploads/f26f35bd-7d55-5294-804e-c3a8249fe52a"
HOLD = "/home/claude/holdout"


def mod(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p70 = mod(f"{HERE}/70_panel_coherence.py", "p70")
m78 = mod(f"{HERE}/78_icca_validation.py", "m78")

MODULE = ["CPS1", "OTC", "ARG1", "TAT", "G6PC1", "PCK1"]
PLASMA = ["ALB", "TTR", "TF", "SERPINA1", "AHSG", "APOH", "FGA", "FGB", "FGG", "F2"]
TESTED = MODULE + PLASMA[:]


def main() -> int:
    out = {"what": "numbers quoted in the manuscript that are provenance or "
                   "robustness rather than primary results"}

    # ---- 1. Yi 2023: completeness per sheet and the primary-sheet analysis
    yi = {}
    for sheet in ("EncyclopeDIA", "Spectronaut", "DIANN"):
        m = pd.read_csv(f"{HOLD}/yi2023_{sheet}_proteins_A.tsv.gz", sep="\t", index_col=0)
        yi[sheet] = {"n_symbols": int(m.shape[0]),
                     "completeness": float(m.notna().mean().mean())}
    out["yi2023_matrix_completeness"] = yi
    # accession map rate on the primary Yi matrix, quoted in Experimental Procedures
    amap0 = {}
    with gzip.open(f"{HOLD}/uniprot_to_symbol.tsv.gz", "rt") as fh:
        r0 = csv.reader(fh, delimiter="\t"); next(r0)
        for a0, g0 in r0:
            amap0[a0] = g0
    import openpyxl as _ox
    wb0 = _ox.load_workbook(
        "/mnt/user-data/uploads/1-s2.0-S1535947623001159-mmc4.xlsx", read_only=True)
    it0 = wb0["DIANN_proteins"].iter_rows(min_col=1, max_col=1, values_only=True)
    next(it0)
    tot0 = hit0 = 0
    for row0 in it0:
        if row0[0] is None:
            continue
        tot0 += 1
        if any(a in amap0 for a in
               [x.strip() for x in str(row0[0]).split(";")]):
            hit0 += 1
    wb0.close()
    out["yi2023_accession_map_rate"] = {
        "protein_groups": tot0, "mapped": hit0, "rate": hit0 / tot0}

    hcc = set(json.load(open(f"{HOLD}/hcc_paired_patients.json")))
    m = pd.read_csv(f"{HOLD}/yi2023_EncyclopeDIA_proteins_A.tsv.gz",
                    sep="\t", index_col=0).groupby(level=0).max()
    tc = [f"T{p}" for p in sorted(hcc, key=int)]
    pc = [f"N{p}" for p in sorted(hcc, key=int)]

    d = p70.member_shifts(m, tc, pc, MODULE)
    loo = {}
    for drop in list(d):
        z = [v["z"] for g, v in d.items() if g != drop]
        loo[f"without_{drop}"] = {"n": len(z), "mean_z": float(np.mean(z)),
                                  "surviving_share": p70.coherence(z)["cancellation"]}
    out["yi2023_leave_one_out"] = {
        "per_member": loo,
        "mean_z_min": min(v["mean_z"] for v in loo.values()),
        "mean_z_max": max(v["mean_z"] for v in loo.values()),
        "share_min": min(v["surviving_share"] for v in loo.values()),
        "share_max": max(v["surviving_share"] for v in loo.values())}

    sub = m.loc[[g for g in MODULE if g in d]]
    z = sub.sub(sub[tc + pc].mean(axis=1), axis=0).div(
        sub[tc + pc].std(axis=1, ddof=1), axis=0).mean(axis=0)
    delta = np.array([z[t] - z[p] for t, p in zip(tc, pc)])
    out["yi2023_per_pair_module"] = {
        "n_pairs": int(len(delta)), "n_pairs_falling": int((delta < 0).sum()),
        "mean_paired_delta": float(delta.mean()),
        "wilcoxon_p": float(stats.wilcoxon(delta).pvalue)}

    # ---- 2. Yi 2023 technical replicates
    import openpyxl
    amap = {}
    with gzip.open(f"{HOLD}/uniprot_to_symbol.tsv.gz", "rt") as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for a, g in r:
            amap[a] = {"G6PC": "G6PC1"}.get(g, g)
    wb = openpyxl.load_workbook(
        f"{U}/95db7fac-41467_2026_70817_MOESM4_ESM.xlsx", read_only=True) \
        if False else openpyxl.load_workbook(
        "/mnt/user-data/uploads/1-s2.0-S1535947623001159-mmc4.xlsx", read_only=True)
    ws = wb["EncyclopeDIA_proteins"]
    it = ws.iter_rows(values_only=True); hdr = list(next(it))
    cols = hdr[4:]
    pat = re.compile(r"_DIA_([NT])_(\d+)([AB])$")
    idx = {}
    for i, c in enumerate(cols):
        mm = pat.search(str(c))
        if mm:
            idx[(mm.group(1), mm.group(2), mm.group(3))] = i
    vals = {}
    for row in it:
        if row[0] is None:
            continue
        g = next((amap[a] for a in
                  [re.sub(r"-\d+$", "", x.strip()) for x in re.split(r"[;,]", str(row[0]))]
                  if a in amap), None)
        if g in TESTED:
            v = np.array([row[4 + i] if isinstance(row[4 + i], (int, float)) else np.nan
                          for i in range(len(cols))], float)
            vals[g] = np.fmax(vals[g], v) if g in vals else v
    wb.close()
    B = sorted({k[1] for k in idx if k[2] == "B"} & hcc, key=int)
    pairsAB = [(vals[g][idx[(arm, p, "A")]], vals[g][idx[(arm, p, "B")]])
               for g in vals for arm in ("N", "T") for p in B
               if (arm, p, "A") in idx and (arm, p, "B") in idx
               and np.isfinite(vals[g][idx[(arm, p, "A")]])
               and np.isfinite(vals[g][idx[(arm, p, "B")]])]
    a = np.array([x[0] for x in pairsAB]); b = np.array([x[1] for x in pairsAB])
    out["yi2023_technical_replicates"] = {
        "n_hcc_patients_with_a_replicate": len(B),
        "n_proteins_examined": len(vals),
        "n_value_pairs": int(len(a)),
        "pearson_r": float(np.corrcoef(a, b)[0, 1]),
        "median_absolute_difference_log2": float(np.median(np.abs(a - b)))}

    # ---- 3. the secreted annotation's coverage of the comparison sets
    comp = {}
    with gzip.open(f"{ROOT}/data/Human.GRCh38.p13.annot.tsv.gz", "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            comp[row["Symbol"]] = {i for i in
                                   (row.get("GOComponentID") or "").split("///")
                                   if i.startswith("GO:")}
    sets = json.load(open(f"{ROOT}/genesets/eligible.json"))
    universe = set()
    for v in sets.values():
        universe |= set(v)
    inann = [g for g in universe if g in comp]
    out["annotation_coverage"] = {
        "genes_in_the_eligible_sets": len(universe),
        "present_in_the_annotation": len(inann),
        "share_present": len(inann) / len(universe),
        "with_at_least_one_component_term":
            sum(1 for g in inann if comp[g]),
        "share_with_a_component_term":
            sum(1 for g in inann if comp[g]) / len(universe)}

    # ---- 4. cholangiocarcinoma: pair counts and the raw scale
    amap2 = m78.acc_map(f"{HOLD}/uniprot_to_symbol.tsv.gz", paper_tables=[
        (f"{U}/f0fd543f-41467_2026_70817_MOESM5_ESM.xlsx", "Limmaresult_Tumor_TANM"),
        (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx", "Coxph_TTR"),
        (f"{U}/26629c2b-41467_2026_70817_MOESM6_ESM.xlsx", "Coxph_OS")])
    ic = {}
    for label, mx, sheet, sx, tl, nl in [
        ("MSKCC-ICC", f"{U}/95db7fac-41467_2026_70817_MOESM4_ESM.xlsx",
         "TimsTOF_ICC_ExprMat_log2_median",
         f"{U}/a7ddcca0-41467_2026_70817_MOESM3_ESM.xlsx", "Tumor", "normal"),
        ("UKF-ICC", f"{U}/8ae8e955-41467_2026_70817_MOESM7_ESM.xlsx",
         "UKF_ICC_ExprMat_log2_median",
         f"{U}/f53c7329-41467_2026_70817_MOESM11_ESM.xlsx", "Tumor", "TANM")]:
        mat, t, p, prov = m78.load_cohort(mx, sheet, sx, tl, nl, amap2)
        rec = {"provenance": prov, "completeness": float(mat.notna().mean().mean()),
               "members": {}}
        for g in MODULE:
            if g not in mat.index:
                continue
            row = mat.loc[g]
            if isinstance(row, pd.DataFrame):
                row = row.max()
            tv = np.array([row[c] for c in t], float)
            av = np.array([row[c] for c in p], float)
            ok = np.isfinite(tv) & np.isfinite(av)
            if ok.sum() < p70.MIN_PAIRS:
                continue
            rec["members"][g] = {
                "n_pairs": int(ok.sum()),
                "tumor_median": float(np.median(tv[ok])),
                "adjacent_median": float(np.median(av[ok])),
                "mean_raw_difference_log2": float(np.mean(tv[ok] - av[ok])),
                "n_pairs_falling": int((tv[ok] < av[ok]).sum())}
        for g in ("ALB", "FGA"):
            if g in mat.index:
                row = mat.loc[g]
                if isinstance(row, pd.DataFrame):
                    row = row.max()
                tv = np.array([row[c] for c in t], float)
                av = np.array([row[c] for c in p], float)
                ok = np.isfinite(tv) & np.isfinite(av)
                rec["members"][g] = {"n_pairs": int(ok.sum()),
                                     "n_pairs_rising": int((tv[ok] > av[ok]).sum())}
        raws = [v["mean_raw_difference_log2"] for g, v in rec["members"].items()
                if g in MODULE]
        rec["module_raw_difference_range_log2"] = [min(raws), max(raws)] if raws else None
        # exact composition of the deposited matrix, quoted in the paper
        import pandas as _pd
        smp = _pd.read_excel(sx)
        rec["sample_table"] = {
            "n_samples": int(len(smp)),
            "n_patients": int(smp["Patient_No"].nunique()),
            "n_tumor": int((smp["Group"] == tl).sum()),
            "n_adjacent": int((smp["Group"] == nl).sum()),
            "n_patients_with_both_arms": prov["patients_with_both_arms"]}
        ic[label] = rec
    out["cholangiocarcinoma_extras"] = ic
    out["cholangiocarcinoma_module_raw_difference_range_log2"] = [
        min(v["module_raw_difference_range_log2"][0] for v in ic.values()),
        max(v["module_raw_difference_range_log2"][1] for v in ic.values())]

    with open(f"{ROOT}/results/MANUSCRIPT_EXTRAS.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("cholangiocarcinoma_extras",)}, indent=1)[:2600])
    print("\ncholangiocarcinoma module raw differences (log2):",
          out["cholangiocarcinoma_module_raw_difference_range_log2"])
    for L, v in ic.items():
        print(f"  {L}: completeness {v['completeness']:.3f}; "
              + ", ".join(f"{g} n={x['n_pairs']}" for g, x in v["members"].items()))
    print(f"\nwrote {ROOT}/results/MANUSCRIPT_EXTRAS.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
