#!/usr/bin/env python3
"""76_secretome_verify.py

75_secretome_transfer.py found that within published liver gene sets the
secreted members sit below the rest at transcript level and above the rest at
protein level. That test treats 114 MSigDB sets as if they were independent,
which they are not, and it does not rule out abundance or plasma exudation.

This script attacks the finding four ways:

  A  gene level, no gene sets at all: every quantified gene in each matrix,
     split by annotation. One test per dataset, so set overlap cannot matter.
  B  abundance matched, in the matrices that are on an abundance scale.
  C  liver-made secreted proteins against secreted proteins the liver does not
     make, which separates hepatocyte secretion from plasma exudation.
  D  how much the 114 sets overlap, and the set-level test repeated on a
     greedily chosen gene-disjoint subfamily.
  E  per gene, the agreement between its transcript-level and protein-level
     shift, computed separately for secreted and other genes.

Writes results/SECRETOME_VERIFY.json.
"""
import csv, gzip, importlib.util, json, os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location(
    "panel_coherence", os.path.join(HERE, "70_panel_coherence.py"))
p70 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p70)

ES = "GO:0005615"
# proteins in plasma that hepatocytes do not make (exudation reporters)
NOT_LIVER_MADE = ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHM", "IGKC",
                  "IGLC1", "IGLC2", "JCHAIN", "HBB", "HBA1", "HBD", "CA1",
                  "SLC4A1", "PRDX2"]


def components():
    path = os.path.join(ROOT, "data", "Human.GRCh38.p13.annot.tsv.gz")
    out = {}
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ids = (row.get("GOComponentID") or "").split("///")
            out[row["Symbol"]] = {i for i in ids if i.startswith("GO:")}
    return out


def all_gene_shifts(mat, tc, pc, min_pairs=10):
    """paired shift for every gene in the matrix, not just a set's members."""
    return p70.member_shifts(mat, tc, pc, list(mat.index))


def mw(a, b):
    if len(a) < 5 or len(b) < 5:
        return None
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    # rank-biserial correlation, a bounded effect size
    rb = 2.0 * u.statistic / (len(a) * len(b)) - 1.0
    return {"n_secreted": len(a), "n_other": len(b),
            "mean_secreted": float(np.mean(a)), "mean_other": float(np.mean(b)),
            "median_secreted": float(np.median(a)), "median_other": float(np.median(b)),
            "difference_of_means": float(np.mean(a) - np.mean(b)),
            "mannwhitney_p": float(u.pvalue), "rank_biserial": float(rb)}


def main() -> int:
    comp = components()
    secreted = {g for g, s in comp.items() if ES in s}

    datasets = {
        "GSE14520 (microarray)":
            ("transcriptome",) + p70.load_tx(
                f"{p70.IR_DATA}/GSE14520_symbols.tsv.gz",
                f"{p70.IR_RESULTS}/GSE14520/phenotype.tsv"),
        "GSE76427 (microarray)":
            ("transcriptome",) + p70.load_tx(
                f"{p70.IR_DATA}/GSE76427_symbols.tsv.gz",
                f"{p70.IR_RESULTS}/GSE76427/phenotype.tsv"),
        "TCGA-LIHC (RNA sequencing)":
            ("transcriptome",) + p70.load_tx(
                f"{p70.IR_DATA}/TCGA_LIHC_symbols.tsv.gz",
                f"{p70.IR_RESULTS}/TCGA_LIHC/phenotype.tsv", "tumor"),
        "Gao 2019 proteome":   ("proteome",) + p70.load_gao(),
        "Jiang 2019 proteome": ("proteome",) + p70.load_jiang("raw"),
    }

    out = {"what": "post hoc; not pre-registered",
           "secreted_term": ES, "A_gene_level": {}, "B_abundance_matched": {},
           "C_liver_made_vs_not": {}, "E_cross_level_agreement": {}}
    shifts = {}

    for name, (level, mat, tc, pc) in datasets.items():
        d = all_gene_shifts(mat, tc, pc)
        shifts[name] = d
        sec = [v["z"] for g, v in d.items() if g in secreted]
        oth = [v["z"] for g, v in d.items() if g not in secreted]
        r = mw(sec, oth) or {}
        r["level"] = level
        r["n_genes_scored"] = len(d)
        out["A_gene_level"][name] = r
        print(f"A  {name:32s} {level[:5]}  n={len(d):6d}  "
              f"secreted {r.get('mean_secreted', float('nan')):+.3f} vs other "
              f"{r.get('mean_other', float('nan')):+.3f}  "
              f"diff {r.get('difference_of_means', float('nan')):+.3f}  "
              f"P={r.get('mannwhitney_p', float('nan')):.2g}")

        # ---- B: abundance matched, only where the matrix is on an abundance scale
        diag = p70.matrix_diagnostic(mat, list(tc) + list(pc))
        if diag.get("on_an_abundance_scale"):
            ab = mat[list(tc) + list(pc)].mean(axis=1)
            ab = ab.reindex([g for g in d]).dropna()
            pct = ab.rank(pct=True) * 100.0
            zser = pd.Series({g: d[g]["z"] for g in ab.index})
            sec_idx = [g for g in ab.index if g in secreted]
            oth_idx = [g for g in ab.index if g not in secreted]
            # for every secreted gene, draw the non-secreted genes within 1.5
            # percentile points of it, and average their shift
            op = np.array([pct[h] for h in oth_idx]); oz = np.array([zser[h] for h in oth_idx])
            order = np.argsort(op); op = op[order]; oz = oz[order]
            cum = np.concatenate([[0.0], np.cumsum(oz)])
            matched = []
            for g in sec_idx:
                lo = np.searchsorted(op, pct[g] - 1.5, "left")
                hi = np.searchsorted(op, pct[g] + 1.5, "right")
                if hi - lo >= 10:
                    matched.append(zser[g] - (cum[hi] - cum[lo]) / (hi - lo))
            w = stats.wilcoxon(matched) if len(matched) > 5 else None
            out["B_abundance_matched"][name] = {
                "level": level,
                "median_secreted_abundance_percentile":
                    float(np.median([pct[g] for g in sec_idx])),
                "median_other_abundance_percentile":
                    float(np.median([pct[g] for g in oth_idx])),
                "n_secreted_matched": len(matched),
                "median_difference_from_abundance_matched_neighbors":
                    float(np.median(matched)) if matched else None,
                "mean_difference_from_abundance_matched_neighbors":
                    float(np.mean(matched)) if matched else None,
                "wilcoxon_p": float(w.pvalue) if w is not None else None}
            b = out["B_abundance_matched"][name]
            print(f"B  {name:32s} matched n={b['n_secreted_matched']:4d}  "
                  f"median diff {b['median_difference_from_abundance_matched_neighbors']:+.3f}  "
                  f"P={b['wilcoxon_p']:.2g}")
        else:
            out["B_abundance_matched"][name] = {
                "level": level, "skipped": "matrix is not on an abundance scale",
                "fraction_of_values_below_zero": diag.get("fraction_below_zero")}

        # ---- C: liver-made secreted against secreted the liver does not make
        nl = [v["z"] for g, v in d.items() if g in NOT_LIVER_MADE]
        lm = [v["z"] for g, v in d.items()
              if g in secreted and g not in NOT_LIVER_MADE]
        out["C_liver_made_vs_not"][name] = {
            "level": level,
            "n_not_liver_made": len(nl),
            "mean_not_liver_made": float(np.mean(nl)) if nl else None,
            "n_liver_made_secreted": len(lm),
            "mean_liver_made_secreted": float(np.mean(lm)) if lm else None,
            "members_not_liver_made": {g: round(d[g]["z"], 4)
                                       for g in NOT_LIVER_MADE if g in d}}

    # ---- D: how independent are the 114 sets, and a disjoint subfamily
    sets = json.load(open(f"{p70.IR_GENESETS}/eligible.json"))
    scored = [s for s in sets
              if all(sum(1 for g in sets[s] if g in shifts[n])
                     >= p70.MIN_MEMBERS for n in datasets)]
    js = []
    for i, a in enumerate(scored):
        for b in scored[i + 1:]:
            A, B = set(sets[a]), set(sets[b])
            js.append(len(A & B) / len(A | B))
    js = np.array(js)
    # greedy gene-disjoint subfamily, largest set first
    chosen, used = [], set()
    for s in sorted(scored, key=lambda s: -len(sets[s])):
        if not (set(sets[s]) & used):
            chosen.append(s); used |= set(sets[s])
    out["D_set_overlap"] = {
        "n_sets": len(scored),
        "median_pairwise_jaccard": float(np.median(js)),
        "mean_pairwise_jaccard": float(np.mean(js)),
        "share_of_pairs_sharing_any_gene": float(np.mean(js > 0)),
        "max_pairwise_jaccard": float(js.max()),
        "gene_disjoint_subfamily": chosen}
    print(f"D  114-set overlap: median Jaccard {np.median(js):.3f}, "
          f"{np.mean(js > 0):.1%} of pairs share a gene; "
          f"gene-disjoint subfamily has {len(chosen)} sets")

    st = json.load(open(f"{ROOT}/results/SECRETOME_TRANSFER.json"))
    dis = {}
    for name in datasets:
        d = [st["per_dataset_per_set"][name][s]["mean_z_secreted"]
             - st["per_dataset_per_set"][name][s]["mean_z_other"]
             for s in chosen
             if s in st["per_dataset_per_set"][name]
             and st["per_dataset_per_set"][name][s]["n_secreted"] >= 3
             and st["per_dataset_per_set"][name][s]["n_other"] >= 3]
        dis[name] = {"level": datasets[name][0], "n_sets": len(d),
                     "median_secreted_minus_other": float(np.median(d)) if d else None,
                     "all_values": [round(v, 4) for v in d]}
    out["D_disjoint_subfamily_result"] = dis

    # ---- E: per gene, transcript-level against protein-level shift
    TX = [n for n in datasets if datasets[n][0] == "transcriptome"]
    PR = [n for n in datasets if datasets[n][0] == "proteome"]
    common = set(shifts[TX[0]])
    for n in TX[1:] + PR:
        common &= set(shifts[n])
    common = sorted(common)
    tz = {g: float(np.mean([shifts[n][g]["z"] for n in TX])) for g in common}
    pz = {g: float(np.mean([shifts[n][g]["z"] for n in PR])) for g in common}
    for lab, genes in (("secreted", [g for g in common if g in secreted]),
                       ("other", [g for g in common if g not in secreted])):
        a = [tz[g] for g in genes]; b = [pz[g] for g in genes]
        rho, p = stats.spearmanr(a, b)
        out["E_cross_level_agreement"][lab] = {
            "n_genes": len(genes), "spearman_rho": float(rho), "p": float(p),
            "mean_transcript_shift": float(np.mean(a)),
            "mean_protein_shift": float(np.mean(b)),
            "share_same_sign": float(np.mean(
                [np.sign(x) == np.sign(y) for x, y in zip(a, b)]))}
        e = out["E_cross_level_agreement"][lab]
        print(f"E  {lab:9s} n={len(genes):5d}  tx {e['mean_transcript_shift']:+.3f}  "
              f"pr {e['mean_protein_shift']:+.3f}  rho={rho:+.3f} (P={p:.2g})  "
              f"same sign {e['share_same_sign']:.3f}")

    with open(f"{ROOT}/results/SECRETOME_VERIFY.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"\nwrote {ROOT}/results/SECRETOME_VERIFY.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
