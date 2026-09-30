#!/usr/bin/env python3
"""75_secretome_transfer.py

Does the secretome content of a liver gene set predict whether it transfers
from transcript level to protein level?

The 22-gene hepatocyte identity panel splits in both deposited liver proteomes:
its secreted plasma proteins rise in tumor while its metabolic enzymes fall,
and the two halves cancel. That observation was made on one hand-built panel.
This script asks whether the same thing is a general property of liver gene
sets, using the 114 published liver signatures that are already scored in all
five datasets.

Annotation source: the NCBI GRCh38.p13 gene annotation table already deposited
in this archive (data/Human.GRCh38.p13.annot.tsv.gz), GO cellular component
column. A gene is called secreted when it carries GO:0005615, extracellular
space. GO:0070062, extracellular exosome, is deliberately NOT used: it is
carried by cytosolic enzymes such as PCK1 and ALDOB and does not separate
secreted from intracellular proteins.

Writes results/SECRETOME_TRANSFER.json.
"""
import csv, gzip, importlib.util, json, math, os, sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

spec = importlib.util.spec_from_file_location(
    "panel_coherence", os.path.join(HERE, "70_panel_coherence.py"))
pc70 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc70)

ES = "GO:0005615"   # extracellular space
ER = "GO:0005576"   # extracellular region (broader; sensitivity analysis only)
MIN_SIDE = 3        # a set contributes to the within-set test only with 3+ per side


# ---------------------------------------------------------------- annotation

def load_components():
    path = os.path.join(ROOT, "data", "Human.GRCh38.p13.annot.tsv.gz")
    comp = {}
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ids = (row.get("GOComponentID") or "").split("///")
            comp[row["Symbol"]] = {i for i in ids if i.startswith("GO:")}
    return comp


def annotation_check(comp):
    """Report what the annotation contains before it is used for anything."""
    positives = ["ALB", "TTR", "APOH", "F2", "SERPINA1", "AHSG", "FGA", "FGB",
                 "FGG", "TF", "APOA1", "HP", "SERPINC1", "ITIH4", "A2M", "AMBP",
                 "APCS", "ORM1", "C3", "HPX", "APOA2", "CFB"]
    negatives = ["HNF4A", "HNF1A", "FOXA1", "FOXA2", "NR1H4", "CPS1", "OTC",
                 "TAT", "PCK1", "ALDOB", "G6PC1", "CYP2E1", "SLC2A2", "ARG1",
                 "ASGR1"]
    def call(g, broad=False):
        s = comp.get(g, set())
        return (ES in s) or (broad and ER in s)
    out = {"source": "data/Human.GRCh38.p13.annot.tsv.gz, GOComponentID",
           "secreted_term": ES, "genes_in_annotation": len(comp)}
    for lab, broad in (("strict_GO_0005615", False), ("broad_incl_GO_0005576", True)):
        tp = [g for g in positives if call(g, broad)]
        fp = [g for g in negatives if call(g, broad)]
        out[lab] = {
            "known_plasma_proteins_called_secreted": f"{len(tp)}/{len(positives)}",
            "missed": [g for g in positives if not call(g, broad)],
            "known_intracellular_called_secreted": f"{len(fp)}/{len(negatives)}",
            "false_positives": fp}
    return out


# ------------------------------------------------------------------ analysis

def side_means(shifts, secreted):
    """mean z of the secreted and the non-secreted members of one scored set."""
    s = [v["z"] for g, v in shifts.items() if g in secreted]
    n = [v["z"] for g, v in shifts.items() if g not in secreted]
    return s, n


def main() -> int:
    comp = load_components()
    check = annotation_check(comp)
    print(json.dumps(check, indent=2))

    secreted = {g for g, s in comp.items() if ES in s}
    secreted_broad = {g for g, s in comp.items() if ES in s or ER in s}

    datasets = {
        "GSE14520 (microarray)":
            ("transcriptome",) + pc70.load_tx(
                f"{pc70.IR_DATA}/GSE14520_symbols.tsv.gz",
                f"{pc70.IR_RESULTS}/GSE14520/phenotype.tsv"),
        "GSE76427 (microarray)":
            ("transcriptome",) + pc70.load_tx(
                f"{pc70.IR_DATA}/GSE76427_symbols.tsv.gz",
                f"{pc70.IR_RESULTS}/GSE76427/phenotype.tsv"),
        "TCGA-LIHC (RNA sequencing)":
            ("transcriptome",) + pc70.load_tx(
                f"{pc70.IR_DATA}/TCGA_LIHC_symbols.tsv.gz",
                f"{pc70.IR_RESULTS}/TCGA_LIHC/phenotype.tsv", "tumor"),
        "Gao 2019 proteome":   ("proteome",) + pc70.load_gao(),
        "Jiang 2019 proteome": ("proteome",) + pc70.load_jiang("raw"),
    }

    sets = json.load(open(f"{pc70.IR_GENESETS}/eligible.json"))
    sets = dict(sets)
    sets["__PANEL__"] = list(pc70.D1)

    per = {name: {} for name in datasets}
    for name, (level, mat, tc, pc) in datasets.items():
        for sname, members in sets.items():
            sh = pc70.member_shifts(mat, tc, pc, members)
            if len(sh) < pc70.MIN_MEMBERS:
                continue
            zs = [v["z"] for v in sh.values()]
            s, n = side_means(sh, secreted)
            sb, nb = side_means(sh, secreted_broad)
            per[name][sname] = {
                "level": level, "n": len(sh),
                "coherence": pc70.coherence(zs),
                "mean_z": float(np.mean(zs)),
                "n_secreted": len(s), "n_other": len(n),
                "mean_z_secreted": float(np.mean(s)) if s else None,
                "mean_z_other": float(np.mean(n)) if n else None,
                "n_secreted_broad": len(sb),
                "mean_z_secreted_broad": float(np.mean(sb)) if sb else None,
                "mean_z_other_broad": float(np.mean(nb)) if nb else None,
            }

    TX = [k for k, v in datasets.items() if v[0] == "transcriptome"]
    PR = [k for k, v in datasets.items() if v[0] == "proteome"]
    scored_everywhere = [s for s in sets
                         if all(s in per[n] for n in datasets)]
    print(f"\nsets scored in all five datasets: "
          f"{len([s for s in scored_everywhere if s != '__PANEL__'])} "
          f"(+ the panel)")

    # ---- H1: does secreted content predict the coherence drop?
    rows = []
    for s in scored_everywhere:
        ctx = float(np.mean([per[n][s]["coherence"]["cancellation"] for n in TX]))
        cpr = float(np.mean([per[n][s]["coherence"]["cancellation"] for n in PR]))
        frac = float(np.mean([per[n][s]["n_secreted"] / per[n][s]["n"] for n in PR]))
        rows.append({"set": s, "drop": ctx - cpr,
                     "coherence_transcriptome": ctx, "coherence_proteome": cpr,
                     "secreted_fraction_in_proteomes": frac,
                     "n_scored_proteome": float(np.mean([per[n][s]["n"] for n in PR]))})
    body = [r for r in rows if r["set"] != "__PANEL__"]
    panel = [r for r in rows if r["set"] == "__PANEL__"][0]

    x = np.array([r["secreted_fraction_in_proteomes"] for r in body])
    y = np.array([r["drop"] for r in body])
    rho, prho = stats.spearmanr(x, y)
    pear, ppear = stats.pearsonr(x, y)
    size = np.log10([r["n_scored_proteome"] for r in body])
    X = np.column_stack([np.ones_like(x), x, size])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (resid @ resid) / dof))
    tstat = beta / se
    pvals = [float(2 * stats.t.sf(abs(t), dof)) for t in tstat]

    h1 = {"n_sets": len(body),
          "spearman_rho": float(rho), "spearman_p": float(prho),
          "pearson_r": float(pear), "pearson_p": float(ppear),
          "ols_with_log10_set_size": {
              "beta_secreted_fraction": float(beta[1]),
              "p_secreted_fraction": pvals[1],
              "beta_log10_size": float(beta[2]), "p_log10_size": pvals[2]},
          "panel": {"secreted_fraction": panel["secreted_fraction_in_proteomes"],
                    "drop": panel["drop"]},
          "secreted_fraction_percentile_of_panel": float(
              100.0 * np.mean(x <= panel["secreted_fraction_in_proteomes"]))}
    print("\nH1  drop ~ secreted fraction: "
          f"rho={rho:+.3f} (P={prho:.2g}), r={pear:+.3f} (P={ppear:.2g}); "
          f"OLS beta={beta[1]:+.3f} (P={pvals[1]:.2g})")

    # ---- H2: within a set, do secreted members sit above the others?
    h2 = {}
    for name in datasets:
        diffs, names = [], []
        for s, v in per[name].items():
            if s == "__PANEL__":
                continue
            if v["n_secreted"] >= MIN_SIDE and v["n_other"] >= MIN_SIDE:
                diffs.append(v["mean_z_secreted"] - v["mean_z_other"])
                names.append(s)
        d = np.array(diffs)
        w = stats.wilcoxon(d) if len(d) > 5 else None
        pv = per[name].get("__PANEL__")
        h2[name] = {
            "level": datasets[name][0], "n_sets_with_both_sides": len(d),
            "median_secreted_minus_other": float(np.median(d)) if len(d) else None,
            "mean_secreted_minus_other": float(np.mean(d)) if len(d) else None,
            "share_positive": float(np.mean(d > 0)) if len(d) else None,
            "wilcoxon_p": float(w.pvalue) if w is not None else None,
            "panel_secreted_minus_other":
                (pv["mean_z_secreted"] - pv["mean_z_other"])
                if pv and pv["n_secreted"] and pv["n_other"] else None,
            "panel_percentile_among_sets": float(100.0 * np.mean(
                d <= (pv["mean_z_secreted"] - pv["mean_z_other"])))
                if pv and len(d) else None,
        }
        print(f"H2  {name:32s} {h2[name]['level'][:5]}  "
              f"n={h2[name]['n_sets_with_both_sides']:3d}  "
              f"median={h2[name]['median_secreted_minus_other']:+.3f}  "
              f"share>0={h2[name]['share_positive']:.3f}  "
              f"P={h2[name]['wilcoxon_p']:.2g}")

    # ---- H2 sensitivity under the broad annotation
    h2b = {}
    for name in datasets:
        d = np.array([v["mean_z_secreted_broad"] - v["mean_z_other_broad"]
                      for s, v in per[name].items()
                      if s != "__PANEL__" and v["n_secreted_broad"] >= MIN_SIDE
                      and v["n"] - v["n_secreted_broad"] >= MIN_SIDE])
        w = stats.wilcoxon(d) if len(d) > 5 else None
        h2b[name] = {"n": len(d), "median": float(np.median(d)) if len(d) else None,
                     "share_positive": float(np.mean(d > 0)) if len(d) else None,
                     "wilcoxon_p": float(w.pvalue) if w is not None else None}

    out = {"what": "post hoc; not pre-registered",
           "annotation_check": check,
           "n_secreted_genes_in_annotation": len(secreted),
           "min_members_per_side": MIN_SIDE,
           "h1_drop_vs_secreted_fraction": h1,
           "h1_per_set": sorted(body, key=lambda r: -r["drop"]),
           "h2_within_set_secreted_minus_other": h2,
           "h2_broad_annotation_sensitivity": h2b,
           "per_dataset_per_set": per}
    os.makedirs(f"{ROOT}/results", exist_ok=True)
    with open(f"{ROOT}/results/SECRETOME_TRANSFER.json", "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"\nwrote {ROOT}/results/SECRETOME_TRANSFER.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
