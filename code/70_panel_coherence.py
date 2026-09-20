#!/usr/bin/env python3
"""70_panel_coherence.py -- is a tumor-adjacent covariate panel one axis at both
measurement levels?

Post hoc. Nothing here was pre-registered. It was written after GigaScience
declined the parent manuscript, to re-derive from the raw deposits the one claim
the parent manuscript rested on: that one published liver proteome violates the
hepatocyte-identity premise while another satisfies it. The re-derivation does
not support that claim, and what it finds instead is reported here.

Every number is written to PANEL_COHERENCE.json. Nothing is retyped anywhere.
"""
from __future__ import annotations
import csv, json, os, sys
import numpy as np, pandas as pd
from scipy import stats

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_cands = [_here, _os.path.join(_here, "rsi", "code"), _os.path.join(_here, "code"),
          _os.path.join(_os.path.dirname(_here), "code")]
if _os.environ.get("IR_ROOT"):
    _cands.insert(0, _os.path.join(_os.environ["IR_ROOT"], "code"))
for _c in _cands:
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        if _c not in _sys.path:
            _sys.path.insert(0, _c)
        break
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, GENESETS as IR_GENESETS, \
    DATA as IR_DATA, CODE as IR_CODE, FIGURES as IR_FIGURES, DOCS as IR_DOCS

csv.field_size_limit(10_000_000)
OUT = f"{IR_RESULTS}/PANEL_COHERENCE.json"

CLASS = {
    "secreted plasma protein": ["ALB","TTR","TF","SERPINA1","AHSG","APOH",
                                "FGA","FGB","FGG","F2"],
    "metabolic enzyme":        ["CPS1","OTC","ARG1","TAT","G6PC1","PCK1"],
    "transcription factor":    ["HNF4A","HNF1A","FOXA1","FOXA2","NR1H4"],
    "surface receptor":        ["ASGR1"],
}
D1 = [g for v in CLASS.values() for g in v]
OF_CLASS = {g: c for c, gs in CLASS.items() for g in gs}
PLASMA_OFF_PANEL = ["HP","HPX","APOA1","APOA2","C3","ORM1","A2M","SERPINC1",
                    "AMBP","ITIH4","APCS","CFB"]
ERYTHROCYTE = ["HBB","HBA1","HBD","CA1","SLC4A1","PRDX2"]
MIN_PAIRS = 10          # a member is used only if quantified in both arms of 10+ pairs
MIN_MEMBERS = 10        # a gene set is scored only if 10+ members clear that bar


# ---------------------------------------------------------------- loaders ---
def load_jiang(mode="raw"):
    """PRIDE PXD006512, MaxQuant proteinGroups.txt, iBAQ columns.

    mode: raw (log2 of iBAQ, zeros treated as not quantified), median (log2 then
    per-sample median centering over quantified proteins), quantile (the parent
    study's pipeline: quantile normalization with tied zeros, then log2)."""
    src = (f"{IR_DOCS}/jiang/MaxQuant results and its relative supplementary "
           "materials/proteinGroups.txt")
    with open(src, newline="", encoding="latin-1") as fh:
        rdr = csv.reader(fh, delimiter="\t"); hdr = next(rdr)
        idx = {c: i for i, c in enumerate(hdr)}
        ib = [c for c in hdr if c.startswith("iBAQ ") and c != "iBAQ peptides"]
        fl = [k for k in ("Only identified by site", "Reverse",
                          "Potential contaminant") if k in idx]
        rows, genes = [], []
        for r in rdr:
            if any(r[idx[k]].strip() == "+" for k in fl):
                continue
            g = r[idx["Gene names"]].split(";")[0].strip()
            if g:
                genes.append(g); rows.append([r[idx[c]] for c in ib])
    m = pd.DataFrame(rows, index=genes, columns=[c[5:] for c in ib])
    m = m.apply(pd.to_numeric, errors="coerce").fillna(0.0).groupby(level=0).max()
    m = m.loc[:, m.sum(axis=0) > 0]
    m = m[(m > 0).sum(axis=1) > 0]
    if mode == "quantile":
        ranks = m.rank(method="average", axis=0)
        ms = np.sort(m.values, axis=0).mean(axis=1)
        grid = np.arange(1, m.shape[0] + 1)
        q = pd.DataFrame(np.column_stack([np.interp(ranks[c].values, grid, ms)
                                          for c in m.columns]),
                         index=m.index, columns=m.columns)
        lg = np.log2(q + 1.0).replace(0, np.nan)
    else:
        lg = np.log2(m.replace(0, np.nan))
        if mode == "median":
            lg = lg.sub(lg.median(axis=0), axis=1)
    ids = sorted(set(c[:-1] for c in lg.columns if c.endswith("T")) &
                 set(c[:-1] for c in lg.columns if c.endswith("P")))
    return lg, [i + "T" for i in ids], [i + "P" for i in ids]


def load_gao():
    """Gao et al. 2019 Table S1 as the authors published it, 6,478 proteins."""
    m = pd.read_csv(f"{IR_DOCS}/gao_proteins.tsv", sep="\t", index_col=0,
                    low_memory=False).groupby(level=0).max()
    ids = sorted(set(c[1:] for c in m.columns if c.startswith("T")) &
                 set(c[1:] for c in m.columns if c.startswith("N")))
    return m, ["T" + i for i in ids], ["N" + i for i in ids]


def load_tx(matrix, pheno, tumor_label=None):
    m = pd.read_csv(matrix, sep="\t", index_col=0).groupby(level=0).max()
    ph = pd.read_csv(pheno, sep="\t", index_col=0)
    ph = ph[ph.index.isin(m.columns)]
    labs = list(ph["tissue"].unique())
    tl = tumor_label or next(l for l in labs if "tumor" in l.lower()
                             and "non" not in l.lower() and "adjacent" not in l.lower())
    al = next(l for l in labs if l != tl)
    t, a = ph[ph.tissue == tl], ph[ph.tissue == al]
    ids = sorted(set(t.patient_id) & set(a.patient_id))
    return (m,
            [t[t.patient_id == p].index[0] for p in ids],
            [a[a.patient_id == p].index[0] for p in ids])


# ----------------------------------------------------------------- shifts ---
def member_shifts(mat, tc, pc, genes):
    """Paired tumor minus adjacent, as a z on the pooled spread of both arms."""
    out = {}
    for g in genes:
        if g not in mat.index:
            continue
        v = mat.loc[g]
        if isinstance(v, pd.DataFrame):
            v = v.iloc[0]
        x = v[tc].values.astype(float); y = v[pc].values.astype(float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < MIN_PAIRS:
            continue
        sd = np.concatenate([x[ok], y[ok]]).std(ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        out[g] = {"z": float((x[ok] - y[ok]).mean() / sd), "n_pairs": int(ok.sum())}
    return out


def composite_sample_test(mat, tc, pc, genes):
    """The test the parent manuscript used: score each sample on the panel, then
    a paired test over samples. This is what turns a near-zero panel into a
    significant verdict whose sign is arbitrary."""
    present = [g for g in genes if g in mat.index]
    sub = mat.loc[present]
    sub = sub[~sub.index.duplicated()]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
    score = z.mean(axis=0, skipna=True)
    x = score[tc].values.astype(float); y = score[pc].values.astype(float)
    ok = ~(np.isnan(x) | np.isnan(y))
    d = x[ok] - y[ok]
    w = stats.wilcoxon(d) if ok.sum() > 10 else None
    return {"n_members_scored": len(present), "n_pairs": int(ok.sum()),
            "mean_paired_delta": float(d.mean()),
            "direction": "HIGHER in tumor" if d.mean() > 0 else "LOWER in tumor",
            "wilcoxon_p": float(w.pvalue) if w else None}


def coherence(zs):
    """Two ways of saying whether the members of a set point the same way.

    sign_agreement is the share of members on the majority side.
    cancellation is |mean| over mean|.|: 1 when nothing cancels, 0 when the set
    averages to nothing because its members oppose each other."""
    a = np.array(zs, float)
    if a.size == 0:
        return None
    pos = int((a > 0).sum()); neg = int((a < 0).sum())
    return {"n": int(a.size),
            "sign_agreement": float(max(pos, neg) / a.size),
            "cancellation": float(abs(a.mean()) / np.abs(a).mean())
            if np.abs(a).mean() > 0 else None,
            "mean_z": float(a.mean())}


def main() -> int:
    datasets = {
        "GSE14520 (array transcriptome)":
            ("transcriptome",) + load_tx(f"{IR_DATA}/GSE14520_symbols.tsv.gz",
                                         f"{IR_RESULTS}/GSE14520/phenotype.tsv",
                                         "HCC tumor"),
        "GSE76427 (array transcriptome)":
            ("transcriptome",) + load_tx(f"{IR_DATA}/GSE76427_symbols.tsv.gz",
                                         f"{IR_RESULTS}/GSE76427/phenotype.tsv"),
        "TCGA-LIHC (RNA sequencing)":
            ("transcriptome",) + load_tx(f"{IR_DATA}/TCGA_LIHC_symbols.tsv.gz",
                                         f"{IR_RESULTS}/TCGA_LIHC/phenotype.tsv",
                                         "tumor"),
        "Gao 2019 proteome":   ("proteome",) + load_gao(),
        "Jiang 2019 proteome": ("proteome",) + load_jiang("raw"),
    }

    res = {"what": "post hoc; not pre-registered",
           "min_pairs_per_member": MIN_PAIRS,
           "min_members_per_set": MIN_MEMBERS,
           "panel_classes": CLASS, "datasets": {}}

    sets = json.load(open(f"{IR_GENESETS}/eligible.json"))

    for name, (level, mat, tc, pc) in datasets.items():
        d = member_shifts(mat, tc, pc, D1)
        by_class = {}
        for cls, gs in CLASS.items():
            v = [d[g]["z"] for g in gs if g in d]
            if v:
                by_class[cls] = {"n": len(v), "mean_z": float(np.mean(v)),
                                 "members": {g: round(d[g]["z"], 4)
                                             for g in gs if g in d}}
        sec = [d[g]["z"] for g in CLASS["secreted plasma protein"] if g in d]
        met = [d[g]["z"] for g in CLASS["metabolic enzyme"] if g in d]
        split = None
        if len(sec) > 2 and len(met) > 2:
            u = stats.mannwhitneyu(sec, met, alternative="two-sided")
            split = {"U": float(u.statistic), "p": float(u.pvalue),
                     "separation_complete": bool(min(sec) > max(met)
                                                 or max(sec) < min(met))}
        allz = [v["z"] for v in d.values()]
        t = stats.ttest_1samp(allz, 0.0)

        # the two controls that would explain the split if they were the cause
        ery = member_shifts(mat, tc, pc, ERYTHROCYTE)
        off = member_shifts(mat, tc, pc, PLASMA_OFF_PANEL)

        # abundance, the obvious alternative
        X = mat[tc].astype(float).values; Y = mat[pc].astype(float).values
        ok = ~(np.isnan(X) | np.isnan(Y)); keep = ok.sum(axis=1) >= MIN_PAIRS
        dd = np.where(ok, X - Y, np.nan)[keep]
        sd = np.nanstd(np.concatenate([np.where(ok, X, np.nan)[keep],
                                       np.where(ok, Y, np.nan)[keep]], axis=1),
                       axis=1, ddof=1)
        zz = np.nanmean(dd, axis=1) / np.where(sd == 0, np.nan, sd)
        ab = np.nanmean(np.where(ok, (X + Y) / 2, np.nan)[keep], axis=1)
        gi = mat.index[keep]
        good = ~(np.isnan(zz) | np.isnan(ab))
        rho = stats.spearmanr(ab[good], zz[good])
        pct = (pd.Series(ab[good], index=gi[good]).rank(pct=True) * 100)

        # every eligible liver set, at this level
        coh = {}
        for sname, members in sets.items():
            ds = member_shifts(mat, tc, pc, members)
            if len(ds) >= MIN_MEMBERS:
                coh[sname] = coherence([v["z"] for v in ds.values()])

        res["datasets"][name] = {
            "level": level, "n_pairs": len(tc),
            "panel_members": {g: {**d[g], "class": OF_CLASS[g]} for g in d},
            "panel_by_class": by_class,
            "secreted_vs_metabolic": split,
            "panel_mean_z": float(np.mean(allz)),
            "panel_t_over_members": {"t": float(t.statistic), "p": float(t.pvalue)},
            "panel_composite_sample_test": composite_sample_test(mat, tc, pc, D1),
            "erythrocyte_mean_z": float(np.mean([v["z"] for v in ery.values()]))
            if ery else None,
            "off_panel_plasma_mean_z": float(np.mean([v["z"] for v in off.values()]))
            if off else None,
            "abundance": {
                "spearman_shift_vs_abundance": float(rho.statistic),
                "p": float(rho.pvalue),
                "secreted_percentile_median":
                    float(np.median([pct[g] for g in CLASS["secreted plasma protein"]
                                     if g in pct.index])),
                "metabolic_percentile_median":
                    float(np.median([pct[g] for g in CLASS["metabolic enzyme"]
                                     if g in pct.index]))},
            "gene_set_coherence": coh,
        }
        print(f"  {name}: {len(tc)} pairs, {len(d)}/{len(D1)} panel members, "
              f"{len(coh)} gene sets scored")

    # Jiang under three preprocessings, to rule the split out as an artifact
    res["jiang_normalization_sensitivity"] = {}
    for mode in ("raw", "median", "quantile"):
        m, tc, pc = load_jiang(mode)
        d = member_shifts(m, tc, pc, D1)
        sec = [d[g]["z"] for g in CLASS["secreted plasma protein"] if g in d]
        met = [d[g]["z"] for g in CLASS["metabolic enzyme"] if g in d]
        res["jiang_normalization_sensitivity"][mode] = {
            "secreted_mean_z": float(np.mean(sec)), "metabolic_mean_z": float(np.mean(met)),
            "panel_mean_z": float(np.mean([v["z"] for v in d.values()])),
            "separation_complete": bool(min(sec) > max(met))}

    # the paired comparison of coherence between levels, over the sets both carry
    tx = [n for n, v in res["datasets"].items() if v["level"] == "transcriptome"]
    pr = [n for n, v in res["datasets"].items() if v["level"] == "proteome"]
    shared = set.intersection(*[set(res["datasets"][n]["gene_set_coherence"])
                                for n in tx + pr])
    def mean_over(ns, sname, key):
        return float(np.mean([res["datasets"][n]["gene_set_coherence"][sname][key]
                              for n in ns]))
    rows = [{"set": s,
             "transcriptome_cancellation": mean_over(tx, s, "cancellation"),
             "proteome_cancellation": mean_over(pr, s, "cancellation"),
             "transcriptome_sign_agreement": mean_over(tx, s, "sign_agreement"),
             "proteome_sign_agreement": mean_over(pr, s, "sign_agreement")}
            for s in sorted(shared)]
    a = np.array([r["transcriptome_cancellation"] for r in rows])
    b = np.array([r["proteome_cancellation"] for r in rows])
    sa = np.array([r["transcriptome_sign_agreement"] for r in rows])
    sb = np.array([r["proteome_sign_agreement"] for r in rows])
    res["coherence_by_level"] = {
        "n_sets_scored_at_both_levels": len(rows),
        "cancellation_median_transcriptome": float(np.median(a)),
        "cancellation_median_proteome": float(np.median(b)),
        "cancellation_wilcoxon_p": float(stats.wilcoxon(a, b).pvalue),
        "n_sets_cancellation_lower_in_proteome": int((b < a).sum()),
        "sign_agreement_median_transcriptome": float(np.median(sa)),
        "sign_agreement_median_proteome": float(np.median(sb)),
        "sign_agreement_wilcoxon_p": float(stats.wilcoxon(sa, sb).pvalue),
        "per_set": rows,
    }
    json.dump(res, open(OUT, "w"), indent=1)
    c = res["coherence_by_level"]
    print(f"\n  {c['n_sets_scored_at_both_levels']} gene sets scored at both levels")
    print(f"  cancellation ratio: transcriptome {c['cancellation_median_transcriptome']:.3f}"
          f"  proteome {c['cancellation_median_proteome']:.3f}"
          f"  (P = {c['cancellation_wilcoxon_p']:.3g},"
          f" lower in proteome for {c['n_sets_cancellation_lower_in_proteome']})")
    print(f"  sign agreement:     transcriptome {c['sign_agreement_median_transcriptome']:.3f}"
          f"  proteome {c['sign_agreement_median_proteome']:.3f}"
          f"  (P = {c['sign_agreement_wilcoxon_p']:.3g})")
    print(f"\n  wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
