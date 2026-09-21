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
# Erythrocyte proteins report the CELLULAR blood compartment: how much red cell
# content a piece of tissue carries. The alternative this paper has to exclude is
# plasma EXUDATION, a different compartment, so it is tested with plasma proteins
# the liver does not make. Immunoglobulins and JCHAIN are made by plasma cells.
IMMUNOGLOBULIN = ["IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHM",
                  "IGKC","IGLC1","IGLC2","JCHAIN"]
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


def composite_sample_test(mat, tc, pc, genes, restrict=None, scale="paired"):
    """The test the parent manuscript used: score each sample on the panel, then
    a paired test over samples. This is what turns a near-zero panel into a
    significant verdict whose sign is arbitrary.

    restrict, when given, is the member list that cleared MIN_PAIRS, so the
    composite is built from exactly the members every other number here uses.
    Without it the composite silently admits members quantified in a handful of
    pairs, which is how the parent manuscript scored them."""
    present = [g for g in genes if g in mat.index]
    if restrict is not None:
        present = [g for g in present if g in set(restrict)]
    sub = mat.loc[present]
    sub = sub[~sub.index.duplicated()]
    # Which samples define the scale each member is standardized on. "paired"
    # uses the two arms of the pairs and nothing else, which is the only choice
    # that means the same thing in every matrix here: GSE14520 carries 445
    # columns for 213 pairs and GSE76427 carries 167 for 52, so standardizing
    # over every column puts unpaired samples into the scale of three datasets
    # and not the other two, in a paper that compares that quantity across them.
    ref = sub[list(tc) + list(pc)] if scale == "paired" else sub
    z = sub.sub(ref.mean(axis=1), axis=0).div(ref.std(axis=1, ddof=1), axis=0)
    score = z.mean(axis=0, skipna=True)
    x = score[tc].values.astype(float); y = score[pc].values.astype(float)
    ok = ~(np.isnan(x) | np.isnan(y))
    d = x[ok] - y[ok]
    w = stats.wilcoxon(d) if ok.sum() > 10 else None
    return {"n_members_scored": len(present), "n_pairs": int(ok.sum()),
            "standardized_on": scale,
            "mean_paired_delta": float(d.mean()),
            "direction": "HIGHER in tumor" if d.mean() > 0 else "LOWER in tumor",
            "wilcoxon_p": float(w.pvalue) if w else None}


def missingness(mat, tc, pc):
    """How much is unquantified, and is it the same in the two arms?

    A complete-case paired difference conditions on the protein being seen in
    both arms. When one arm detects less than the other, that conditioning is
    not symmetric, and a reviewer of a mass spectrometry reuse paper will ask.
    """
    t = mat[list(tc)].astype(float).values
    p = mat[list(pc)].astype(float).values
    return {"fraction_unquantified_tumor": float(np.isnan(t).mean()),
            "fraction_unquantified_adjacent": float(np.isnan(p).mean()),
            "mean_features_per_tumor_sample": float((~np.isnan(t)).sum(axis=0).mean()),
            "mean_features_per_adjacent_sample": float((~np.isnan(p)).sum(axis=0).mean())}


def matrix_diagnostic(mat, cols):
    """Is this matrix on an abundance scale at all?

    A per-protein centered log ratio has about half its values below zero and a
    row median near zero for most rows; an abundance matrix does not. The mean
    of a row is then a ratio, not an amount, and no statement about which
    proteins are abundant can be read off it."""
    v = mat[cols].astype(float)
    arr = v.values
    fin = np.isfinite(arr)
    rowmed = np.nanmedian(np.where(fin, arr, np.nan), axis=1)
    rowmed = rowmed[np.isfinite(rowmed)]
    frac_neg = float((arr[fin] < 0).mean()) if fin.any() else float("nan")
    return {"fraction_of_values_below_zero": frac_neg,
            "median_of_per_protein_medians": float(np.median(rowmed)),
            "share_of_proteins_with_median_within_0.5_of_zero":
                float((np.abs(rowmed) < 0.5).mean()),
            "on_an_abundance_scale": bool(frac_neg < 0.2)}


def unfiltered_classes(mat, tc, pc, min_pairs=1):
    """Class means over every panel member present in the matrix."""
    out = {}
    for cls, gs in CLASS.items():
        z = []
        for g in gs:
            if g not in mat.index:
                continue
            v = mat.loc[g]
            if isinstance(v, pd.DataFrame):
                v = v.iloc[0]
            x = v[tc].values.astype(float); y = v[pc].values.astype(float)
            ok = ~(np.isnan(x) | np.isnan(y))
            if ok.sum() < min_pairs:
                continue
            sd = np.concatenate([x[ok], y[ok]]).std(ddof=1)
            if np.isfinite(sd) and sd > 0:
                z.append(float((x[ok] - y[ok]).mean() / sd))
        if z:
            out[cls] = {"n": len(z), "mean_z": float(np.mean(z)),
                        "min": float(min(z)), "max": float(max(z))}
    sec = out.get("secreted plasma protein"); met = out.get("metabolic enzyme")
    return {"by_class": out,
            "separation_complete": bool(sec and met and sec["min"] > met["max"])}


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
            import math as _m
            floor = 2.0 / _m.comb(len(sec) + len(met), min(len(sec), len(met)))
            split = {"U": float(u.statistic), "p": float(u.pvalue),
                     "n_secreted": len(sec), "n_metabolic": len(met),
                     "min_attainable_p": float(floor),
                     "p_is_at_floor": bool(abs(u.pvalue - floor) < 1e-12),
                     "separation_complete": bool(min(sec) > max(met)
                                                 or max(sec) < min(met))}
        # Post hoc, and named as such: the sign partition the members actually
        # fall into, which is not the contrast the class labels anticipated.
        UP_CLASSES = {"secreted plasma protein", "transcription factor"}
        up = [g for g in d if OF_CLASS[g] in UP_CLASSES]
        dn = [g for g in d if OF_CLASS[g] not in UP_CLASSES]
        tab = [[sum(1 for g in up if d[g]["z"] > 0), sum(1 for g in up if d[g]["z"] <= 0)],
               [sum(1 for g in dn if d[g]["z"] > 0), sum(1 for g in dn if d[g]["z"] <= 0)]]
        fe = stats.fisher_exact(tab)
        import math as _m2
        _n_up, _n_dn = len(up), len(dn)
        _floor = 1.0 / _m2.comb(_n_up + _n_dn, min(_n_up, _n_dn))
        partition = {
            "note": "post hoc grouping, defined after inspecting the member shifts; "
                    "it merges two of the four pre-specified classes against the "
                    "other two, so the choice was made among seven possible merges",
            "group_up": sorted(up), "group_down": sorted(dn),
            "table_up_pos_neg_down_pos_neg": tab,
            "fisher_p": float(fe[1]),
            "min_attainable_p": float(_floor),
            "p_is_at_floor": bool(abs(fe[1] - _floor) < 1e-12),
            "partition_is_perfect": bool(tab[0][1] == 0 and tab[1][0] == 0)}

        allz = [v["z"] for v in d.values()]
        t = stats.ttest_1samp(allz, 0.0)
        _a = np.array(allz, float)
        _se = _a.std(ddof=1) / np.sqrt(_a.size)
        _t = stats.t.ppf(0.975, _a.size - 1)
        member_ci = [float(_a.mean() - _t * _se), float(_a.mean() + _t * _se)]

        # the two controls that would explain the split if they were the cause
        ery = member_shifts(mat, tc, pc, ERYTHROCYTE)
        off = member_shifts(mat, tc, pc, PLASMA_OFF_PANEL)
        igg = member_shifts(mat, tc, pc, IMMUNOGLOBULIN)

        # abundance, the obvious alternative. Read the scale first: a matrix of
        # per-protein centered ratios carries no abundance information at all,
        # and the percentiles below would then rank ratios, not amounts.
        diag = matrix_diagnostic(mat, list(tc) + list(pc))
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
        # A correlation alone cannot say whether the gap in abundance BETWEEN the
        # two classes is large enough to produce the gap in shift between them.
        # This asks it directly: take every protein sitting where each class
        # sits, and see how far apart those two neighborhoods move.
        zser = pd.Series(zz[good], index=gi[good])
        band = None
        if len(pct) > 50:
            sp = np.median([pct[g] for g in CLASS["secreted plasma protein"]
                            if g in pct.index])
            mp = np.median([pct[g] for g in CLASS["metabolic enzyme"]
                            if g in pct.index])
            def nbhd(centre, half=1.5):
                sel = pct[(pct >= centre - half) & (pct <= centre + half)].index
                return zser.reindex(sel).dropna()
            ns, nm = nbhd(sp), nbhd(mp)
            if len(ns) > 5 and len(nm) > 5:
                band = {"secreted_neighborhood_percentile": float(sp),
                        "metabolic_neighborhood_percentile": float(mp),
                        "secreted_neighborhood_mean_z": float(ns.mean()),
                        "metabolic_neighborhood_mean_z": float(nm.mean()),
                        "n_secreted_neighborhood": int(len(ns)),
                        "n_metabolic_neighborhood": int(len(nm)),
                        "gap_expected_from_abundance":
                            float(ns.mean() - nm.mean())}
                obs_gap = (by_class.get("secreted plasma protein", {}).get("mean_z")
                           or 0.0) - (by_class.get("metabolic enzyme", {}).get("mean_z")
                                      or 0.0)
                if obs_gap:
                    band["observed_class_gap"] = float(obs_gap)
                    band["share_of_class_gap_explained"] = float(
                        band["gap_expected_from_abundance"] / obs_gap)

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
            "class_difference_secreted_minus_metabolic":
                float(by_class["secreted plasma protein"]["mean_z"]
                      - by_class["metabolic enzyme"]["mean_z"])
                if "secreted plasma protein" in by_class and "metabolic enzyme" in by_class
                else None,
            "observed_sign_partition": partition,
            "panel_mean_z": float(np.mean(allz)),
            "panel_coherence": coherence(allz),
            # each class scored on its own, because the question the split raises
            # is whether the halves are incoherent or merely opposed
            "class_coherence": {c: coherence([d[g]["z"] for g in gs if g in d])
                                for c, gs in CLASS.items()
                                if any(g in d for g in gs)},
            "panel_t_over_members": {"t": float(t.statistic), "p": float(t.pvalue),
                                     "mean": float(np.mean(allz)),
                                     "ci95": member_ci,
                                     "note": "members are not independent; this "
                                             "interval is descriptive"},
            "panel_composite_sample_test":
                composite_sample_test(mat, tc, pc, D1, restrict=list(d)),
            "panel_composite_variants": {
                f"{sc}_{'min_pairs' if r else 'all_present'}":
                    composite_sample_test(mat, tc, pc, D1,
                                          restrict=list(d) if r else None, scale=sc)
                for sc in ("paired", "all_columns") for r in (True, False)},
            "erythrocyte_mean_z": float(np.mean([v["z"] for v in ery.values()]))
            if ery else None,
            "erythrocyte_members": {g: {"z": round(v["z"], 4),
                                        "n_pairs": v["n_pairs"]} for g, v in ery.items()},
            "immunoglobulin_mean_z": float(np.mean([v["z"] for v in igg.values()]))
            if igg else None,
            "immunoglobulin_members": {g: {"z": round(v["z"], 4),
                                           "n_pairs": v["n_pairs"]} for g, v in igg.items()},
            "off_panel_plasma_mean_z": float(np.mean([v["z"] for v in off.values()]))
            if off else None,
            "off_panel_plasma_members": {g: {"z": round(v["z"], 4),
                                             "n_pairs": v["n_pairs"]} for g, v in off.items()},
            "missingness": missingness(mat, tc, pc),
            "matrix_scale": diag,
            "abundance": {
                "interpretable": diag["on_an_abundance_scale"],
                "spearman_shift_vs_abundance": float(rho.statistic),
                "p": float(rho.pvalue),
                "n_proteins": int(good.sum()),
                "secreted_percentile_median":
                    float(np.median([pct[g] for g in CLASS["secreted plasma protein"]
                                     if g in pct.index])),
                "metabolic_percentile_median":
                    float(np.median([pct[g] for g in CLASS["metabolic enzyme"]
                                     if g in pct.index])),
                "neighborhood_test": band},
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
            "panel_coherence": coherence([v["z"] for v in d.values()]),
            "separation_complete": bool(min(sec) > max(met)),
            "members": {g: round(v["z"], 4) for g, v in d.items()},
            # The composite test belongs inside this loop. Leaving it outside is
            # how a preprocessing-dependent verdict gets reported as a fact.
            "composite_sample_test":
                composite_sample_test(m, tc, pc, D1, restrict=list(d)),
            "composite_sample_test_unfiltered":
                composite_sample_test(m, tc, pc, D1),
            # The class means over every member present, not only those clearing
            # the 10-pair threshold, so that "the split survives all six
            # pipelines" is a statement the result file supports rather than one
            # the reader has to take on trust.
            "unfiltered_member_set": unfiltered_classes(m, tc, pc)}

    # Jiang is 41% unquantified and the two arms do not detect equally, so the
    # split is recomputed on the members that need no complete-case rule at all:
    # those quantified in both arms of every pair.
    mj, tcj, pcj = load_jiang("raw")
    dj = member_shifts(mj, tcj, pcj, D1)
    full = [g for g, v in dj.items() if v["n_pairs"] == len(tcj)]
    secf = [dj[g]["z"] for g in CLASS["secreted plasma protein"] if g in full]
    metf = [dj[g]["z"] for g in CLASS["metabolic enzyme"] if g in full]
    uf = stats.mannwhitneyu(secf, metf, alternative="two-sided") if (
        len(secf) > 2 and len(metf) > 2) else None
    res["jiang_complete_quantification"] = {
        "note": "members quantified in both arms of all pairs, so no member is "
                "scored on a subset of patients",
        "n_pairs": len(tcj),
        "n_members": len(full), "members": sorted(full),
        "members_dropped": sorted(set(dj) - set(full)),
        "secreted_mean_z": float(np.mean(secf)) if secf else None,
        "metabolic_mean_z": float(np.mean(metf)) if metf else None,
        "n_secreted": len(secf), "n_metabolic": len(metf),
        "separation_complete": bool(min(secf) > max(metf)) if secf and metf else None,
        "mannwhitney_p": float(uf.pvalue) if uf else None,
        "panel_coherence": coherence([dj[g]["z"] for g in full]),
        "composite_sample_test": composite_sample_test(mj, tcj, pcj, D1, restrict=full),
        "per_member_detection": {
            g: {"n_pairs_both_arms": dj[g]["n_pairs"],
                "n_tumor": int((~np.isnan(mj.loc[g, tcj].astype(float).values)).sum())
                if g in mj.index else None,
                "n_adjacent": int((~np.isnan(mj.loc[g, pcj].astype(float).values)).sum())
                if g in mj.index else None}
            for g in sorted(dj)},
    }

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
    # Where the panel sits among the sets that were scored at both levels, and
    # how its drop in surviving share compares with theirs. This is the scope
    # answer: whether a fall this large is ordinary or particular to this panel.
    panel_pct = {}
    for n in tx + pr:
        c = res["datasets"][n]["gene_set_coherence"]
        vals = [c[s2]["cancellation"] for s2 in sorted(shared)
                if c[s2]["cancellation"] is not None]
        pv = res["datasets"][n]["panel_coherence"]["cancellation"]
        panel_pct[n] = {
            "panel_surviving_share": float(pv),
            "percentile_within_shared_sets":
                float(100.0 * sum(1 for v in vals if v < pv) / len(vals)),
            "n_shared_sets_compared": len(vals)}
    panel_tx = float(np.mean([panel_pct[n]["panel_surviving_share"] for n in tx]))
    panel_pr = float(np.mean([panel_pct[n]["panel_surviving_share"] for n in pr]))
    drops = np.array([r["transcriptome_cancellation"] - r["proteome_cancellation"]
                      for r in rows])
    panel_drop = panel_tx - panel_pr
    # The panel's surviving share at transcript level is computed over 22
    # members and at protein level over the 18 the proteomes quantify. Comparing
    # those two directly compares two different sets, so the same drop is given
    # again over the members both levels share.
    shared_members = sorted(set.intersection(
        *[set(res["datasets"][n]["panel_members"]) for n in tx + pr]))
    def share_on(nm, members):
        z = np.array([res["datasets"][nm]["panel_members"][g]["z"]
                      for g in members if g in res["datasets"][nm]["panel_members"]])
        return float(abs(z.mean()) / np.abs(z).mean())
    m_tx = {n: share_on(n, shared_members) for n in tx}
    m_pr = {n: share_on(n, shared_members) for n in pr}
    res["panel_on_shared_members"] = {
        "members": shared_members, "n_members": len(shared_members),
        "surviving_share_transcriptome": m_tx,
        "surviving_share_proteome": m_pr,
        "transcriptome_mean": float(np.mean(list(m_tx.values()))),
        "proteome_mean": float(np.mean(list(m_pr.values()))),
        "drop": float(np.mean(list(m_tx.values())) - np.mean(list(m_pr.values())))}

    # where each class, scored alone, would sit among the comparison sets
    class_pct = {}
    for n in tx + pr:
        cc = res["datasets"][n]["gene_set_coherence"]
        vals = [cc[s2]["cancellation"] for s2 in sorted(shared)
                if cc[s2]["cancellation"] is not None]
        class_pct[n] = {}
        for c, co in res["datasets"][n]["class_coherence"].items():
            if co and co["cancellation"] is not None:
                class_pct[n][c] = {
                    "surviving_share": co["cancellation"],
                    "mean_z": co["mean_z"],
                    "n": co["n"],
                    "percentile_within_shared_sets":
                        float(100.0 * sum(1 for v in vals if v < co["cancellation"])
                              / len(vals))}
    res["class_within_shared_sets"] = class_pct
    res["class_sign_across_levels"] = {
        c: {"transcriptome_mean_z": [res["datasets"][n]["class_coherence"][c]["mean_z"]
                                     for n in tx if c in res["datasets"][n]["class_coherence"]],
            "proteome_mean_z": [res["datasets"][n]["class_coherence"][c]["mean_z"]
                                for n in pr if c in res["datasets"][n]["class_coherence"]],
            "keeps_its_sign_in_all_five": bool(len({
                int(np.sign(res["datasets"][n]["class_coherence"][c]["mean_z"]))
                for n in tx + pr if c in res["datasets"][n]["class_coherence"]}) == 1)}
        for c in CLASS}

    res["panel_within_shared_sets"] = {
        "per_dataset": panel_pct,
        "panel_surviving_share_transcriptome_mean": panel_tx,
        "panel_surviving_share_proteome_mean": panel_pr,
        "panel_drop": panel_drop,
        "set_drop_median": float(np.median(drops)),
        "set_drop_max": float(drops.max()),
        "n_sets_dropping_at_least_as_much": int((drops >= panel_drop).sum()),
        "n_sets_compared": int(drops.size),
        "n_sets_with_proteome_share_below_panel":
            int(sum(1 for r in rows if r["proteome_cancellation"] < panel_pr)),
        # the same count taken inside each proteome rather than on the average,
        # which is the number the percentiles in the text describe
        "n_sets_below_panel_within_each_proteome": {
            n: int(sum(1 for s2 in sorted(shared)
                       if res["datasets"][n]["gene_set_coherence"][s2]["cancellation"]
                       is not None
                       and res["datasets"][n]["gene_set_coherence"][s2]["cancellation"]
                       < panel_pct[n]["panel_surviving_share"]))
            for n in pr},
    }

    res["coherence_by_level"] = {
        "n_sets_scored_at_both_levels": len(rows),
        "cancellation_median_transcriptome": float(np.median(a)),
        "cancellation_median_proteome": float(np.median(b)),
        "cancellation_wilcoxon_p": float(stats.wilcoxon(a, b).pvalue),
        "n_sets_cancellation_lower_in_proteome": int((b < a).sum()),
        "sign_agreement_median_transcriptome": float(np.median(sa)),
        "sign_agreement_median_proteome": float(np.median(sb)),
        "sign_agreement_wilcoxon_p": float(stats.wilcoxon(sa, sb).pvalue),
        # "does not change" is an acceptance of a null unless the interval is
        # given, so it is given.
        "cancellation_mean_difference": float((a - b).mean()),
        "cancellation_difference_ci95": [
            float((a - b).mean() - stats.t.ppf(0.975, a.size - 1) * (a - b).std(ddof=1) / np.sqrt(a.size)),
            float((a - b).mean() + stats.t.ppf(0.975, a.size - 1) * (a - b).std(ddof=1) / np.sqrt(a.size))],
        "sign_agreement_mean_difference": float((sa - sb).mean()),
        "sign_agreement_difference_ci95": [
            float((sa - sb).mean() - stats.t.ppf(0.975, sa.size - 1) * (sa - sb).std(ddof=1) / np.sqrt(sa.size)),
            float((sa - sb).mean() + stats.t.ppf(0.975, sa.size - 1) * (sa - sb).std(ddof=1) / np.sqrt(sa.size))],
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
