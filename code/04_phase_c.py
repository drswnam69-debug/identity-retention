#!/usr/bin/env python3
"""04_phase_c.py -- H2 and H3: the HCC transition and the field effect.

H2  RSI is higher in NASH-HCC tumor than in adjacent non-tumor liver.
H3  RSI in adjacent non-tumor NASH liver is already higher than in NASH liver
    from patients WITHOUT HCC. This is the field-effect test and the most
    original comparison in the study.

Tumor status is an exposure here, not a survival outcome; PHASE E outcomes
stay unopened.

Example
-------
python3 code/04_phase_c.py --cohort GSE164760 --outdir results
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rsi_config as cfg      # noqa: E402
import stats_lite as sl       # noqa: E402

# GEO tissue labels, in increasing order of proximity to cancer
TISSUE_ORDER = ["Healthy liver", "NASH liver", "Cirrhotic liver",
                "Non-tumoral NASH liver adjacent to HCC", "NASH-HCC tumor"]
SHORT = {"Healthy liver": "Healthy", "NASH liver": "NASH",
         "Cirrhotic liver": "Cirrhotic",
         "Non-tumoral NASH liver adjacent to HCC": "Adjacent",
         "NASH-HCC tumor": "Tumor"}


def hodges_lehmann(x, y) -> float:
    """Median of all pairwise differences: the shift Mann-Whitney tests."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    if len(x) * len(y) > 4_000_000:          # subsample huge products
        rng = np.random.default_rng(0)
        x = rng.choice(x, 2000, replace=False)
        y = rng.choice(y, 2000, replace=False)
    return float(np.median(x[:, None] - y[None, :]))


def cliffs_delta(x, y) -> float:
    """Nonparametric effect size in [-1, 1]; 0 means complete overlap."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return float("nan")
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float((gt - lt) / (len(x) * len(y)))


def contrast(a, b, label_a: str, label_b: str) -> dict:
    u, z, p = sl.mannwhitney_u(a, b)
    return {
        "comparison": f"{label_a} vs {label_b}",
        "n_a": int(len(a)), "n_b": int(len(b)),
        "median_a": round(float(np.median(a)), 4) if len(a) else None,
        "median_b": round(float(np.median(b)), 4) if len(b) else None,
        "hodges_lehmann_shift": round(hodges_lehmann(a, b), 4),
        "cliffs_delta": round(cliffs_delta(a, b), 4),
        "mannwhitney_z": round(z, 4), "p_two_sided": p,
    }


def signed_rank(diff: np.ndarray) -> tuple[float, float]:
    """Wilcoxon signed-rank on paired differences, normal approximation."""
    d = np.asarray([v for v in diff if np.isfinite(v) and v != 0], float)
    n = len(d)
    if n < 6:
        return float("nan"), float("nan")
    r = sl.rankdata(np.abs(d))
    w = float(r[d > 0].sum())
    mu = n * (n + 1) / 4.0
    _, counts = np.unique(np.abs(d), return_counts=True)
    tie = ((counts ** 3 - counts).sum()) / 48.0
    var = n * (n + 1) * (2 * n + 1) / 24.0 - tie
    if var <= 0:
        return float("nan"), float("nan")
    z = (w - mu) / math.sqrt(var)
    return z, sl.norm_two_sided(z)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default="GSE164760")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()

    out = os.path.join(args.outdir, args.cohort)
    scores = pd.read_csv(os.path.join(out, "rsi.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(out, "phenotype.tsv"), sep="\t",
                        index_col=0).loc[scores.index]
    tissue = pheno["tissue"].astype(str).str.strip()

    def g(label, col="RSI"):
        return scores.loc[tissue == label, col].to_numpy(float)

    res: dict = {
        "cohort": args.cohort, "n": int(len(scores)),
        "lock_hash": cfg.lock_hash(),
        "platform_note": "Affymetrix HG-U219 microarray; analyzed as an "
                         "independent arm, never merged with the RNA-seq "
                         "cohorts",
        "group_sizes": {SHORT.get(t, t): int((tissue == t).sum())
                        for t in TISSUE_ORDER},
    }

    # ---- H2: tumor vs adjacent non-tumor -------------------------------
    res["H2"] = contrast(g("NASH-HCC tumor"),
                         g("Non-tumoral NASH liver adjacent to HCC"),
                         "tumor", "adjacent")
    res["H2"]["hypothesis"] = cfg.HYPOTHESES["H2"]

    # ---- H3: the field effect ------------------------------------------
    res["H3"] = contrast(g("Non-tumoral NASH liver adjacent to HCC"),
                         g("NASH liver"), "adjacent", "NASH without HCC")
    res["H3"]["hypothesis"] = cfg.HYPOTHESES["H3"]

    # ---- context contrasts ---------------------------------------------
    res["context"] = [
        contrast(g("NASH liver"), g("Healthy liver"), "NASH", "healthy"),
        contrast(g("Cirrhotic liver"), g("NASH liver"), "cirrhotic", "NASH"),
        contrast(g("NASH-HCC tumor"), g("NASH liver"), "tumor",
                 "NASH without HCC"),
    ]

    # ---- ordered trend across the whole ladder --------------------------
    groups = [g(t) for t in TISSUE_ORDER if (tissue == t).any()]
    labels = [SHORT[t] for t in TISSUE_ORDER if (tissue == t).any()]
    jt, z, p = sl.jonckheere_terpstra(groups, alternative="increasing")
    res["ordered_trend"] = {
        "ladder": labels, "n_per_group": [len(x) for x in groups],
        "median_per_group": [round(float(np.median(x)), 4) for x in groups],
        "JT": round(jt, 1), "z": round(z, 4), "p_increasing": p,
        "note": "exploratory: the ladder mixes tissue types, it is not a "
                "within-patient progression",
    }

    # ---- components and genes for the two primary contrasts -------------
    res["components"] = {}
    for comp in ("z_supply", "z_reduction", "z_drain"):
        res["components"][comp] = {
            "H2": contrast(g("NASH-HCC tumor", comp),
                           g("Non-tumoral NASH liver adjacent to HCC", comp),
                           "tumor", "adjacent"),
            "H3": contrast(g("Non-tumoral NASH liver adjacent to HCC", comp),
                           g("NASH liver", comp), "adjacent", "NASH"),
        }

    gene_cols = [c for c in scores.columns if c.startswith("gene_")]
    genes = {}
    for c in gene_cols:
        genes[c[5:]] = {
            "H2": contrast(g("NASH-HCC tumor", c),
                           g("Non-tumoral NASH liver adjacent to HCC", c),
                           "tumor", "adjacent"),
            "H3": contrast(g("Non-tumoral NASH liver adjacent to HCC", c),
                           g("NASH liver", c), "adjacent", "NASH"),
        }
    for key in ("H2", "H3"):
        names = list(genes)
        q = sl.benjamini_hochberg([genes[n][key]["p_two_sided"] for n in names])
        for n, qv in zip(names, q):
            genes[n][key]["q_bh"] = float(qv)
    res["genes"] = genes

    # ---- sensitivity: paired H2 on inferred patient pairs ---------------
    # GEO declares no patient identifier. Source names suggest consecutive
    # numbering for tumor / adjacent from one patient. That is an inference,
    # so this is a sensitivity analysis, never the primary test.
    pairs = []
    src_col = next((c for c in ("source_name_ch1", "source_name", "title")
                    if c in pheno.columns), None)
    if src_col is not None:
        src = pheno[src_col].astype(str)
        num = src.str.extract(r"(\d+)")[0]
        tum = {int(n): i for i, (n, t) in enumerate(zip(num, tissue))
               if pd.notna(n) and t == "NASH-HCC tumor"}
        adj = {int(n): i for i, (n, t) in enumerate(zip(num, tissue))
               if pd.notna(n) and t.startswith("Non-tumoral")}
        for n, i in sorted(tum.items()):
            if n + 1 in adj:
                pairs.append((scores["RSI"].iloc[i],
                              scores["RSI"].iloc[adj[n + 1]]))
    if len(pairs) >= 6:
        d = np.array([a - b for a, b in pairs], float)
        z_sr, p_sr = signed_rank(d)
        res["H2_paired_sensitivity"] = {
            "n_pairs": len(pairs),
            "median_difference": round(float(np.median(d)), 4),
            "signed_rank_z": round(z_sr, 4), "p_two_sided": p_sr,
            "caveat": "pairs inferred from consecutive source identifiers; "
                      "GEO declares no patient ID. Sensitivity only.",
        }
    else:
        res["H2_paired_sensitivity"] = {
            "n_pairs": len(pairs), "note": "too few inferred pairs to test"}

    with open(os.path.join(out, "phase_c.json"), "w") as fh:
        json.dump(res, fh, indent=2, default=float)

    # ---- console summary ------------------------------------------------
    print(f"\n=== {args.cohort} — PHASE C ===")
    print("  groups: " + ", ".join(f"{k}={v}" for k, v in
                                   res["group_sizes"].items()))
    for key, title in (("H2", "H2  tumor vs adjacent non-tumor"),
                       ("H3", "H3  adjacent vs NASH without HCC  [field effect]")):
        r = res[key]
        print(f"\n  {title}")
        print(f"    medians {r['median_a']} vs {r['median_b']}   "
              f"shift {r['hodges_lehmann_shift']:+.3f}")
        print(f"    Cliff's delta {r['cliffs_delta']:+.3f}   "
              f"z={r['mannwhitney_z']:.3f}   P={r['p_two_sided']:.3g}")
    ps = res["H2_paired_sensitivity"]
    if "signed_rank_z" in ps:
        print(f"\n  H2 paired sensitivity ({ps['n_pairs']} inferred pairs): "
              f"median diff {ps['median_difference']:+.3f}, "
              f"P={ps['p_two_sided']:.3g}")
    ot = res["ordered_trend"]
    print(f"\n  exploratory ladder {' < '.join(ot['ladder'])}")
    print(f"    medians {ot['median_per_group']}")
    print(f"    JT z={ot['z']:.3f}  P={ot['p_increasing']:.3g}")
    sig = {k: v for k, v in genes.items() if v["H3"].get("q_bh", 1) < 0.05}
    print(f"\n  genes with q<0.05 for the field-effect contrast: "
          f"{len(sig)}/{len(genes)}")
    print(f"  wrote -> {out}/phase_c.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
