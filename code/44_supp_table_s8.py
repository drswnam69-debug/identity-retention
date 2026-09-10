#!/usr/bin/env python3
"""Supplementary Table S8: every enumerated gene set in the three tissues,
with its size, platform coverage, unadjusted paired effect, inclusion status
and the reason for any exclusion.

Reads only archived artifacts; computes nothing new.
"""
import csv, json, os, re

ROOT = "/home/claude/rsi"
GS   = os.path.join(ROOT, "genesets")
RES  = os.path.join(ROOT, "results")
OUT  = "/home/claude/SupplementaryTable_S8_enumeration.csv"

GMT  = os.path.join(GS, "c2.cgp.v2026.1.Hs.symbols.gmt")
KEYS = {"liver":  ("LIVER", "HEPAT"),
        "lung":   ("LUNG", "PULMONARY", "ALVEOLAR"),
        "kidney": ("KIDNEY", "RENAL")}
COHORT = {"liver": "GSE14520", "lung": "TCGA_LUAD", "kidney": "TCGA_KIRC"}
SIZE_MIN, SIZE_MAX, MIN_ON_PLATFORM = 15, 1000, 10

# full collection
sets = {}
with open(GMT) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) > 2:
            sets[p[0]] = sorted(set(g for g in p[2:] if g))

rows = []
for tissue, keys in KEYS.items():
    bench = json.load(open(os.path.join(RES, f"SIGNATURE_BENCHMARK_{COHORT[tissue]}.json")))
    byname = {s["name"]: s for s in bench["signatures"]}
    names = sorted(n for n in sets if any(k in n for k in keys))
    for n in names:
        size = len(sets[n])
        r = {"tissue": tissue, "cohort": COHORT[tissue], "set_name": n,
             "n_genes_in_set": size, "n_genes_on_platform": "",
             "unadjusted_mean_delta": "", "wilcoxon_p": "",
             "retention_D1": "", "retention_joint": "",
             "status": "", "exclusion_reason": ""}
        if not (SIZE_MIN <= size <= SIZE_MAX):
            r["status"] = "excluded"
            r["exclusion_reason"] = (f"set size {size} outside {SIZE_MIN} to {SIZE_MAX}")
            rows.append(r); continue
        s = byname.get(n)
        if s is None:
            r["status"] = "excluded"
            r["exclusion_reason"] = f"fewer than {MIN_ON_PLATFORM} symbols present on the platform"
            rows.append(r); continue
        r["n_genes_on_platform"] = s.get("n_on_platform", "")
        r["unadjusted_mean_delta"] = s.get("unadjusted_mean_delta", "")
        r["wilcoxon_p"] = s.get("wilcoxon_p", "")
        st = s.get("status", "")
        if st == "ok":
            r["status"] = "included"
            r["retention_D1"] = s.get("retention_D1", "")
            r["retention_joint"] = s.get("retention_joint", "")
        elif st in ("null_effect_excluded", "null_effect"):
            r["status"] = "excluded"
            r["exclusion_reason"] = "no unadjusted paired effect (two-sided Wilcoxon P >= 0.05)"
        elif st in ("not_on_platform", "not_evaluable_on_platform"):
            r["status"] = "excluded"
            r["exclusion_reason"] = f"fewer than {MIN_ON_PLATFORM} symbols present on the platform"
        else:
            r["status"] = "excluded"
            r["exclusion_reason"] = st
        rows.append(r)

cols = ["tissue", "cohort", "set_name", "n_genes_in_set", "n_genes_on_platform",
        "unadjusted_mean_delta", "wilcoxon_p", "retention_D1", "retention_joint",
        "status", "exclusion_reason"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)

from collections import Counter
print(OUT, len(rows), "rows")
for t in KEYS:
    sub = [r for r in rows if r["tissue"] == t]
    print(f"  {t:7s} enumerated {len(sub):4d}  included {sum(r['status']=='included' for r in sub):4d}")
print(Counter(r["exclusion_reason"] for r in rows if r["status"] == "excluded"))
