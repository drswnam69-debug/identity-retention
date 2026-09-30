#!/usr/bin/env python3
"""80_verify_MCP.py

Two checks on manuscript_MCP.md.

A. EXHAUSTIVE NUMERIC MATCHING. Every number printed in the abstract, the In
   Brief, the body and the figure legends must be reproducible from
   results/MANUSCRIPT_NUMBERS.json. A checker that looks only for sentences I
   thought to look for finds only what it was told to look for, so the default
   is that an unmatched number is an error and has to be either traced or
   whitelisted with a reason.

B. TARGETED CLAIMS. Statements whose error would not show up as a wrong
   number, for example a direction, an ordering or a count of members.

Exit code is non-zero if anything fails.
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MS = f"{ROOT}/manuscript/manuscript_MCP.md"
NUM = f"{ROOT}/results/MANUSCRIPT_NUMBERS.json"

# numbers that are not results: structural constants, identifiers, thresholds
WHITELIST = {
    # structural / definitional
    "22": "the panel has 22 genes, by definition",
    "10": "the 10-pair threshold and the ten plasma proteins",
    "6": "six metabolic enzymes in the panel definition",
    "5": "five transcription factors in the panel definition",
    "1": "one surface receptor; and 'a surviving share of 1'",
    "4": "four classes; four enzymes; four limits",
    "3": "three transcriptomes; Python 3; at least three gate members",
    "2": "two levels, two classes, two blood compartments, two cohorts",
    "8": "eight datasets",
    "0": "zero, as a reference value",
    "149": "sets retrieved by name before filtering",
    "126": "sets kept at 15 to 1,000 genes",
    "15": "lower size bound for eligibility; 15 proteins in the replicate check",
    "1000": "upper size bound for eligibility",
    "1.5": "the abundance neighborhood half-width in percentile points",
    "0.9": "the pre-registered surviving-share threshold",
    "0.5": "the within-0.5-of-zero diagnostic",
    "20": "the 20% below-zero diagnostic threshold",
    "0.05": "not used; guard",
    "445": "columns in the GSE14520 GPL3921 series matrix",
    "05": "date fragment of the Xena file version",
    "09": "date fragment of the Xena file version",
    "2024": "date fragment of the Xena file version",
    "0005615": "GO term for extracellular space",
    "0005576": "GO term for extracellular region",
    "0070062": "GO term for extracellular exosome",
    "38": "GRCh38",
    "13": "GRCh38.p13",
    "2026": "MSigDB release year",
    "9477": "size of the accession map",
    "764": "total pairs; checked separately against the sum",
    "2019": "part of the dataset names Gao 2019 and Jiang 2019",
    "2023": "part of the dataset name Yi 2023",
    "000095336": "the MassIVE accession of the cholangiocarcinoma deposit",
}
SKIP_SECTIONS = ("## References", "## Abbreviations")


def body_text(s):
    """the parts whose numbers must all be traceable"""
    out = []
    for name, nxt in [("Abstract", "In Brief"), ("In Brief", "Introduction"),
                      ("Introduction", "Experimental Procedures"),
                      ("Experimental Procedures", "Results"),
                      ("Results", "Discussion"),
                      ("Discussion", "Acknowledgments"),
                      ("Figure legends", None)]:
        i = s.index(f"## {name}")
        j = s.index(f"## {nxt}") if nxt else len(s)
        out.append(s[i:j])
    return "\n".join(out)


def producible(o, acc):
    """every leaf number in the JSON, at several roundings"""
    if isinstance(o, dict):
        for v in o.values():
            producible(v, acc)
    elif isinstance(o, list):
        for v in o:
            producible(v, acc)
    elif isinstance(o, bool):
        pass
    elif isinstance(o, (int, float)):
        x = float(o)
        for d in range(0, 5):
            acc.add(f"{abs(x):.{d}f}".rstrip("0").rstrip(".") if d else f"{abs(x):.0f}")
            acc.add(f"{abs(x):.{d}f}")
        for d in range(0, 3):                       # as a percentage
            acc.add(f"{abs(x) * 100:.{d}f}")
        for d in range(0, 3):                       # as a power of ten exponent
            pass
    return acc


def main() -> int:
    s = open(MS, encoding="utf-8").read()
    s_ms = s
    d = json.load(open(NUM))
    ok = set()
    producible(d, ok)
    # values the text may state as a derived difference or sum
    rows = d["datasets"]
    ok.add(f"{sum(r['n_pairs'] for r in rows)}")
    for g in d["gap_secreted_minus_metabolic"].values():
        for dd in range(0, 4):
            ok.add(f"{abs(g):.{dd}f}")
    ok = {t.rstrip("0").rstrip(".") if "." in t else t for t in ok} | ok

    txt = body_text(s)
    txt = re.sub(r"`[^`]*`", " ", txt)              # code spans: file names
    txt = re.sub(r"\[[0-9,\s]+\]", " ", txt)        # citation brackets
    txt = re.sub(r"\*([A-Z0-9]+)\*", " ", txt)      # italicized gene symbols
    txt = re.sub(r"GO:\d+", " ", txt)
    txt = re.sub(r"10\^-?\d+\^", " TENPOW ", txt)
    txt = re.sub(r"\b(GRCh38\.p13|PXD\d+|GSE\d+|MSV\d+|v2026\.1\.Hs)\b", " ", txt)

    tokens = re.findall(r"(?<![A-Za-z0-9.])(\d+(?:,\d{3})*(?:\.\d+)?)", txt)
    fails, used_wl = [], set()
    for t in tokens:
        raw = t.replace(",", "")
        cands = {raw, raw.rstrip("0").rstrip(".") if "." in raw else raw}
        if cands & ok:
            continue
        if raw in WHITELIST:
            used_wl.add(raw); continue
        fails.append(t)

    print("== A. exhaustive numeric matching ==")
    print(f"   numeric tokens examined: {len(tokens)}")
    print(f"   whitelisted constants used: {len(used_wl)}")
    if fails:
        print(f"   UNTRACED: {len(fails)}")
        for t in sorted(set(fails), key=lambda x: -len(x)):
            m = re.search(r".{0,80}\b" + re.escape(t) + r"\b.{0,60}", txt, re.S)
            print(f"     {t:>12s}  ...{(m.group(0) if m else '').strip()[:130]}...")
    else:
        print("   every number traces to the result file")

    # ---------------------------------------------------------------- claims
    print("\n== B. targeted claims ==")
    checks = []
    by = {r["dataset"]: r for r in rows}

    def add(name, cond):
        checks.append((name, bool(cond)))

    add("eight datasets are tabulated", d["n_datasets"] == 8)
    add("total pairs is 764", sum(r["n_pairs"] for r in rows) == 764)
    add("the module is negative in all eight",
        all(r["metabolic_mean_z"] < 0 for r in rows))
    add("the module's surviving share is at least 0.95 in all eight",
        min(r["metabolic_share"] for r in rows) >= 0.95)
    add("the text does not round that minimum up to 'at least 0.954'",
        "at least 0.954" not in s_ms)
    add("plasma is negative in all three transcriptomes",
        all(r["secreted_mean_z"] < 0 for r in rows if r["level"] == "transcriptome"))
    add("plasma is positive in all five proteomes",
        all(r["secreted_mean_z"] > 0 for r in rows if r["level"] == "proteome"))
    add("plasma strength in proteomes runs from +0.076 to +0.514",
        abs(min(r["secreted_mean_z"] for r in rows if r["level"] == "proteome") - 0.076) < 5e-4
        and abs(max(r["secreted_mean_z"] for r in rows if r["level"] == "proteome") - 0.514) < 5e-4)
    add("plasma surviving share in proteomes runs from 0.215 to 1.000",
        abs(min(r["secreted_share"] for r in rows if r["level"] == "proteome") - 0.215) < 5e-3
        and abs(max(r["secreted_share"] for r in rows if r["level"] == "proteome") - 1.0) < 5e-4)

    cs = d["comparison_sets"]
    pdrop = cs["panel_drop"]
    add("no comparison set drops as much as the panel",
        pdrop["n_sets_dropping_at_least_as_much"] == 0)
    add("the largest set drop is smaller than the panel's",
        pdrop["set_drop_max"] < pdrop["panel_drop"])
    add("the 114-set distribution does not move (CI spans zero)",
        cs["difference_ci95"][0] < 0 < cs["difference_ci95"][1])
    add("the panel is above the 50th percentile in every transcriptome",
        all(v["percentile_within_shared_sets"] > 50
            for k, v in pdrop["per_dataset"].items() if "proteome" not in k))
    add("the panel is below the 25th percentile in both first proteomes",
        all(v["percentile_within_shared_sets"] < 25
            for k, v in pdrop["per_dataset"].items() if "proteome" in k))

    sq = d["secretome_question"]
    add("secreted sit below the rest at gene level in all five",
        all(v["difference_of_means"] < 0
            for v in sq["gene_level_secreted_minus_other"].values()))
    add("the secretome-fraction predictor is null",
        sq["h1_drop_vs_secreted_fraction"]["spearman_p"] > 0.05)
    add("the gene-disjoint subfamily has 2 sets",
        sq["n_sets_in_gene_disjoint_subfamily"] == 2)
    add("more than half of set pairs share a gene",
        sq["set_overlap"]["share_of_pairs_sharing_any_gene"] > 0.5)
    add("the abundance-matched difference is null in Jiang only",
        sq["abundance_matched"]["Jiang 2019 proteome"]["wilcoxon_p"] > 0.5
        and all(v["wilcoxon_p"] < 1e-30 for k, v in sq["abundance_matched"].items()
                if v.get("level") == "transcriptome"))

    h = d["held_out_test"]["EncyclopeDIA"]
    add("the held-out gate passes", h["gate"]["passes"])
    add("held-out H1 passes", h["H1_metabolic_module"]["passes"])
    add("held-out H2 passes", h["H2_plasma_proteins"]["passes"])
    add("held-out H3 passes", h["H3_extremity"]["passes"])
    add("held-out composite is not significant",
        h["H4_composite"]["composite"]["wilcoxon_p"] > 0.05)
    add("TAT is the only module member that does not fall in the held-out cohort",
        [g for g, z in h["H1_metabolic_module"]["members_quantified"].items() if z > 0] == ["TAT"])
    add("abundance accounts for less than a fifth of the module's fall",
        h["H5_abundance"]["share_of_module_fall_explained_by_abundance"] < 0.2)
    add("the primary matrix is the most complete of the three",
        d["held_out_primary_sheet"] == "EncyclopeDIA")
    add("DIA-NN disagrees on H1", not d["held_out_test"]["DIANN"]["H1_metabolic_module"]["passes"])
    add("Spectronaut cannot be evaluated",
        not d["held_out_test"]["Spectronaut"]["gate"]["passes"])

    for label, c in d["cholangiocarcinoma"]["cohorts"].items():
        add(f"{label}: pre-registered gate failed", not c["gate_preregistered"]["passes"])
        add(f"{label}: amended gate passed", c["gate_amended"]["passes"])
        add(f"{label}: module falls with share 1.000",
            c["G1_module"]["passes"] and c["G1_module"]["coherence"]["cancellation"] == 1.0)
        add(f"{label}: module deeper than every HCC cohort",
            c["G2_deeper_than_hcc"]["deeper_than_all_three"])
        add(f"{label}: plasma rises, so the prediction failed",
            not c["G3_plasma"]["passes"])
        add(f"{label}: module below the 10th percentile",
            c["G5_extremity"]["passes"])
        add(f"{label}: albumin rises", c["G3_plasma"]["members"]["ALB"] > 0)
        add(f"{label}: the four coagulation proteins fall",
            all(c["G3_plasma"]["members"][g] < 0 for g in ("FGA", "FGB", "FGG", "F2")))
        add(f"{label}: exploratory status recorded", "EXPLORATORY" in c["status"])

    ann = sq["annotation_check"]["strict_GO_0005615"]
    add("the annotation calls 20 of 22 plasma proteins secreted",
        ann["known_plasma_proteins_called_secreted"] == "20/22")
    add("the annotation's only intracellular false positive is ARG1",
        ann["false_positives"] == ["ARG1"])

    bad = [n for n, v in checks if not v]
    for n, v in checks:
        print(f"   {'PASS' if v else 'FAIL'}  {n}")
    print(f"\n{len(checks) - len(bad)}/{len(checks)} claim checks pass")
    if fails or bad:
        print("\nFAILED")
        return 1
    print("\nALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
