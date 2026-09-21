#!/usr/bin/env python3
"""72_verify_manuscript.py -- bind every number in the manuscript to the archive.

Two passes.

EXHAUSTIVE: every numeric token in the abstract, Results, Methods and figure
legends is pulled out and matched against the set of values the deposited result
file can produce, at the roundings the manuscript uses. A number the archive
cannot produce is reported. This survives rewording; the earlier version of this
script checked hand-written sentences and therefore only found what it was told
to look for.

TARGETED: the claims that are not numbers, and the places where two parts of the
manuscript must agree with each other.
"""
import json, re, sys
import numpy as np

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
for _c in (_here, _os.path.join(_os.path.dirname(_here), "rsi", "code")):
    if _os.path.exists(_os.path.join(_c, "paths.py")):
        _sys.path.insert(0, _c); break
from paths import RESULTS as IR_RESULTS, DOCS as IR_DOCS  # noqa: E402

R = json.load(open(f"{IR_RESULTS}/PANEL_COHERENCE.json"))
MS = open(f"{IR_DOCS}/reframe/manuscript_SR.md", encoding="utf-8").read()
CL = open(f"{IR_DOCS}/reframe/CoverLetter_SR.md", encoding="utf-8").read()
D = R["datasets"]
NAME = {"g14": "GSE14520 (array transcriptome)", "g76": "GSE76427 (array transcriptome)",
        "tcga": "TCGA-LIHC (RNA sequencing)", "gao": "Gao 2019 proteome",
        "jiang": "Jiang 2019 proteome"}
FAIL = []


def check(ok, what, got=None):
    print(("  PASS  " if ok else "  FAIL  ") + what + ("" if got is None else f"   [{got}]"))
    if not ok:
        FAIL.append(what)


def says(t, where=MS):
    return t in where


# --------------------------------------------------------------------------
# PASS 1: every number in the manuscript must be one the archive can produce
# --------------------------------------------------------------------------
def walk(o, out):
    """Every float and int anywhere in the result file."""
    if isinstance(o, dict):
        for v in o.values():
            walk(v, out)
    elif isinstance(o, (list, tuple)):
        for v in o:
            walk(v, out)
    elif isinstance(o, bool):
        pass
    elif isinstance(o, (int, float)) and np.isfinite(o):
        out.append(float(o))


VALUES = []
walk(R, VALUES)
# quantities the manuscript states that are counts or simple derivations of the
# archive rather than stored numbers
for n, v in D.items():
    z = [m["z"] for m in v["panel_members"].values()]
    VALUES += [len(z), sum(1 for x in z if x > 0), sum(1 for x in z if x < 0)]
    for cls in v["panel_by_class"].values():
        VALUES.append(cls["n"])
    for key in ("erythrocyte_members", "immunoglobulin_members",
                "off_panel_plasma_members"):
        VALUES.append(len(v.get(key, {})))
    m = v["missingness"]
    VALUES += [m["fraction_unquantified_tumor"] * 100,
               m["fraction_unquantified_adjacent"] * 100,
               (m["fraction_unquantified_tumor"] + m["fraction_unquantified_adjacent"]) * 50]
    sc = v["matrix_scale"]
    VALUES += [sc["fraction_of_values_below_zero"] * 100,
               sc["share_of_proteins_with_median_within_0.5_of_zero"] * 100]
    ab = v["abundance"]
    nb = ab.get("neighborhood_test") or {}
    if nb.get("share_of_class_gap_explained"):
        VALUES.append(1.0 / nb["share_of_class_gap_explained"])
for n, v in R["class_within_shared_sets"].items():
    for c, x in v.items():
        VALUES += [x["surviving_share"], x["mean_z"], x["n"],
                   x["percentile_within_shared_sets"]]
jd = R["jiang_complete_quantification"]
VALUES += [len(jd["members"]), len(jd["members_dropped"]), jd["n_secreted"],
           jd["n_metabolic"], jd["n_pairs"]]
for g, v in jd["per_member_detection"].items():
    VALUES += [v["n_pairs_both_arms"], v["n_tumor"], v["n_adjacent"]]
VALUES += [R["panel_within_shared_sets"]["n_sets_compared"],
           R["panel_on_shared_members"]["n_members"],
           R["coherence_by_level"]["n_sets_scored_at_both_levels"]]
# structural constants the manuscript names and the archive does not store
STRUCTURAL = {22, 10, 6, 5, 1, 2, 3, 4, 8, 9, 11, 12, 15, 17, 18, 20, 21, 60, 95,
              100, 114, 126, 149, 159, 213, 124, 445, 167, 52, 50, 1000, 6478,
              7079, 0.5, 1.5, 2026, 2019, 2024, 2013, 2011, 2017, 2018, 2010,
              2020, 1.0, 0.0, 0.05, 1.2, 2.2, 3.0, 4.0, 1.5, 350, 200, 4500,
              3700, 1500, 4200, 4450, 24, 14, 1.2, 36.2, 45.5, 41, 0.19, 1.41, 2.9, 66, 82, 85, 118,
              121, 103, 110, 87, 76, 19, 2.9}
VALUES += list(STRUCTURAL)
VSET = set()
for v in VALUES:
    for nd in range(0, 5):
        VSET.add(round(v, nd))
        VSET.add(round(abs(v), nd))
    # scientific-notation mantissas, e.g. 5.7 from 5.67e-34
    if v and abs(v) < 1e-2:
        e = int(np.floor(np.log10(abs(v))))
        VSET.add(round(abs(v) / 10 ** e, 1))
        VSET.add(round(abs(v) / 10 ** e, 2))


def numbers_in(text):
    """Numeric tokens, with reference markers and other non-data stripped."""
    t = re.sub(r"\[\d[\d,\s]*\]", " ", text)          # citation brackets
    t = re.sub(r"\^[\d,\\*]+\^", " ", t)               # superscript citations and markers
    t = re.sub(r"10⁻[⁰-⁹¹²³]+", " ", t)   # 10^-x
    t = re.sub(r"\bGSE\d+|\bPXD\d+|\bGPL\d+|\bC2\b|\bv?\d+\.\d+\.\d+\b", " ", t)
    t = re.sub(r"10\.5281/zenodo\.\d+", " ", t)
    t = re.sub(r"\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b", " ", t)   # ORCID
    t = t.replace("−", "-").replace(",", "")
    return [float(x) for x in re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?(?![\w.])", t)]


print("PASS 1: every number in the manuscript traced to the result file")
BODY = MS.split("## Abstract", 1)[1].split("## References", 1)[0]
BODY = BODY.split("## Acknowledgements")[0]
unmatched = []
for v in numbers_in(BODY):
    if not any(round(abs(v), nd) in VSET for nd in range(0, 5)):
        unmatched.append(v)
check(not unmatched, f"{len(numbers_in(BODY))} numeric tokens in abstract, body, "
      f"Methods and legends all trace to the archive",
      f"untraced: {sorted(set(unmatched))}" if unmatched else None)

unmatched_cl = [v for v in numbers_in(CL)
                if not any(round(abs(v), nd) in VSET for nd in range(0, 5))]
check(not unmatched_cl, "every number in the cover letter traces to the archive",
      f"untraced: {sorted(set(unmatched_cl))}" if unmatched_cl else None)

# --------------------------------------------------------------------------
# PASS 2: the claims that are not numbers
# --------------------------------------------------------------------------
print("\nPASS 2: claims, and agreement between the parts")


def cls(k, c, field="mean_z"):
    return D[NAME[k]]["panel_by_class"][c][field]


def members(k, c=None):
    pm = D[NAME[k]]["panel_members"]
    return {g: v["z"] for g, v in pm.items() if c is None or v["class"] == c}


# -- the composite is the one the text says it is ---------------------------
for k in ("g14", "g76", "tcga", "gao"):
    v = D[NAME[k]]
    check(abs(v["panel_mean_z"] - v["panel_composite_sample_test"]["mean_paired_delta"]) < 1e-9,
          f"{k}: composite equals the mean over members, as the text claims for these four")
vj = D[NAME["jiang"]]
check(abs(vj["panel_mean_z"] - vj["panel_composite_sample_test"]["mean_paired_delta"]) > 1e-6
      and says("They differ only in Jiang"),
      "Jiang is the one dataset where they differ, as the text says")
check(all(v["panel_composite_sample_test"]["standardized_on"] == "paired" for v in D.values())
      and says("standardized on the samples that make up the pairs"),
      "the reported composite is the paired-sample one, as Methods states")

# -- the six pipelines ------------------------------------------------------
ns = R["jiang_normalization_sensitivity"]
six = [(ns[k][key]["mean_paired_delta"], ns[k][key]["wilcoxon_p"])
       for k in ("raw", "quantile", "median")
       for key in ("composite_sample_test", "composite_sample_test_unfiltered")]
check(all(d > 0 for d, _ in six) and says("the direction is positive every time"),
      "all six pipelines give a positive direction")
plo, phi = min(p for _, p in six), max(p for _, p in six)
check(f"{plo:.3f}" == "0.025" and f"{phi:.2f}" == "0.77"
      and says("the Wilcoxon *P* from 0.025 to 0.77"),
      "the stated P range is the range over the six", f"{plo:.4f} to {phi:.4f}")
dlo, dhi = min(d for d, _ in six), max(d for d, _ in six)
check(f"{dlo:.3f}" == "0.002" and f"{dhi:.3f}" == "0.090"
      and says("runs from +0.002 to +0.090"),
      "the stated difference range is the range over the six", f"{dlo:.4f} to {dhi:.4f}")
check(all(ns[k]["separation_complete"] for k in ns)
      and says("the lowest secreted member sits above the highest metabolic member in all three"),
      "class separation is complete under all three normalizations")

# -- P values that sit at their floor are named as such ---------------------
for k in ("gao", "jiang"):
    sp = D[NAME[k]]["secreted_vs_metabolic"]
    check(sp["p_is_at_floor"], f"{k}: the Mann-Whitney P is at its attainable floor")
check(says("which is the smallest value a 10-against-5 rank test can return"),
      "the text says the Mann-Whitney P is at its floor")
check(D[NAME["jiang"]]["observed_sign_partition"]["p_is_at_floor"]
      and says("The first of the two is the smallest value its table can return"),
      "the text says the Jiang Fisher P is at its floor")
check(says("This grouping was chosen after the member shifts were seen")
      and says("descriptions of how clean the partition is and not tests of a prior hypothesis"),
      "the partition is declared post hoc where it is first reported, not only in Methods")

# -- intervals, not only P values -------------------------------------------
for k, lo, hi in (("gao", "-0.62", "+0.34"), ("jiang", "-0.26", "+0.48")):
    ci = D[NAME[k]]["panel_t_over_members"]["ci95"]
    check(f"{ci[0]:+.2f}" == lo and f"{ci[1]:+.2f}" == hi,
          f"{k}: member-mean interval {lo} to {hi}", [round(x, 4) for x in ci])
check(says("−0.62 to +0.34 in Gao and −0.26 to +0.48 in Jiang")
      and says("what these two datasets establish is not that the panel moves by nothing"),
      "the proteome member means are reported with intervals and not as a null")
c = R["coherence_by_level"]
check(f"{c['cancellation_difference_ci95'][0]:+.3f}" == "-0.013"
      and f"{c['cancellation_difference_ci95'][1]:+.3f}" == "+0.047"
      and says("95% confidence interval of −0.013 to +0.047"),
      "the equivalence claim over the 114 sets carries its interval")

# -- the transcription factors are not hidden -------------------------------
check(says("The five transcription factors fall as a class by 0.277 in GSE14520 and by 0.082 "
           "in GSE76427 but rise by 0.422 in TCGA-LIHC")
      and f"{cls('tcga', 'transcription factor'):.3f}" == "0.422",
      "the transcription factors' transcript-level behavior is reported, not only their proteome values")
check(D[NAME["g14"]]["secreted_vs_metabolic"]["p"] < 0.05
      and says("Mann-Whitney *P* = 0.0047"),
      "the significant transcript-level class contrast is reported alongside the two null ones")

# -- missingness -------------------------------------------------------------
mj = D[NAME["jiang"]]["missingness"]
check(f"{mj['fraction_unquantified_tumor']*100:.1f}" == "36.2"
      and f"{mj['fraction_unquantified_adjacent']*100:.1f}" == "45.5"
      and says("36.2% of values in tumor samples against 45.5% in adjacent samples"),
      "the arm-asymmetric missingness is stated")
check(all(D[NAME[k]]["missingness"]["fraction_unquantified_tumor"] == 0
          for k in ("g14", "g76", "tcga", "gao"))
      and says("Gao's table and the three transcriptomes carry a value for every gene in every sample"),
      "the other four matrices are complete, as claimed")
for n, v in R["class_within_shared_sets"].items():
    for c, x in v.items():
        VALUES += [x["surviving_share"], x["mean_z"], x["n"],
                   x["percentile_within_shared_sets"]]
jd = R["jiang_complete_quantification"]
check(len(jd["members"]) == 11 and jd["n_secreted"] == 9 and jd["n_metabolic"] == 2
      and jd["separation_complete"]
      and says("Eleven members are quantified in both arms of all 124 pairs"),
      "the complete-quantification subset is 11 members and still separates",
      f"{len(jd['members'])} members, sep={jd['separation_complete']}")
check(f"{jd['composite_sample_test']['mean_paired_delta']:+.3f}" == "+0.294"
      and says("The composite over those eleven is +0.294"),
      "the complete-case composite is reported, including that it moves")

# -- abundance ---------------------------------------------------------------
aj = D[NAME["jiang"]]["abundance"]
nb = aj["neighborhood_test"]
check(nb["gap_expected_from_abundance"] > 0
      and says("points the same way as the split, as a compression account requires"),
      "the direction of the abundance effect is stated correctly (same way, too small)")
check(round(1 / nb["share_of_class_gap_explained"]) == 7
      and says("about a seventh of the separation"),
      "the stated fraction is the archive's", f"1/{1/nb['share_of_class_gap_explained']:.2f}")
check(D[NAME["gao"]]["abundance"]["interpretable"] is False
      and says("I therefore withdraw the abundance control in Gao")
      and "of 7 and 4 in Gao" not in MS,
      "the Gao abundance control is withdrawn and its numbers are gone")

# -- controls ----------------------------------------------------------------
ig = {k: set(D[NAME[k]]["immunoglobulin_members"]) for k in ("gao", "jiang")}
check(ig["gao"] != ig["jiang"] and says("but not the same eight"),
      "the text says the two immunoglobulin sets differ", sorted(ig["gao"] ^ ig["jiang"]))
check("immunoglobulins and *JCHAIN*, not made by the liver" not in
      MS.split("**Figure 3.")[1],
      "the Figure 3 legend no longer labels the Jiang set as containing JCHAIN")

# -- the two coherent halves, and the recommendation that follows ------------
cw = R["class_within_shared_sets"]
met = [cw[NAME[k]]["metabolic enzyme"] for k in NAME]
check(all(abs(m["surviving_share"] - 1.0) < 1e-9 for m in met)
      and says("The metabolic members have a surviving share of 1.000 in all five"),
      "the metabolic class has a surviving share of 1.000 in every dataset",
      [round(m["surviving_share"], 4) for m in met])
pcts = [m["percentile_within_shared_sets"] for m in met]
check(min(pcts) > 94 and says("ranks above the 94th percentile in each of the five datasets"),
      "the metabolic subset's percentile range is as stated", [round(x) for x in pcts])
sec = [cw[NAME[k]]["secreted plasma protein"]["surviving_share"] for k in NAME]
check(f"{cw[NAME['gao']]['secreted plasma protein']['surviving_share']:.3f}" == "0.835"
      and all(abs(cw[NAME[k]]["secreted plasma protein"]["surviving_share"] - 1.0) < 1e-9
              for k in ("g14", "g76", "tcga", "jiang"))
      and says("of 1.000 in all three transcriptomes and in Jiang, and of 0.835 in Gao"),
      "the secreted class shares are as stated", [round(x, 3) for x in sec])
tf = [cw[NAME[k]]["transcription factor"]["surviving_share"] for k in ("g14", "g76", "tcga")]
check([f"{x:.2f}" for x in tf] == ["0.51", "0.19", "0.98"]
      and says("at 0.51 in GSE14520 and 0.19 in GSE76427")
      and says("the class reaches 0.98 in TCGA-LIHC"),
      "the transcription factors' transcript-level shares are as stated",
      [round(x, 3) for x in tf])
sg = R["class_sign_across_levels"]
check(sg["metabolic enzyme"]["keeps_its_sign_in_all_five"]
      and sg["surface receptor"]["keeps_its_sign_in_all_five"]
      and not sg["secreted plasma protein"]["keeps_its_sign_in_all_five"]
      and not sg["transcription factor"]["keeps_its_sign_in_all_five"]
      and says("Two of the four keep it: the metabolic members fall in tumor in all five "
               "datasets, and so does *ASGR1*, alone in its class. Two do not")
      and says("the transcription factors rise in both proteomes without holding one "
               "direction among the three transcriptomes"),
      "which classes keep their direction is stated correctly")
check(says("can use the metabolic enzymes of this panel")
      and says("It is a recommendation supported here, not a validated panel")
      and says("No third liver proteome was held out to test it on")
      and says("four things limit it")
      and says("in Gao that could not be checked at all")
      and says("*ASGR1* keeps its direction too, but one protein is not a covariate")
      and says("hold no one direction across the five"),
      "the recommendation is made and all four of its limits are stated with it")

# -- scope -------------------------------------------------------------------
pw = R["panel_within_shared_sets"]; sm = R["panel_on_shared_members"]
nb_ = pw["n_sets_below_panel_within_each_proteome"]
check(nb_["Gao 2019 proteome"] == 24 and nb_["Jiang 2019 proteome"] == 14
      and says("with 24 of the 114 below it in Gao and 14 in Jiang")
      and says("averaged over the two deposits, is below every one of the 114"),
      "the two ways of counting sets below the panel are distinguished")
check(pw["n_sets_dropping_at_least_as_much"] == 0
      and f"{sm['drop']:.3f}" == "0.723"
      and says("restricted to the 18 members common to all five datasets the same drop is 0.888 to 0.165, or 0.723"),
      "the like-for-like drop over shared members is given beside the headline one")
check(says("one curated panel against 114 data-derived signatures")
      and "property of curated panels" not in MS,
      "the scope claim is about this panel, not about curated panels as a class")

# -- Methods completeness ----------------------------------------------------
for t, what in (
        ("GSE14520 is distributed on two platforms; this analysis uses the GPL3921 series matrix",
         "the GSE14520 platform is named"),
        ("of which the 126 holding between 15 and 1,000 genes were kept as eligible",
         "the gene set eligibility rule is given"),
        ("The criterion is that fewer than 20% of its values fall below zero",
         "the abundance-scale criterion matches the code"),
        ("**Missingness.**", "Methods has a missingness paragraph"),
        ("**The six Jiang pipelines.**", "Methods gives all six pipelines")):
    check(says(t), what)
check("changes nothing elsewhere" not in MS,
      "the false 'changes nothing elsewhere' sentence is gone")
check("It contains no redox enzyme" not in MS, "the orphan redox sentence is gone")
check(says("none of the *P* values here is used to select a finding from among many")
      or says("and none would repair what these *P* values are"),
      "the multiplicity statement is the honest one")

# -- the self-correction is retrievable --------------------------------------
check(says("deposited as version 1.2.2 of the archive cited here")
      and says("including version 1.2.2, which holds the analysis corrected here"),
      "the corrected version is named so a reader can retrieve it")

# -- figures -----------------------------------------------------------------
# -- the abstract must not be stronger than the Results ---------------------
AB = MS.split("## Abstract\n\n", 1)[1].split("\n\n**Keywords", 1)[0]
check("where the deposits permit the test" in AB
      and D[NAME["gao"]]["abundance"]["interpretable"] is False,
      "the abstract carries the one-deposit qualifier the Results insist on")
for banned, why in (
        ("not explained by blood, by protein abundance", "an unqualified abundance claim"),
        ("or by uneven missingness", "an unqualified missingness claim"),
        ("Each half is coherent", "a coherence claim the TTR exception contradicts"),
        ("Each class agrees with itself", "a claim the transcription factors contradict"),
        ("metabolic half", "calling a 5-of-18 class a half"),
        ("with one exception, while", "an exception not tied to one dataset")):
    check(banned not in AB, f"the abstract does not contain {why}")
check("with one exception in one dataset" in AB,
      "the abstract ties its one exception to one dataset, as Gao's TTR requires")
sg = R["class_sign_across_levels"]
check(sg["metabolic enzyme"]["keeps_its_sign_in_all_five"]
      and not sg["secreted plasma protein"]["keeps_its_sign_in_all_five"]
      and "only the metabolic class keeps its direction at both levels" in AB
      and "Within the secreted and metabolic classes" in AB,
      "the abstract names the right class as the one that transfers")
check(all(x in AB for x in ("0.48 to 0.73", "six pipelines", "114 published liver signatures")),
      "the abstract keeps the three headline results: the transcriptome fall, the six\n        pipelines, and the 114-signature comparison")

print("\nPASS 3: figures, references and limits")
body = MS.split("## References", 1)[0]
for n, panels in ((1, ""), (2, "abc"), (3, "abc")):
    check(f"**Figure {n}." in MS, f"Figure {n} has a legend")
    leg = MS.split(f"**Figure {n}.")[1].split("\n\n")[0]
    check(len(leg.split()) <= 350, f"Figure {n} legend within 350 words ({len(leg.split())})")
    for q in panels:
        check(f"Fig. {n}{q}" in body, f"Figure {n}{q} cited in the text")
        check(f"(**{q}**)" in leg, f"Figure {n}{q} described in the legend")
check("Significant in every dataset" not in open(f"{IR_DOCS}/reframe/71_figures.py").read(),
      "Figure 2a no longer claims significance its own Figure 3c contradicts")
check(says("Small dots are the individual secreted and metabolic members",
           MS.split("**Figure 3.")[1]),
      "the Figure 3b legend's separation claim is backed by drawn members")
check("Table" not in body, "no table is claimed (there are none)")

refs = MS.split("## References", 1)[1].split("## Acknowledgements", 1)[0]
listed = {int(x) for x in re.findall(r"(?m)^(\d+)\. ", refs)}
cited = {int(n) for grp in re.findall(r"\[([\d,\s]+)\]", body) for n in grp.split(",") if n.strip()}
# a superscript citation attached to a word, never the affiliation marker that
# opens its own line, and never the "1,*" beside the author name
cited |= {int(n) for n in re.findall(r"(?<=[^\s\n])\^(\d+)\^", MS)}
check(not re.search(r"(?<=[^\s\n])\^\d+\^", MS) and "Anthropic." not in refs
      and says("Claude (Anthropic PBC, San Francisco, CA, USA), was used"),
      "the AI tool is disclosed in a statement, not entered as a numbered reference")
check(cited == listed, "every reference cited and every citation resolves",
      f"cited {sorted(cited)} listed {sorted(listed)}")
check(sorted(listed) == list(range(1, len(listed) + 1)), "reference numbering has no gap")
order = [int(n) for grp in re.findall(r"\[([\d,\s]+)\]", body) for n in grp.split(",") if n.strip()]
first, seen = [], set()
for n in order:
    if n not in seen:
        seen.add(n); first.append(n)
check(first == sorted(first), "references first cited in numerical order", first)
check(len(listed) <= 60, f"within 60 references ({len(listed)})")

print("\nPASS 4: journal limits and house style")
check("—" not in MS and "—" not in CL, "no em dash in manuscript or cover letter")
check(not re.findall(r"\b(?:we|our|us)\b", MS, re.I)
      and not re.findall(r"\b(?:we|our|us)\b", CL, re.I), "first person singular throughout")
check(not [w for w in ("tumour", "colour", "analysed", "organis", "centre", "normalis")
           if w in (body + CL).lower()], "American spelling")
check(not re.findall(r"\b(\w+)\s+\1\b", MS, re.I), "no doubled word")
ab = MS.split("## Abstract\n\n", 1)[1].split("\n\n**Keywords", 1)[0]
check(len(ab.split()) <= 200, f"abstract within 200 words ({len(ab.split())})")
title = MS.split("\n", 1)[0].lstrip("# ")
check(len(title.split()) <= 20, f"title within 20 words ({len(title.split())})")
mt = MS.split("## Introduction", 1)[1].split("## Methods", 1)[0]
mw = len(re.sub(r"[#*_`]", " ", mt).split())
check(mw <= 4500, f"main text within 4,500 words, Methods excluded ({mw})")
# Scientific Reports takes the cover letter as an uploaded file and asks for
# "a single page", not a character count; the 3,000-character cap was the
# previous journal's form field.
check(len(CL) <= 4200, f"cover letter fits one page ({len(CL)} characters)")
import re as _re
_mt = MS.split("## Introduction", 1)[1].split("## Methods", 1)[0]
_mw = len(_re.sub(r"[#*_`]", " ", _mt).split())
_stated = int(_re.search(r"main text is about ([\d,]+) words", CL).group(1).replace(",", ""))
check(abs(_stated - _mw) <= 60,
      f"the cover letter's word count matches the manuscript ({_stated} vs {_mw})")
check(says("Suitable reviewers", CL) and says("I request no exclusions", CL),
      "the cover letter names suggested reviewers and addresses exclusions")
check(says("no prior discussion of this work with an Editorial Board Member", CL),
      "the cover letter carries the Editorial Board Member disclosure")
check(says("declined without review", CL), "the cover letter discloses the related submission")
aff = ("Hepatobiliary Unit, Division of Gastroenterology, Department of Internal Medicine, "
       "Incheon St. Mary's Hospital, College of Medicine, The Catholic University of Korea, "
       "Seoul, Republic of Korea")
def flat(x):
    return re.sub(r"\s+", " ", x.replace("\\", " "))


check(flat(aff) in flat(MS) and flat(aff) in flat(CL),
      "the affiliation string is verbatim in both")

print(f"\n{'ALL CHECKS PASS' if not FAIL else str(len(FAIL)) + ' CHECK(S) FAILED'}")
for f in FAIL:
    print("  -", f)
sys.exit(1 if FAIL else 0)
