#!/usr/bin/env python3
"""Anchored claim checking.

`50_consistency_check.py` asks whether a number exists anywhere in the archive.
With forty thousand archived values almost any number does, so that test has no
teeth. This module instead binds each claim in the manuscript to the exact
archived quantity it is about, and fails when the two differ.

Each anchor is (name, regex with one capture group, a function of the archive).
"""
from __future__ import annotations

import json
import os
import re
import statistics as st
import numpy as np

RESULTS = "/home/claude/rsi/results"



SUPMAP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5",
          "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-"}


def sci(text: str) -> float:
    """Parse '1.0 × 10⁻¹⁹' or '0.0024' into a float."""
    t = text.strip().replace("−", "-")
    m = re.match(r"([\d.]+)\s*×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", t)
    if m:
        exp = "".join(SUPMAP.get(c, c) for c in m.group(2))
        return float(m.group(1)) * 10 ** int(exp)
    return float(t)


def J(name: str):
    return json.load(open(os.path.join(RESULTS, name), encoding="utf-8"))


def sigs(f: str):
    return [s for s in J(f)["signatures"] if s.get("status") == "ok"]


def med(vals):
    return st.median(vals)


def build_anchors():
    """Return a list of dicts: name, pattern, expected, tol."""
    B14 = sigs("SIGNATURE_BENCHMARK_GSE14520.json")
    BLI = sigs("SIGNATURE_BENCHMARK_TCGA_LIHC.json")
    BLU = sigs("SIGNATURE_BENCHMARK_TCGA_LUAD.json")
    BKI = sigs("SIGNATURE_BENCHMARK_TCGA_KIRC.json")
    B76 = sigs("SIGNATURE_BENCHMARK_GSE76427.json")
    r14 = [s["retention_joint"] for s in B14]
    C14 = J("GSE14520_composition.json")
    C76 = J("GSE76427_composition.json")
    U14 = J("COVARIATE_COMPARISON_6u_GSE14520.json")
    ULU = J("COVARIATE_COMPARISON_6u_TCGA_LUAD.json")
    UKI = J("COVARIATE_COMPARISON_6u_TCGA_KIRC.json")
    V14 = J("RETENTION_INTERVALS_6v_GSE14520.json")
    V76 = J("RETENTION_INTERVALS_6v_GSE76427.json")
    VLI = J("RETENTION_INTERVALS_6v_TCGA_LIHC.json")
    VLU = J("RETENTION_INTERVALS_6v_TCGA_LUAD.json")
    VKI = J("RETENTION_INTERVALS_6v_TCGA_KIRC.json")
    SIM = J("ESTIMATOR_SIMULATION_6v.json")
    PW = J("TCGA_LIHC/PREMISE_6w.json")
    PX = J("TCGA_KIRC/PREMISE_6x.json")
    PY = J("TCGA_LUAD/PREMISE_6y.json")
    KW = J("TCGA_LIHC/CONSOLIDATED_6w.json")
    GAO = J("PROTEIN_VALIDATION_6s_Gao2019.json")
    JIA = J("PROTEIN_VALIDATION_6s_Jiang2019.json")
    MR = J("GAO_MRNA_PROTEIN_6aa.json")
    HE = J("GAO_HE_PURITY_6aa.json")
    MS_ = J("MS_ARTIFACT_6ab.json")
    PZ = J("PROGNOSTIC_6z_TCGA_LIHC.json")
    DIR = J("DIRECTION_6ad.json")
    FLOW = J("ENUMERATION_FLOW.json")
    GAP = J("ARCHIVE_GAPS_2026-09-09.json")
    AUD = J("AUDIT_RECOMPUTE_2026-09-09.json")

    A = []

    def add(name, pattern, expected, tol=0.0006, scale=1.0):
        A.append(dict(name=name, pattern=pattern, expected=float(expected),
                      tol=tol, scale=scale))

    # ---- the benchmark distribution in liver
    add("liver median retention", r"identity-retention fraction in the joint model had a median of \*\*([\d.]+)\*\*", med(r14))
    add("liver IQR low", r"interquartile range ([\d.]+) to [\d.]+; Supplementary Tables", np.percentile(r14, 25))
    add("liver IQR high", r"interquartile range [\d.]+ to ([\d.]+); Supplementary Tables", np.percentile(r14, 75))
    add("liver n below 50%", r"\*\*(\d+) of 119 \(25%\)\*\* fell below 50%", sum(v < .5 for v in r14), tol=0.5)
    add("liver n evaluable", r"Across those (\d+) signatures the identity-retention", len(B14), tol=0.5)
    add("liver range low", r"ranged from ([\d.]+) to [\d.]+, median 0\.758", min(r14))
    add("liver range high", r"ranged from [\d.]+ to ([\d.]+), median 0\.758", max(r14))

    # ---- this study's own modules
    add("REDUCTION joint GSE14520", r"REDUCTION retained \*\*(\d+)%\*\* \(intercept \+0\.764",
        C14["module_reduction_joint"]["retained_fraction"] * 100, tol=0.6)
    add("REDUCTION joint intercept", r"intercept \+([\d.]+), 95% CI 0\.614", C14["module_reduction_joint"]["estimate"])
    add("DRAIN joint GSE14520", r"DRAIN's \*\*(\d+)%\*\* \(−0\.248",
        C14["module_drain_joint"]["retained_fraction"] * 100, tol=0.6)
    add("REDUCTION D1 GSE14520", r"After adjustment for D1 alone, REDUCTION retained \*\*(\d+)%\*\*",
        U14["own_modules"]["REDUCTION"]["retention_D1"] * 100, tol=0.6)
    add("DRAIN D1 GSE14520", r"REDUCTION retained \*\*84%\*\* and DRAIN \*\*(\d+)%\*\* in GSE14520",
        U14["own_modules"]["DRAIN"]["retention_D1"] * 100, tol=0.7)
    add("REDUCTION reproduction", r"\(REDUCTION ([\d.]+)%, DRAIN 31\.03%\)",
        U14["own_modules"]["REDUCTION"]["retention_joint"] * 100, tol=0.01)
    add("DRAIN reproduction", r"\(REDUCTION 89\.05%, DRAIN ([\d.]+)%\)",
        U14["own_modules"]["DRAIN"]["retention_joint"] * 100, tol=0.01)
    add("GSE76427 DRAIN joint", r"−0\.364 \(−0\.659 to −0\.069\), \*P\* = 0\.017, \*\*([\d.]+)%\*\*",
        C76["module_drain_joint"]["retained_fraction"] * 100, tol=0.06)

    # ---- 6u
    add("composition-only median liver", r"the median signature retained \*\*([\d.]+)\*\* of its unadjusted shift",
        U14["retention_composition_only"]["median"])
    add("composition-only below 50 liver", r"only \*\*(\d+) of 119\*\* fell below 50%",
        U14["retention_composition_only"]["n_below_50pct"], tol=0.5)
    add("identity-only below 50 liver", r"against (\d+) of 119 under the identity-only model",
        U14["retention_identity_only"]["n_below_50pct"], tol=0.5)
    add("6u Spearman covariates", r"Spearman ρ = −([\d.]+), \*P\* = 1\.2 × 10⁻⁴",
        abs(U14["covariate_association"]["spearman_rho"]))
    add("lung identity-only median", r"leaving ([\d.]+) against 0\.845", ULU["retention_composition_only"]["median"])
    add("lung composition-only median", r"leaving 0\.699 against ([\d.]+)", ULU["retention_identity_only"]["median"])
    add("lung reclassified", r"(\d+) of 36 signatures change side", ULU["agreement"]["n_reclassified_at_50pct"], tol=0.5)
    add("kidney identity-only median", r"([\d.]+) under identity against 0\.651", UKI["retention_identity_only"]["median"])
    add("kidney composition-only median", r"0\.666 under identity against ([\d.]+)", UKI["retention_composition_only"]["median"])
    add("kidney reclassified", r"(\d+) of the 15 signatures change side", UKI["agreement"]["n_reclassified_at_50pct"], tol=0.5)

    # ---- 6v intervals
    for lab, D, pat in [("GSE14520", V14, r"Median interval widths are ([\d.]+) in the 213-pair"),
                        ("TCGA-LIHC", VLI, r"are 0\.327 in the 213-pair liver cohort, ([\d.]+) and"),
                        ("GSE76427", V76, r"0\.392 and ([\d.]+) in the two small liver cohorts"),
                        ("TCGA-KIRC", VKI, r"at 50 and 52 pairs, ([\d.]+) in kidney"),
                        ("TCGA-LUAD", VLU, r"([\d.]+) in lung, and the number of signatures")]:
        add(f"median interval width {lab}", pat, D["median_ci_width"])
    add("supported counts", r"entirely below 50% is (\d+), 4, 0, 3 and 0", V14["n_below_50_interval"], tol=0.5)

    # ---- simulation
    cov = [c["coverage_95"] for c in SIM["cells"]]
    add("simulation cells", r"pre-registered grid of (\d+) conditions", len(SIM["cells"]), tol=0.5)
    add("coverage low", r"ranged from ([\d.]+) to 0\.980", min(cov))
    add("coverage high", r"ranged from 0\.918 to ([\d.]+)", max(cov))
    bias = sorted((abs(c["bias"]) for c in SIM["cells"]), reverse=True)
    add("largest bias", r"by 0\.054 and ([\d.]+), the only two cells", bias[0])
    add("second bias", r"by ([\d.]+) and 0\.061, the only two cells", bias[1])

    # ---- premise checks
    add("kidney identity fall", r"Identity falls by ([\d.]+) z units in kidney", abs(PX["identity_paired_delta"]))
    add("lung identity fall", r"and by ([\d.]+) in lung", abs(PY["identity_paired_delta"]))
    add("C1 lung", r"C1 falls in lung \(−([\d.]+),", abs(PY["C1_paired_delta"]))
    add("C1 liver", r"does not move in liver \(\+([\d.]+),", PW["C1_paired_delta"])
    add("C1 kidney", r"rises sharply in kidney \(\+([\d.]+),", PX["C1_paired_delta"])

    # ---- third cohort
    add("TCGA-LIHC median", r"Median retention is \*\*([\d.]+)\*\*", med([s["retention_joint"] for s in BLI]))
    add("TCGA-LIHC evaluable", r"(\d+) of 103 sets fall below half", sum(s["retention_joint"] < .5 for s in BLI), tol=0.5)
    add("median difference 6w", r"a difference of ([\d.]+) and inside the 0\.15 margin",
        abs(KW["benchmark"]["median_difference_vs_GSE14520"]))
    add("shared spearman", r"correlate at Spearman ρ = \*\*\+([\d.]+)\*\*", KW["benchmark"]["spearman_shared"])
    add("REDUCTION TCGA retention", r"retains \*\*([\d.]+)%\*\* \(interval 0\.851",
        KW["own_modules"]["REDUCTION_4gene"]["retention_joint"] * 100, tol=0.06)
    add("DRAIN prime TCGA retention", r"with a retention of \*\*([\d.]+)%\*\*",
        KW["own_modules"]["DRAIN_without_POR_6l"]["retention_joint"] * 100, tol=0.06)

    # ---- tissues
    add("lung median", r"Median retention is ([\d.]+) in lung", med([s["retention_joint"] for s in BLU]))
    add("kidney median", r"and ([\d.]+) in kidney across 15", med([s["retention_joint"] for s in BKI]))

    # ---- protein
    add("Gao D1 premise", r"\*\*holds\*\* \(−([\d.]+), \*P\* = 3\.3", abs(GAO["positive_control_D1"]["mean_paired_delta"]))
    add("Jiang D1 premise", r"\*\*fails\*\* \(\+([\d.]+) higher", JIA["positive_control_D1"]["mean_paired_delta"])
    add("Gao DRAIN retention", r"retained \*\*([\d.]+)%\*\* after adjustment for D1",
        GAO["modules"]["DRAIN"]["D1"]["retention"] * 100, tol=0.06)
    add("Gao REDUCTION shift", r"REDUCTION showed \*\*no shift at all\*\* \(−([\d.]+),",
        abs(GAO["modules"]["REDUCTION"]["unadjusted_mean_delta"]))
    add("mRNA-protein rho", r"\*\*ρ = \+([\d.]+) \(\*P\* = 0\.36\)\*\*", MR["spearman_rho"])
    add("purity vs C1", r"Spearman ρ = \*\*−([\d.]+)\*\*, \*P\* = 1\.2 × 10⁻¹³", abs(HE["purity_vs_C1"]["rho"]))
    add("purity vs D1", r"\(ρ = \+([\d.]+), \*P\* = 0\.15\)", HE["purity_vs_D1"]["rho"])
    add("MS abundance rho", r"gives Spearman ρ = \*\*−([\d.]+)\*\*", abs(MS_["predictors"]["mean_abundance"]["spearman_rho"]))
    add("MS coverage rho", r"proteome coverage ρ = \*\*\+([\d.]+)\*\*", MS_["predictors"]["coverage"]["spearman_rho"], tol=0.0006)
    add("CYB5R3 percentile", r"\*\*(\d+)th percentile of abundance\*\*", MS_["own_modules"]["CYB5R3"]["percentile"], tol=0.5)

    # ---- survival
    add("6z rho", r"\*\*Spearman ρ = \+([\d.]+) \(\*P\* = 0\.0024\)\*\*", PZ["spearman_rho"])
    add("6z HR identity", r"hazard ratio of ([\d.]+) per standard deviation \(0\.781", PZ["D1"]["HR_per_SD"])

    # ---- direction, 6ad
    add("6ad rho liver", r"\(ρ = \+([\d.]+), \*P\* = 1\.2 × 10⁻¹⁰\)", DIR["GSE14520"]["rho_retention_vs_signed_shift"])
    add("6ad rise median", r"retain a median ([\d.]+) against 0\.532", DIR["GSE14520"]["median_rise"])
    add("6ad fall median", r"median 0\.887 against ([\d.]+) for the 59", DIR["GSE14520"]["median_fall"])
    add("6ad size rho", r"Set size is unrelated to retention \(Spearman ρ = \+([\d.]+)", DIR["GSE14520"]["rho_retention_vs_set_size"], tol=0.006)
    add("6ad lung rho", r"absent in lung \(ρ = \+([\d.]+)\)", abs(DIR["TCGA_LUAD"]["rho_retention_vs_signed_shift"]), tol=0.0015)

    # ---- enumeration
    add("liver keyword hits", r"That enumeration returned (\d+) sets, of which 126", FLOW["liver"]["keyword"], tol=0.5)
    add("liver size-eligible", r"of which (\d+) carry between 15 and 1000 symbols", FLOW["liver"]["size_eligible"], tol=0.5)
    add("kidney enumerated", r"the enumeration returned (\d+) sets, of which 15", FLOW["kidney"]["keyword"], tol=0.5)

    # ---- Table 4 and the archive gaps
    g = GAP["GSE14520"]["per_gene"]
    add("Table 4 MTARC2 GSE14520", r"\| \*MTARC2\* \| −([\d.]+) \(1\.2 × 10⁻²⁴\)", abs(g["MTARC2"]["paired_mean_delta"]))
    add("Table 4 POR GSE14520", r"\| \*POR\* \| −([\d.]+) \(8\.8 × 10⁻²³\)", abs(g["POR"]["paired_mean_delta"]))
    add("Table 4 CYB5R1 GSE14520", r"\| \*CYB5R1\* \| \+([\d.]+) \(3\.3 × 10⁻²⁷\)", g["CYB5R1"]["paired_mean_delta"])
    add("Table 4 NQO1 GSE14520", r"\| \*NQO1\* \| \+([\d.]+) \(5\.4 × 10⁻²⁰\)", g["NQO1"]["paired_mean_delta"])
    add("SUPPLY retention", r"it retained ([\d.]+)% of its shift after joint adjustment",
        GAP["GSE14520"]["SUPPLY_joint"]["retention"] * 100, tol=0.06)
    add("restricted D1 REDUCTION", r"REDUCTION from 83\.6% to ([\d.]+)%",
        AUD["gse14520_restricted_D1_sensitivity"]["D1_19_protein_measured"]["REDUCTION"]["retention"] * 100, tol=0.06)
    add("restricted D1 DRAIN", r"moves DRAIN retention only from 24\.7% to ([\d.]+)%",
        AUD["gse14520_restricted_D1_sensitivity"]["D1_19_protein_measured"]["DRAIN"]["retention"] * 100, tol=0.06)

    # ---- sample counts
    add("TCGA-LIHC samples", r"the UCSC Xena GDC hub in log2\(transcripts per million \+ 1\), (\d+) samples",
        GAP["tcga_total_samples"]["TCGA_LIHC"], tol=0.5)
    add("TCGA-LUAD samples", r"TCGA-LUAD, (\d+) samples giving 58 pairs", GAP["tcga_total_samples"]["TCGA_LUAD"], tol=0.5)
    add("TCGA-KIRC samples", r"TCGA-KIRC, (\d+) samples giving 72 pairs", GAP["tcga_total_samples"]["TCGA_KIRC"], tol=0.5)

    # ---- §6ae: the matched-random negative control and the curvature check --
    NC = J("NEGATIVE_CONTROL_6ae.json")
    nc, cv = NC["negative_control"], NC["curvature"]
    add("6ae random matched median",
        r"leaves a median of \*\*([\d.]+)\*\* against D1's",
        nc["random_median_of_medians"])
    add("6ae D1 median in the control",
        r"against D1's \*\*([\d.]+)\*\*, and puts", nc["D1_median"])
    add("6ae random below 50",
        r"puts \*\*(\d+)\*\* signatures below half against D1's",
        nc["random_n_below_50_median"], tol=0.5)
    add("6ae D1 below 50",
        r"below half against D1's \*\*(\d+)\*\*", nc["D1_n_below_50"], tol=0.5)
    add("6ae correlation",
        r"correlate at median \*r\* = \*\*(\d+\.\d+)\*\*",
        nc["median_r_with_D1_retention"])
    add("6ae fraction at or below",
        r"and (\d+\.\d+)% of matched random covariates remove at least as much",
        nc["fraction_random_at_or_below_D1_median"] * 100, tol=0.06)
    add("6ae no-shift median",
        r"does not move between tumor and adjacent tissue leaves a median of (\d+\.\d+)",
        nc["other_pools"]["no_tumor_shift"]["median_of_medians"])
    add("6ae rising median",
        r"one drawn from genes rising in tumor leaves (\d+\.\d+)",
        nc["other_pools"]["rises_in_tumor"]["median_of_medians"])
    ps = J("POOL_SHIFT_SUMMARY_6ae.json")["pools"]
    add("6ae unrestricted own shift",
        r"whose own paired shift has a median of \+(\d+\.\d+) because",
        ps["unrestricted"]["median_own_paired_shift"])
    add("6ae unrestricted median",
        r"whose own paired shift has a median of \+\d+\.\d+ because most genes "
        r"hardly move between the two tissues, leaves (\d+\.\d+)",
        nc["other_pools"]["unrestricted"]["median_of_medians"])
    add("6ae reproduction",
        r"reproduced all 119 archived values to within (\d+\.\d+) first",
        NC["reproduction_of_archive"]["max_abs_diff"], tol=0.00006)
    add("6ae quadratic median",
        r"moves the median from \d+\.\d+ to \*\*(\d+\.\d+)\*\*", cv["quadratic_median"])
    add("6ae quadratic below 50",
        r"count below half from 30 to \*\*(\d+)\*\* of 119",
        cv["quadratic_n_below_50"], tol=0.5)
    add("6ae quadratic spearman",
        r"rankings correlated at Spearman (\d+\.\d+)",
        cv["spearman_linear_quadratic"])

    # ---- P values, compared on the log scale because they are written to one digit
    def addp(name, pattern, expected):
        A.append(dict(name=name, pattern=pattern, expected=float(expected),
                      tol=0.0, scale=1.0, kind="p"))

    addp("P REDUCTION joint", r"95% CI 0\.614 to 0\.913, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)",
         C14["module_reduction_joint"]["p"])
    addp("P DRAIN joint", r"−0\.390 to −0\.106, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)",
         C14["module_drain_joint"]["p"])
    addp("P Gao premise", r"\*\*holds\*\* \(−0\.172, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)",
         GAO["positive_control_D1"]["p"])
    addp("P purity vs C1", r"ρ = \*\*−0\.544\*\*, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)",
         abs(HE["purity_vs_C1"]["p"]))
    addp("P 6z", r"\*\*Spearman ρ = \+0\.296 \(\*P\* = ([\d.]+)\)\*\*", PZ["spearman_p"])
    addp("P 6ad liver", r"ρ = \+0\.548, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)",
         DIR["GSE14520"]["p_retention_vs_signed_shift"])
    addp("P kidney identity", r"in kidney \(\*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)", PX["identity_p"])
    addp("P lung identity", r"1\.230 in lung \(\*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)", PY["identity_p"])
    addp("P C1 lung", r"C1 falls in lung \(−0\.492, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)", PY["C1_p"])
    addp("P C1 kidney", r"kidney \(\+0\.956, \*P\* = ([\d.]+ × 10[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)\)", PX["C1_p"])
    addp("P mRNA-protein", r"ρ = \+0\.094 \(\*P\* = ([\d.]+)\)", MR["spearman_p"])
    return A



def run(ms_text: str):
    """Return (checked, failures, unmatched_patterns)."""
    anchors = build_anchors()
    t = ms_text.replace("−", "−")
    fails, missing, n = [], [], 0
    for a in anchors:
        m = re.search(a["pattern"], t)
        if not m:
            missing.append(a["name"])
            continue
        n += 1
        if a.get("kind") == "p":
            got = sci(m.group(1))
            want = a["expected"]
            # a P value is written to one or two digits; compare on the log scale
            import math
            if want <= 0 or got <= 0:
                ok = got == want
            else:
                ok = abs(math.log10(got) - math.log10(want)) < 0.05
            if not ok:
                fails.append(f"{a['name']}: manuscript {got:.3g}, archive {want:.3g}")
            continue
        got = float(m.group(1).replace("−", "-"))
        want = a["expected"] * a["scale"]
        if abs(got - want) > a["tol"]:
            fails.append(f"{a['name']}: manuscript {got}, archive {round(want, 6)}")
    return n, fails, missing


# ---------------------------------------------------------------------------
# Co-occurrence checking.
#
# Writing an anchor for every number is not practical, but numbers that appear
# together in a sentence came from one place in the pipeline and must appear
# together in the archive. An estimate quoted with an interval that belongs to a
# different model is exactly the error this catches, and a coincidental match on
# three or four values at once does not happen.
# ---------------------------------------------------------------------------
import math


def _q(x):
    """Keep a scalar comparable without destroying small P values.

    Rounding to eight decimals turns 1.16e-10 into 0.0, which silently removed
    every strongly significant P value from the archive and made the claims
    quoting them look unsupported."""
    x = float(x)
    return x if abs(x) < 1e-4 else round(x, 8)


def _own(obj):
    """Scalars written directly on this dict (a list of scalars counts as direct)."""
    out = set()
    if not isinstance(obj, dict):
        return out
    for v in obj.values():
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            out.add(_q(v))
        elif isinstance(v, list) and v and all(
                isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
            out.update(_q(x) for x in v)
    return out


def _nodes(obj, acc):
    """One candidate set per dict node, holding only the scalars written on it.

    Depth is deliberately not merged. Pulling a child's scalars up would let two
    numbers from sibling models support each other, and pulling a list of
    per-signature records up would pool hundreds of unrelated values into one
    set, which is exactly the error class this check exists to find."""
    if isinstance(obj, dict):
        own = _own(obj)
        if own:
            acc.append(own)
        for v in obj.values():
            if isinstance(v, (dict, list)):
                _nodes(v, acc)
        return own
    if isinstance(obj, list):
        for v in obj:
            _nodes(v, acc)
    return set()


_MD_NUM = re.compile(r"[+−-]?\d+\.?\d*(?:[eE][+-]?\d+)?")


def _md_nodes(path, acc):
    """Rows of the archived Markdown reports are nodes too.

    Some quantities, notably the Phase B module trend statistics, were written
    only into a results report and never into a JSON file. Reading each table
    row and each line as a node lets those claims be checked as well."""
    for line in open(path, encoding="utf-8", errors="replace"):
        if not any(ch.isdigit() for ch in line):
            continue
        vals = set()
        for tok in _MD_NUM.findall(line):
            try:
                vals.add(_q(float(tok.replace("−", "-"))))
            except ValueError:
                pass
        if len(vals) > 1:
            acc.append(vals)


def archive_nodes():
    acc = []
    for dp, dn, fn in os.walk(RESULTS):
        for f in sorted(fn):
            full = os.path.join(dp, f)
            if f.endswith(".json"):
                try:
                    d = json.load(open(full, encoding="utf-8"))
                except Exception:
                    continue
                _nodes(d, acc)
            elif f.endswith(".md"):
                _md_nodes(full, acc)
    return acc


def _dec(tok: str) -> int:
    return len(tok.split(".")[1]) if "." in tok else 0


def _match(u, v, d, is_sci=False):
    """Does archived value u carry manuscript value v, written to d decimals?

    Values quoted in scientific notation are compared on the log scale exactly
    as the P value anchors are; an absolute tolerance at the sixth decimal
    would reject 6.8e-4 written for an archived 6.8034e-4, and would also let
    an archived exact zero stand in for any small P value at all."""
    au, av = abs(u), abs(v)
    if is_sci or av < 1e-3:
        if u == 0 or v == 0:
            return u == v
        return abs(math.log10(au) - math.log10(av)) < 0.03
    if u == 0 or v == 0:
        return abs(u - v) <= 0.6 * 10 ** (-d)
    tol = 0.6 * 10 ** (-d)
    return (abs(u - v) <= tol or abs(au - av) <= tol
            or round(u, d) == round(v, d) or round(au, d) == round(av, d))


def _supported(vals, nodes) -> bool:
    """Is there one archived node carrying every value, at its written precision?"""
    for node in nodes:
        if all(any(_match(u, v, min(_dec(tok), 6), sci_) for u in node)
               for tok, v, sci_ in vals):
            return True
    return False


SUP = "[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]"
NUM = r"[+−-]?\d+\.\d+"
PV = r"\d+\.?\d*(?: × 10" + SUP + r"+)?"

TUPLE_PATTERNS = [
    # an estimate, its interval and its P quoted together
    (rf"\(({NUM}), 95% CI ({NUM}) to ({NUM}), \*P\* = ({PV})\)", 4),
    (rf"\(({NUM}), ({NUM}) to ({NUM}), \*P\* = ({PV})\)", 4),
    (rf"\(intercept ({NUM}), 95% CI ({NUM}) to ({NUM}), \*P\* = ({PV})\)", 4),
    # an estimate with its interval
    (rf"\(interval ({NUM}) to ({NUM})", 2),
    (rf"interquartile range ({NUM}) to ({NUM})", 2),
    (rf"an interval of ({NUM}) to ({NUM})", 2),
    (rf"interval for its denominator of ({NUM}) to ({NUM})", 2),
    (rf"is ({NUM}) with an interval of ({NUM}) to ({NUM})", 3),
    (rf"({NUM}), interval ({NUM}) to ({NUM})", 3),
    (rf"95% CI ({NUM}) to ({NUM})", 2),
    (rf"bootstrap interval .{{0,40}}?of ({NUM}) to ({NUM})", 2),
    # a correlation, or a shift, with its P
    (rf"ρ = \*?\*?({NUM})\*?\*?,? \(?\*P\* = ({PV})\)?", 2),
    (rf"\(({NUM}), \*P\* = ({PV})\)", 2),
    (rf"\(\*?\*?({NUM})\*?\*? \(({PV})\)\)", 2),
    (rf"\*z\* = ({NUM}), \*P\* = ({PV})", 2),
    (rf"({NUM}) with \*P\* = ({PV})", 2),
    (rf"\*r\* = ({NUM}), \*R\*² = ({NUM})", 2),
    # table cells: estimate then P in parentheses
    (rf"\| \*?\*?({NUM})\*?\*? \(({PV})\)", 2),
    (rf"\*\*({NUM})\*\* \(({PV})\)", 2),
    (rf"({NUM}) \(({PV})\) \|", 2),
    # table cells: estimate, interval, P and a retention percentage
    (rf"({NUM}) \(({NUM}) to ({NUM})\), \*P\* = ({PV})", 4),
    (rf"({NUM}), \*P\* = ({PV})", 2),
]


_NODE_CACHE = []


def run_tuples(ms_text: str, whole: bool = False):
    """Check every group of numbers quoted together against the archive.

    whole=True scans the text as given, for documents such as the
    pre-registration protocol that have no Abstract-to-References body."""
    if not _NODE_CACHE:
        _NODE_CACHE.extend(archive_nodes())
    nodes = _NODE_CACHE
    body = (ms_text if whole else
            ms_text[ms_text.index("## Abstract"):ms_text.index("\n## References\n")])
    checked, bad = 0, []
    for pat, k in TUPLE_PATTERNS:
        for m in re.finditer(pat, body):
            toks = [m.group(i + 1) for i in range(k)]
            vals = []
            good = True
            for t in toks:
                t2 = t.replace("−", "-").strip()
                try:
                    v = sci(t2) if "×" in t2 else float(t2)
                except ValueError:
                    good = False
                    break
                vals.append((t2.lstrip("+-"), v, "×" in t2))
            if not good:
                continue
            checked += 1
            if not _supported(vals, nodes):
                ctx = body[max(0, m.start() - 70):m.end() + 10].replace("\n", " ")
                bad.append(f"{m.group(0)}  ...{ctx.strip()}...")
    return checked, bad


if __name__ == "__main__":
    txt = open("/home/claude/manuscript_v4.md", encoding="utf-8").read()
    n, fails, missing = run(txt)
    print(f"{n} anchored claims checked against the archive")
    for f in fails:
        print("  FAIL  " + f)
    for m in missing:
        print("  LOST  anchor no longer matches the manuscript: " + m)
    print("\n" + ("ALL ANCHORED CLAIMS AGREE" if not fails and not missing
                  else f"{len(fails)} disagreement(s), {len(missing)} lost anchor(s)"))
    tn, tbad = run_tuples(txt)
    print(f"\n{tn} co-occurring groups checked; {len(tbad)} unsupported")
    for t in tbad:
        print("  UNSUPPORTED  " + t)
