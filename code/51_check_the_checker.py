#!/usr/bin/env python3
"""Mutation test for 50_consistency_check.py.

A checker that only ever passes is worthless. This injects known defects into
copies of the manuscript and confirms each one is caught.
"""
import os, re, shutil, subprocess, sys, tempfile

HOME = "/home/claude"
SRC = f"{HOME}/manuscript_v4.md"
CHK = f"{HOME}/rsi/code/50_consistency_check.py"

MUTATIONS = [
    # companion-file defects: the checker reads these files too
    ("alt text describing the wrong figure",
     "**Figure 8.** Four panels on precision and estimator behavior",
     "**Figure 9.** Four panels on precision and estimator behavior",
     "/home/claude/submission_GigaScience/06_Figure_alt_text.md"),
    ("the cover letter claiming a benchmark the manuscript declines",
     "The premise check and the ratio-stability guard live inside that function",
     "Adjusting for identity is not what adjusting for tumor purity does. "
     "The premise check and the ratio-stability guard live inside that function",
     "/home/claude/CoverLetter_GigaScience.md"),
    ("a literature-sourced percentage silently changed",
     "cirrhosis was recorded in **54%** of the 115 patients",
     "cirrhosis was recorded in **58%** of the 115 patients"),
    ("a British spelling the old word list missed",
     "the dotted gray lines mark", "the dotted grey lines mark"),
    ("Acknowledgements in the British form",
     "**Acknowledgments.** None.", "**Acknowledgements.** None."),
    ("Table 1 calling lung the second tissue",
     "Lung adenocarcinoma; 58 pairs | Third tissue (\u00a76y)",
     "Lung adenocarcinoma; 58 pairs | Second tissue (\u00a76x)"),
    ("a gene set attributed to the protocol that the protocol does not name",
     "I substituted the two directional Gene Ontology sets under the clause of "
     "\u00a76d that requires the set actually used to be recorded with the result, "
     "and I record here that the substitution is mine and was not pre-specified.",
     "These are the two sets \u00a76d named."),
    ("a British spelling left in the cover letter",
     "A tumor specimen and the adjacent tissue",
     "A tumour specimen and the adjacent tissue",
     "/home/claude/CoverLetter_GigaScience.md"),
    ("the cover letter left disagreeing with the manuscript's amendment count",
     "Five amendments in all are not pre-registered",
     "Six amendments in all are not pre-registered"),
    ("the submission metadata file left stale after an abstract revision",
     "median 0.758, 10 keeping an interval below half",
     "median 0.758, 12 keeping an interval below half"),
    # co-occurrence defects: each number is real, only the pairing is wrong
    ("an interval bound taken from a sibling model",
     "(−0.248, −0.390 to −0.106, *P* = 6.8 × 10⁻⁴)",
     "(−0.248, −0.390 to −0.109, *P* = 6.8 × 10⁻⁴)"),
    ("two genes' P values swapped",
     "*CYB5R3* fell (−0.387, *P* = 3.0 × 10⁻⁴) while *AIFM2* rose (+0.308, *P* = 4.9 × 10⁻³)",
     "*CYB5R3* fell (−0.387, *P* = 4.9 × 10⁻³) while *AIFM2* rose (+0.308, *P* = 3.0 × 10⁻⁴)"),
    ("a correlation paired with another correlation's P",
     "ρ = −0.261, *P* = 1.2 × 10⁻⁴", "ρ = −0.261, *P* = 1.2 × 10⁻¹³"),
    ("an estimate replaced by another dataset's estimate",
     "DRAIN fell strongly (−0.944, *P* = 2.3 × 10⁻²¹)",
     "DRAIN fell strongly (−0.716, *P* = 2.3 × 10⁻²¹)"),
    ("a table cell's interval bound moved",
     "+0.373 (0.031 to 0.715), *P* = 0.033",
     "+0.373 (0.031 to 0.735), *P* = 0.033"),
    ("a number changed to one nothing supports",
     "median of **0.758**", "median of **0.799**"),
    ("a retention percentage altered",
     "REDUCTION retained **89%**", "REDUCTION retained **77%**"),
    ("a P value altered",
     "95% CI 0.614 to 0.913, *P* = 1.0 × 10⁻¹⁹", "95% CI 0.614 to 0.913, *P* = 1.0 × 10⁻¹⁴"),
    ("a sample count altered",
     "TCGA-LUAD, 589 samples giving 58 pairs", "TCGA-LUAD, 578 samples giving 58 pairs"),
    ("an interval width altered",
     "Median interval widths are 0.327", "Median interval widths are 0.427"),
    ("a per-gene table value altered",
     "| *MTARC2* | −0.923 (1.2 × 10⁻²⁴)", "| *MTARC2* | −0.823 (1.2 × 10⁻²⁴)"),
    ("a figure legend pushed over 300 words",
     "**Figure 11. What the source proteogenomic study settles (§6aa).**",
     "**Figure 11. What the source proteogenomic study settles (§6aa).** " + "filler word " * 160),
    ("a correlation altered",
     "ρ = +0.548", "ρ = +0.648"),
    ("an em dash introduced",
     "The two arms respond to it differently", "The two arms respond — differently"),
    ("a British spelling introduced",
     "hematoxylin and eosin staining", "haematoxylin and eosin staining"),
    ("first person plural introduced",
     "I defined the identity-retention fraction", "We defined the identity-retention fraction"),
    ("a citation to a reference that does not exist",
     "collection of the Molecular Signatures Database (MSigDB) [", "collection of the Molecular Signatures Database (MSigDB) [99"),
    ("an amendment cited that the protocol lacks",
     "Under §6z, written before any survival value", "Under §6zz, written before any survival value"),
    ("the amendment count made wrong",
     "Thirty-one amendments (§6a to §6ae)", "Thirty amendments (§6a to §6ae)"),
    ("the lock digest altered",
     "`7a2bf934…b680046`", "`7a2bf934…b680047`"),
    ("a table title made too long",
     "**Table 1. Nine public transcriptome series and two published proteomes, with role and pairing.**",
     "**Table 1. Nine public transcriptome series and two published proteomes, each identified here by its accession or its project name, with role and pairing.**"),
    ("the abstract pushed over the word limit",
     "**Conclusions.** How much", "**Conclusions.** " + "padding word " * 40 + "How much"),
]


def run(path: str) -> str:
    env = dict(os.environ)
    src = open(CHK, encoding="utf-8").read().replace(
        f'MS = f"{{HOME}}/manuscript_v4.md"', f'MS = "{path}"')
    tmp = path + ".checker.py"
    open(tmp, "w", encoding="utf-8").write(src)
    out = subprocess.run([sys.executable, tmp], capture_output=True, text=True, env=env)
    os.remove(tmp)
    return out.stdout + out.stderr


def main() -> int:
    base = open(SRC, encoding="utf-8").read()
    tmpdir = tempfile.mkdtemp()
    clean = os.path.join(tmpdir, "manuscript_v4.md")
    shutil.copy(SRC, clean)
    print("baseline (unmutated manuscript)")
    out = run(clean)
    ok0 = "ALL CHECKS PASS" in out
    print(("  PASS  " if ok0 else "  FAIL  ") + "the clean manuscript passes every check")
    caught, missed = 0, []
    print("\nmutations")
    for mut in MUTATIONS:
        # A mutation may name a fourth element: a companion file to damage
        # instead of the manuscript. The checker reads the cover letter, the
        # alt text and the submission metadata by absolute path, so those are
        # mutated in place and restored afterwards.
        name, old, new = mut[0], mut[1], mut[2]
        target = mut[3] if len(mut) > 3 else None
        src = open(target, encoding="utf-8").read() if target else base
        if src.count(old) < 1:
            missed.append(f"{name} (anchor not found)")
            print(f"  SKIP  {name}: anchor not found")
            continue
        p = os.path.join(tmpdir, "manuscript_v4.md")
        open(p, "w", encoding="utf-8").write(base if target else
                                             src.replace(old, new, 1))
        if target:
            open(target, "w", encoding="utf-8").write(src.replace(old, new, 1))
        try:
            out = run(p)
        finally:
            if target:
                open(target, "w", encoding="utf-8").write(src)
        got = "CHECK(S) FAILED" in out
        print(("  PASS  " if got else "  FAIL  ") + f"caught: {name}")
        if got:
            caught += 1
        else:
            missed.append(name)
    shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"\n{caught} of {len(MUTATIONS)} defects caught")
    if missed:
        print("not caught:")
        for m in missed:
            print("  -", m)
    return 0 if (ok0 and not missed) else 1


if __name__ == "__main__":
    sys.exit(main())
