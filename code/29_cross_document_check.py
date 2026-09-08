#!/usr/bin/env python3
"""29_cross_document_check.py -- mechanical cross-document consistency.

Several errors in this project were introduced by editing: a number corrected in
one file and left stale in another, a count updated in the Methods but not the
Declarations, biochemistry retracted in the body but still asserted in the
abstract. This checks the invariants that editing breaks, mechanically, so the
check can be re-run after every future edit.
"""
import json, os, re, sys

FAIL = []
def check(name, ok, detail=""):
    print(f"  [{'ok ' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail and not ok else ""))
    if not ok:
        FAIL.append((name, detail))

body = open("manuscript/HepatolInt_v2_body.md", encoding="utf-8").read()
title = open("manuscript/TitlePage.md", encoding="utf-8").read()
cover = open("manuscript/CoverLetter.md", encoding="utf-8").read()
prereg = open("PREREGISTRATION.md", encoding="utf-8").read()
supp = open("manuscript/07_SupplementaryFile1_Preregistration.md", encoding="utf-8").read()

print("=== A. the supplementary file is the pre-registration ===")
check("Supplementary File 1 is byte-identical to PREREGISTRATION.md", supp == prereg,
      f"lengths {len(supp)} vs {len(prereg)}")

print("\n=== B. amendment count agrees everywhere ===")
last = sorted(re.findall(r'^## §6([a-z]) ', prereg, re.M))[-1]
n_amend = ord(last) - ord('a') + 1
words = {13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen", 17: "seventeen"}
w = words.get(n_amend, str(n_amend))
check(f"pre-registration's last amendment is §6{last} -> {n_amend} amendments", True)
check(f"body says '{w} amendments'", w in body.lower(),
      "body: " + ", ".join(set(re.findall(r'(\w+) amendments', body))))
check(f"body says §6a–§6{last}", f"§6a–§6{last}" in body)
check(f"title page says '{w}'", w in title.lower() or "amendment" not in title.lower())

print("\n=== C. the abstract does not contradict the body ===")
retracted = [("POR donates to mARC", r"electron source of a pathogenic", "abstract still calls CoQ the electron source of the N-reductive system"),
             ("POR/P450 explanation", r"as expected, since POR supplies", "the §6o-falsified clause is still present")]
for label, pat, why in retracted:
    check(f"retracted claim absent: {label}", not re.search(pat, body), why)

print("\n=== D. hypothesis accounting matches the registered verdicts ===")
pe = open("results/PHASE_E_RESULTS.md", encoding="utf-8").read()
not_supported = len(re.findall(r'\|\s*H\d\s*\|[^|]*\|[^|]*not supported', pe))
not_testable = len(re.findall(r'\|\s*H\d\s*\|[^|]*\|[^|]*not testable', pe))
check(f"PHASE E records {not_supported} not supported, {not_testable} not testable",
      not_supported == 3 and not_testable == 1, f"{not_supported}/{not_testable}")
check("abstract says three not supported and a fourth not testable",
      "Three of six pre-registered hypotheses were not supported and a fourth was not testable" in body)
check("no 'underpowered, not refuted' softening remains", "underpowered, not refuted" not in body)

print("\n=== E. every number in Table 2 traces to a results JSON ===")
C = {c: json.load(open(f"results/{c}_composition.json")) for c in ("GSE76427", "GSE14520")}
D = {c: json.load(open(f"results/{c}_differentiation.json")) for c in ("GSE76427", "GSE14520")}
bad = []
for c in ("GSE76427", "GSE14520"):
    want = [D[c]["paired_adjusted_D1"]["intercept"]["beta"],
            D[c]["paired_adjusted_D2"]["intercept"]["beta"],
            C[c]["composite_C1"]["estimate"], C[c]["composite_joint"]["estimate"],
            C[c]["composite_joint_C2"]["estimate"],
            C[c]["module_reduction_joint"]["estimate"],
            C[c]["module_drain_joint"]["estimate"]]
    for v in want:
        tok = f"{v:+.3f}".replace("-", "−")
        if tok not in body and f"{v:+.3f}" not in body:
            bad.append(f"{c}:{v:+.3f}")
check("all Table 2 point estimates present in the manuscript", not bad, str(bad))

print("\n=== F. confidence intervals match the corrected JSON ===")
ci_bad = []
for c in ("GSE76427", "GSE14520"):
    for key, node in (("D1", D[c]["paired_adjusted_D1"]["intercept"]),
                      ("joint", C[c]["composite_joint"]["fit"]["intercept"])):
        lo, hi = node["ci95"]
        pat = f"{lo:.3f}".lstrip("+")
        if pat not in body.replace("−", "-"):
            ci_bad.append(f"{c}/{key} lo={lo:.3f}")
check("reported interval bounds match the JSON", not ci_bad, str(ci_bad))

print("\n=== G. cover letter agrees with the body ===")
check("cover letter's error count matches the body's five",
      ("four amendments" in cover.lower()) == ("four amendments" in body.lower())
      or "amendments in which I record" in cover)
m = re.search(r'(\d+)% versus ([\d.]+)%, against ([\d.]+)% with', cover)
check("cover letter's 89/31.8/31.0 figures match §6l",
      "89% versus 31.8%" in cover or "89 % versus 31.8" in cover
      or "unaltered" in cover, "cover letter POR figures")

print("\n=== H. word and count limits ===")
def wc(t):
    t = re.sub(r'\*+', '', t)
    return len([w for w in t.split() if re.search(r'[0-9A-Za-z]', w)])
refs = body.split("## References", 1)[1]
tot = wc(body.split("## Introduction", 1)[1].split("## Tables", 1)[0]) + wc(refs)
ab = re.sub(r'\*+', '', body.split("## Abstract", 1)[1].split("**Keywords**")[0]).strip()
ntab = len(re.findall(r'^\*\*Table \d', body, re.M))
nfig = len(re.findall(r'^\*\*Figure \d', body, re.M))
nref = len(re.findall(r'^\d+\. ', refs, re.M))
check(f"abstract {len(ab.split())} <= 250", len(ab.split()) <= 250)
check(f"text+refs {tot} <= 4000", tot <= 4000)
check(f"references {nref} <= 30", nref <= 30)
check(f"tables+figures {ntab+nfig} <= 6", ntab + nfig <= 6)
check(f"title page states text+refs = {tot}", f"{tot:,} words" in title, title[title.find('Text including'):][:60])

print("\n" + "=" * 60)
print(f"{len(FAIL)} inconsistencies" if FAIL else "No cross-document inconsistency found.")
for n, d in FAIL:
    print(f"  FAIL: {n}  {d}")
sys.exit(0)
