#!/usr/bin/env python3
"""Cross-check the manuscript against everything it rests on.

Six independent checks, each of which can fail on its own:

  A  style and journal format
  B  the reference list and its citations
  C  every number in the manuscript against the archived result files
  D  the same quantity stated in more than one place inside the manuscript
  E  the manuscript against the locked protocol and its amendments
  F  the numbers hardcoded in the figure scripts
  G  the submission package

Run:  python3 50_consistency_check.py [--verbose]
Exit code 0 when every check passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import unicodedata
from collections import defaultdict


# The manuscript, the protocol and the letter are not in the archive before
# publication; IR_DOCS says where they are and defaults to the directory above
# the archive, which is where they sit in the author's tree.
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
from paths import ROOT as IR_ROOT, RESULTS as IR_RESULTS, CODE as IR_CODE, DOCS as IR_DOCS
HOME = IR_DOCS
MS = f"{HOME}/manuscript_v4.md"
PROTOCOL = f"{HOME}/SupplementaryFile1_v2.md"
LETTER = f"{HOME}/CoverLetter_GigaScience.md"
RESULTS = IR_RESULTS
CODE_DIR = IR_CODE
FIGSCRIPTS = ([f"{HOME}/{n}" for n in (
    "make_figure1_concept.py", "make_figure3_v2.py", "make_figure4.py", "make_figure5.py",
    "make_figure6.py", "make_figure7.py", "make_figure8.py", "make_figure9.py",
    "make_figure11.py", "make_figureS1.py")] + [f"{CODE_DIR}/22_figures_1_2.py"])
FIGSCRIPTS = [p for p in FIGSCRIPTS if os.path.exists(p)] or [
    f"{CODE_DIR}/{n}" for n in sorted(os.path.basename(x) for x in
                                      glob.glob(f"{CODE_DIR}/make_figure*.py"))]
PACKAGE = f"{HOME}/submission_GigaScience"

FAIL: list[str] = []
NOTE: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def head(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


# --------------------------------------------------------------------- helpers
MINUS = "−"
SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
       "⁻": "-"}


def desuper(t: str) -> str:
    return "".join(SUP.get(c, c) for c in t)


def normalize(t: str) -> str:
    """Turn the manuscript's typography into something parseable."""
    t = t.replace(MINUS, "-").replace("–", "-").replace("×", "x")
    t = desuper(t)
    t = re.sub(r"\s*x\s*10\^?(-?\d+)", r"e\1", t)
    return t


def body_of(text: str) -> str:
    return re.sub(r"\n## References\n.*", "", text, flags=re.S)


def numbers_in(text: str) -> list[tuple[str, float]]:
    """Every numeric literal, as (as-written, value)."""
    out = []
    for m in re.finditer(r"(?<![\w.])[+-]?\d+(?:\.\d+)?(?:e[+-]?\d+)?(?![\w])", normalize(text)):
        tok = m.group(0)
        try:
            out.append((tok, float(tok)))
        except ValueError:
            pass
    return out


def walk_values(obj, path=""):
    """Every scalar in a nested JSON, with its key path."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk_values(v, f"{path}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk_values(v, f"{path}[{i}]")
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        yield float(obj), path


def load_archive() -> list[tuple[float, str]]:
    """Every numeric value in the archived results, JSON and markdown alike."""
    vals = []
    for dp, dn, fn in os.walk(RESULTS):
        for f in fn:
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, RESULTS)
            if f.endswith(".json"):
                try:
                    d = json.load(open(p))
                except Exception:
                    continue
                for v, path in walk_values(d):
                    vals.append((v, rel + path))
            elif f.endswith((".md", ".tsv", ".csv")):
                try:
                    t = open(p, encoding="utf-8", errors="ignore").read()
                except Exception:
                    continue
                for tok, v in numbers_in(t):
                    vals.append((v, rel))
    for extra in (f"{HOME}/B1_결과보고_2026-09-07.md", f"{HOME}/B3_결과보고_2026-09-07.md",
                  PROTOCOL):
        if os.path.exists(extra):
            t = open(extra, encoding="utf-8", errors="ignore").read()
            for tok, v in numbers_in(t):
                vals.append((v, os.path.basename(extra)))
    return vals


class Archive:
    """Archived values, indexed so a lookup is a dict hit rather than a scan."""

    def __init__(self, pairs):
        self.n = len(pairs)
        self.raw = {round(v, 6) for v, _ in pairs}
        self.by_dec = [defaultdict(list) for _ in range(7)]
        self.sci = defaultdict(list)
        for v, where in pairs:
            for d in range(7):
                self.by_dec[d][round(v, d)].append(where)
            if v != 0 and abs(v) < 1e-3:
                self.sci[f"{v:.1e}"].append(where)

    def matches(self, value: float) -> list:
        s = repr(float(value))
        if "e" in s:
            dec = 6
        else:
            dec = len(s.split(".")[1]) if "." in s else 0
        dec = min(dec, 6)
        hits = self.by_dec[dec].get(round(value, dec), [])
        if hits:
            return hits
        if value != 0 and abs(value) < 1e-3:
            return self.sci.get(f"{value:.1e}", [])
        return []


# ------------------------------------------------------------------- A. style
BRITISH = (r"\b(behaviour|colour|centre|analysed|analysing|labelled|modelling|normalised|"
           r"summarised|recognised|characterised|utilised|tumour|haematoxylin|favour|"
           r"neighbour|fibre|licence|organisation|generalised|minimised|maximised|"
           r"standardised|emphasised|hypothesised|programme|oesophag\w*|paediatric|artefact|"
           r"grey|acknowledgements|initialise\w*|prioritise\w*|randomise\w*|sensitised|"
           r"optimise\w*|categorise\w*|visualise\w*|ageing|defence|practise[sd]?|"
           r"metre|litre|sulphur\w*|catalogue|dialogue|traveller|labelling|signalling|"
           r"cancelled|modelled|totalled|fuelled|marvellous|sceptic\w*|manoeuvr\w*)\b")
SECTIONS = ["## Abstract", "## Background", "## Data Description", "## Analyses",
            "## Discussion", "## Potential implications", "## Methods",
            "## Availability of source code and requirements", "## Data Availability",
            "## List of abbreviations", "## Declarations", "## References"]


def check_style(ms: str) -> None:
    head("A. Style and journal format")
    body = body_of(ms)
    # The no-em-dash rule governs what I write. A published title is quoted as
    # published, so the reference list is out of scope; altering a title to
    # satisfy a house rule of my own would misquote it.
    _body = ms.split("## References", 1)[0]
    check(_body.count("\u2014") == 0,
          f"no em dash outside the reference list ({_body.count(chr(0x2014))} found)")
    hits = re.findall(BRITISH, body, re.I)
    check(not hits, f"American spelling outside reference titles ({sorted(set(hits))})")
    pl = re.findall(r"\b(we|our|us)\b", body, re.I)
    check(not pl, f"first person singular throughout ({sorted(set(x.lower() for x in pl))})")
    dbl = [m.group(0) for m in re.finditer(r"\b(\w+)\s+\1\b", ms, re.I)
           if m.group(1).lower() not in ("had", "that", "is", "s")]
    check(not dbl, f"no doubled word in the manuscript ({dbl})")
    pos = [ms.index(h) for h in SECTIONS if h in ms]
    check(len(pos) == len(SECTIONS) and pos == sorted(pos), "journal section order")
    ab = re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", ms, re.S).group(1)
    n = len(re.sub(r"\*\*(Background|Results|Conclusions)\.\*\*", "", ab).replace("**", "").split())
    check(n <= 250, f"abstract {n} words (limit 250)")
    abbr = [t for t in re.findall(r"\b([A-Z]{2,}[A-Za-z0-9-]*|[a-z][A-Z]{2,})\b", ab)
            if t not in ("RNA", "SHA-256")]
    check(not abbr, f"abstract carries no manuscript abbreviation ({abbr})")
    for m in re.finditer(r"\*\*Table (\d)\. (.*?)\*\*", ms):
        k = len(m.group(2).rstrip(".").split())
        check(k <= 15, f"Table {m.group(1)} title {k} words (limit 15)")
    for m in re.finditer(r"\*\*Figure (\d+)\.(.*?)(?=\n\n)", ms, re.S):
        k = len(m.group(2).split())
        check(k <= 300, f"Figure {m.group(1)} legend {k} words (limit 300)")
    for bad in ("PLACEHOLDER", "TO BE INSERTED", "TODO", "XXX"):
        check(bad not in ms, f"no '{bad}' left in the manuscript")
    # every abbreviation in the list is used, alphabetical, and defined at first use
    lst = re.search(r"## List of abbreviations\n\n(.*?)\n\n", ms, re.S).group(1)
    keys = [i.split(":")[0].strip() for i in lst.rstrip(".").split(";")]
    check(keys == sorted(keys, key=str.lower), "abbreviation list is alphabetical")
    rest = body.replace(lst, "")
    unused = [k for k in keys if not re.search(r"\b" + re.escape(k) + r"\b", rest)]
    check(not unused, f"every listed abbreviation is used in the text ({unused})")


# -------------------------------------------------------------- B. references
def expand(tok: str) -> list[int]:
    out = []
    for part in tok.split(","):
        part = part.strip().replace("–", "-")
        if "-" in part:
            a, b = part.split("-"); out.extend(range(int(a), int(b) + 1))
        elif part.isdigit():
            out.append(int(part))
    return out


def check_references(ms: str) -> None:
    head("B. References")
    i = ms.index("\n## References\n")
    body, refs = ms[:i], ms[i:]
    entries = {}
    dup = []
    for m in re.finditer(r"^(\d+)\. (.*)$", refs, re.M):
        n = int(m.group(1))
        if n in entries:
            dup.append(n)
        entries[n] = m.group(2).strip()
    check(not dup, f"no duplicate reference numbers ({dup})")
    check(sorted(entries) == list(range(1, len(entries) + 1)),
          f"reference numbers run 1 to {len(entries)} with no gap")
    order = []
    for m in re.finditer(r"\[([0-9][0-9,\s–-]*)\]", body):
        for n in expand(m.group(1)):
            if n not in order:
                order.append(n)
    check(order == sorted(order), f"citations appear in ascending order (first eight {order[:8]})")
    check(not (set(order) - set(entries)), f"every citation exists in the list "
                                           f"({sorted(set(order) - set(entries))})")
    check(not (set(entries) - set(order)), f"every listed reference is cited "
                                           f"({sorted(set(entries) - set(order))})")
    # each entry looks like a Vancouver reference
    for n, e in entries.items():
        ok = bool(re.search(r"\.\s\d{4};", e)) or bool(re.search(r"\.\s\d{4}\.", e))
        check(ok, f"reference {n} carries a year in Vancouver position")
    # a year in the future or before 1990 is almost certainly a typo
    for n, e in entries.items():
        yrs = [int(y) for y in re.findall(r"\b(19[5-9]\d|20[0-2]\d)\b", e)]
        check(all(1950 <= y <= 2026 for y in yrs) if yrs else True,
              f"reference {n} year is plausible ({yrs})")


# ------------------------------------------- C. numbers against the archive
# Numbers a reader would not expect to find in a result file.
IGNORE_EXACT = {
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.5, 0.05, 100.0, 1000.0, 10.0, 15.0, 50.0, 95.0,
    2026.0, 2025.0, 2024.0, 2023.0, 2022.0, 2021.0, 2020.0, 2019.0, 2018.0, 2017.0,
    2016.0, 2015.0, 2014.0, 2013.0, 2011.0, 2010.0, 1997.0, 256.0, 36.0, 6.0, 22.0,
    21.0, 19.0, 17.0, 26.0, 27.0, 30.0, 31.0, 29.0, 28.0, 9.0, 8.0, 7.0, 11.0, 12.0,
    13.0, 14.0, 16.0, 18.0, 20.0, 23.0, 24.0, 25.0,
}
# Regions of the manuscript that are not claims about results.
def analysable(ms: str) -> str:
    t = body_of(ms)
    t = t[t.index("## Abstract"):]            # drop the title page and contact details
    t = re.sub(r"\n## Availability of source code and requirements\n.*?\n## ", "\n## ", t, flags=re.S)
    t = re.sub(r"\n## Declarations\n.*?(\n## |\Z)", r"\1", t, flags=re.S)
    t = re.sub(r"\n## Additional files\n.*?\n## ", "\n## ", t, flags=re.S)
    t = re.sub(r"\[[0-9][0-9,\s–-]*\]", " ", t)          # citation markers
    t = re.sub(r"§6[a-z]+", " ", t)                        # amendment labels
    t = re.sub(r"GSE\d+|TCGA-\w+|PXD\d+|GPL\d+|HG-U\d+\w*|HT-HG-\w+|v?\d+\.\d+\.Hs|"
               r"GENCODE v\d+|C2:CGP|CYB5R\d|MTARC\d|CYP\w+|COQ\d|PDSS\d|FOXA\d|HNF\d\w*|"
               r"NR1H4|G6PC1?|CD\d+\w*|COL\d\w*|KRT\d+|SOX\d|AIFM\d|NQO\d|POR|ITGAM|"
               r"SHA-256|10\.5281/zenodo\.\d+|Figure \d+[a-d]?|Table \d+|S\d+", " ", t)
    return t


def check_anchored(ms: str, proto: str) -> None:
    """The check with teeth: each claim bound to the archived quantity it is about."""
    head("C1. Anchored claims (each bound to its archived quantity)")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "anchored", os.path.join(CODE_DIR, "52_anchored_claims.py"))
    anchored = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(anchored)
    n, fails, missing = anchored.run(ms)
    check(not fails, f"{n} anchored claims agree with the archive ({len(fails)} disagree)")
    for f in fails:
        NOTE.append("    " + f)
    check(not missing, f"every anchor still finds its sentence ({len(missing)} lost)")
    for m in missing:
        NOTE.append("    lost anchor: " + m)

    # Numbers quoted together in one sentence came from one place in the
    # pipeline, so they must appear together in one archived record. This
    # reaches the several hundred quantities that carry no hand-written anchor.
    # The manuscript now describes this checker, so those numbers must be the
    # checker's real numbers and must move when the checker does.
    WORD = {"ninety": 90, "a hundred and two": 102,
            "a hundred and three": 103, "a hundred and four": 104,
            "a hundred and five": 105, "a hundred and six": 106,
            "a hundred and seven": 107, "a hundred and eight": 108,
            "a hundred and nine": 109, "a hundred and ten": 110,
            "a hundred and eleven": 111, "a hundred and twelve": 112,
            "a hundred and thirteen": 113, "a hundred and fourteen": 114,
            "a hundred and fifteen": 115, "a hundred and sixteen": 116,
            "twenty-two": 22, "twenty-three": 23,
            "twenty-four": 24, "twenty-five": 25, "twenty-six": 26,
            "twenty-seven": 27, "twenty-eight": 28, "twenty-nine": 29,
            "thirty": 30, "thirty-one": 31, "thirty-two": 32, "thirty-three": 33, "thirty-four": 34,
            "thirty-five": 35, "thirty-six": 36, "thirty-seven": 37,
            "five": 5, "six": 6}
    m = re.search(r"([\w ]+?) claims are bound one by one", ms)
    check(m is not None and WORD.get(m.group(1).strip()) == n,
          f"the manuscript's stated anchor count matches the {n} anchors that ran")
    m = re.search(r"injects (\S+) known defects", ms)
    live_muts = len(re.findall(r'^    \("', open(
        os.path.join(CODE_DIR, "51_check_the_checker.py"),
        encoding="utf-8").read(), re.M))
    check(m is not None and WORD.get(m.group(1)) == live_muts,
          f"the manuscript's stated defect count matches the {live_muts} "
          f"mutations in the test")

    tn, tbad = anchored.run_tuples(ms)
    check(not tbad, f"{tn} co-occurring groups each sit in one archived record "
                    f"({len(tbad)} do not)")
    for t in tbad:
        NOTE.append("    " + t)

    # A few numbers come from a cited paper rather than from this archive.
    # They cannot be bound to a result file, so they are bound to the table
    # and row they were read from instead.
    ln, lbad = anchored.run_literature(ms)
    check(not lbad, f"{ln} literature-sourced values match their stated source "
                    f"({len(lbad)} do not)")
    for b in lbad:
        NOTE.append("    " + b)


def check_protocol_numbers(proto: str) -> None:
    """The protocol quotes results too, and they must match the same archive."""
    head("E2. Numbers quoted together in the protocol")
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "anchored", os.path.join(CODE_DIR, "52_anchored_claims.py"))
    anchored = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(anchored)
    n, bad = anchored.run_tuples(proto, whole=True)
    check(not bad, f"{n} co-occurring groups in the protocol sit in one "
                   f"archived record ({len(bad)} do not)")
    for t in bad:
        NOTE.append("    " + t)


def check_numbers(ms: str, archive, verbose: bool) -> None:
    head("C2. Every other number has some archived support")
    text = analysable(ms)
    seen = {}
    for tok, val in numbers_in(text):
        if val in IGNORE_EXACT or abs(val) >= 10000:
            continue
        if float(val).is_integer() and abs(val) <= 2500:
            continue           # counts are checked in D and E, not here
        seen.setdefault(round(val, 10), tok)
    unmatched = []
    for val, tok in sorted(seen.items()):
        if archive.matches(val):
            continue
        if archive.matches(round(val / 100.0, 8)):   # a percentage of an archived fraction
            continue
        if archive.matches(round(val * 100.0, 8)):
            continue
        unmatched.append(tok)
    check(not unmatched,
          f"every non-integer quantity has an archived source "
          f"({len(unmatched)} without one)")
    if unmatched:
        for tok in unmatched:
            ctx = re.search(r".{70}" + re.escape(tok.lstrip("+")) + r".{40}", normalize(text))
            NOTE.append(f"    no archived value for {tok}: "
                        f"...{ctx.group(0).strip() if ctx else ''}...")
    if verbose:
        print(f"        {len(seen)} distinct quantities checked")


# ------------------------------------------- D. internal consistency of values
def check_internal(ms: str) -> None:
    head("D. The same quantity in more than one place")
    body = body_of(ms)
    # a quantity introduced with a label must carry the same value everywhere
    pairs = [
        ("liver median retention", [r"median of \*\*0\.758\*\*", r"median 0\.758", r"medians 0\.758"]),
        ("lung median retention", [r"0\.576 in lung", r"medians 0\.758, 0\.576"]),
        ("kidney median retention", [r"0\.467 in kidney", r"0\.576, 0\.467"]),
        ("REDUCTION joint", [r"REDUCTION retained \*\*89%\*\*", r"\*\*89%\*\*"]),
        ("DRAIN joint", [r"DRAIN's \*\*31%\*\*", r"31\.0%"]),
    ]
    for name, pats in pairs:
        found = [p for p in pats if re.search(p, body)]
        check(bool(found), f"{name} is stated somewhere in a recognised form")
    # counts that appear twice must agree
    for label, pat in [("evaluable liver signatures", r"(\d+) (?:published )?liver signatures"),
                       ("pairs in GSE14520", r"GSE14520[^.]{0,40}?(\d{3}) pairs"),
                       ("amendments", r"(thirty|twenty-nine|twenty-eight) (?:timestamped )?amendments")]:
        vals = set(re.findall(pat, body))
        check(len(vals) <= 1, f"{label} stated consistently ({sorted(vals)})")
    # the abstract must not contradict the body on the headline numbers
    ab = re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", ms, re.S).group(1)
    for tok in ("0.758", "0.576", "0.467", "2.345", "0.004"):
        if tok in ab:
            check(tok in body_of(ms).replace(ab, ""),
                  f"abstract value {tok} also appears in the body")


# ------------------------------------------------- E. manuscript vs protocol
def check_protocol(ms: str, proto: str) -> None:
    head("E. Manuscript against the locked protocol")
    # An amendment's own heading says which tissue or subject it governs. The
    # manuscript must not attach an amendment to the wrong one: Table 1 called
    # lung the second tissue and kidney the third, which is the reverse of
    # what 6x and 6y say, and nothing here caught it.
    for sec, word in (("6x", "second"), ("6y", "third")):
        m = re.search(r"^#{1,3} §?" + sec + r"\.\s*(.+)$", proto, re.M)
        if m:
            check(word in m.group(1).lower(),
                  f"the protocol calls \u00a7{sec} the {word} tissue")
    pairs_t = re.findall(r"\|\s*(Second|Third) tissue \(\u00a7(6[a-z]+)\)", ms)
    check(len({w for w, _ in pairs_t}) == len(pairs_t),
          f"Table 1 gives each tissue ordinal once ({[w for w, _ in pairs_t]})")
    check(len({a for _, a in pairs_t}) == len(pairs_t),
          f"Table 1 cites each tissue amendment once ({[a for _, a in pairs_t]})")
    for word, sec in [(w.lower(), a) for w, a in pairs_t]:
        m = re.search(r"^#{1,3} §?" + sec + r"\.\s*(.+)$", proto, re.M)
        check(m is not None and word in m.group(1).lower(),
              f"Table 1's \"{word} tissue\" cites the amendment that calls itself that")
    # Every MSigDB set the manuscript names in backticks must either appear in
    # the protocol, or be declared in the manuscript as a substitution.
    named = set(re.findall(r"`([A-Z][A-Z0-9_]{6,})`", ms))
    for g in sorted(named):
        if g in proto:
            continue
        i = ms.find("`" + g + "`")
        window = ms[max(0, i - 700):i + 700]
        check("substitut" in window or "not pre-specified" in window,
              f"the set {g}, absent from the protocol, is declared a substitution")
    cited = sorted(set(re.findall(r"§(6[a-z]+)", ms)))
    have = set(re.findall(r"§?(6[a-z]+)[\.\s]", proto))
    missing = [a for a in cited if a not in have]
    check(not missing, f"every amendment the manuscript cites exists in the protocol ({missing})")
    letters = sorted(set(re.findall(r"^#{1,3} §?(6[a-z]+)[\.\s]", proto, re.M)),
                     key=lambda x: (len(x), x))
    n_amend = len(letters)
    stated = re.search(r"(Thirty|Twenty-nine|Twenty-eight|Thirty-one)\s+amendments \(§6a to §(6[a-z]+)\)", ms)
    check(stated is not None, "the manuscript states an amendment count and range")
    if stated:
        words = {"Twenty-eight": 28, "Twenty-nine": 29, "Thirty": 30, "Thirty-one": 31}
        check(words[stated.group(1)] == n_amend,
              f"amendment count matches the protocol ({stated.group(1)} vs {n_amend} headings)")
        check(stated.group(2) == letters[-1] if letters else False,
              f"last amendment letter matches ({stated.group(2)} vs {letters[-1] if letters else '?'})")
    # the hash stem
    stem = re.search(r"`(7a2bf934[^`]*)`", ms)
    check(stem is not None and stem.group(1) in proto,
          "the lock digest quoted in the manuscript appears in the protocol")
    # thresholds the manuscript attributes to an amendment must appear in that amendment
    for thr, amend in [("0.30", "6z"), ("0.15", "6w"), ("0.20", "6y"), ("0.70", "6u"), ("−0.30", "6ab")]:
        blk = re.search(rf"#{{1,3}} §?{amend}\b.*?(?=\n#{{1,2}} |\Z)", proto, re.S)
        if blk:
            t = normalize(blk.group(0))
            check(normalize(thr).lstrip("-") in t,
                  f"threshold {thr} appears inside §{amend} of the protocol")
    # a claim about who wrote what when
    WORDNUM = {"three": 3, "four": 4, "five": 5, "six": 6}
    m = re.search(r"(\w+) amendments in all are not pre-registered in any sense", ms)
    n_stated = WORDNUM.get(m.group(1).lower(), 0) if m else 0
    named = set(re.findall(r"§6(?:q|r|ac|ad|ae)\b",
                           ms[m.start():m.start() + 900])) if m else set()
    # every amendment whose own text disclaims pre-registration must be named
    disclaim = set()
    for blk in re.split(r"\n# Amendment ", proto)[1:]:
        tag = re.match(r"(§6[a-z]+)", blk)
        if tag and "not pre-registered" in blk[:1200].lower().replace(
                "nothing here is pre-registered", "not pre-registered"):
            disclaim.add(tag.group(1))
    for tag in ("§6q", "§6r"):
        if tag not in named:
            disclaim.add(tag)
    check(m is not None and n_stated == len(named) and disclaim <= named,
          f"the manuscript names all {n_stated} amendments that are not "
          f"pre-registered (names {sorted(named)}; protocol disclaims "
          f"{sorted(disclaim)})")


# ---------------------------------------------------- F. figures vs the text
def check_figures(ms: str, archive) -> None:
    head("F. Numbers hardcoded in the figure scripts")
    body = normalize(body_of(ms))
    ms_vals = {round(v, 6) for _, v in numbers_in(body)}
    arch_vals = archive.raw
    bad = []
    for p in FIGSCRIPTS:
        if not os.path.exists(p):
            NOTE.append(f"    figure script missing: {p}")
            continue
        src = open(p, encoding="utf-8").read()
        # only literal lists that look like plotted data, not layout coordinates
        for m in re.finditer(r"^\s*(g?vals|tvals|gvals|d1)\s*=\s*\[(.*?)\]", src, re.M | re.S):
            for tok in re.findall(r"[+-]?\d+\.\d+", m.group(2)):
                v = round(float(tok), 6)
                if v in ms_vals or v in arch_vals:
                    continue
                if any(abs(v - a) < 5e-4 for a in arch_vals | ms_vals):
                    continue
                bad.append(f"{os.path.basename(p)}: {tok}")
        # percentage labels drawn on bars
        for m in re.finditer(r"^\s*(g?ret|tret)\s*=\s*\[(.*?)\]", src, re.M):
            for tok in re.findall(r'"(\d+(?:\.\d+)?)%"', m.group(2)):
                v = float(tok) / 100
                if any(abs(v - a) < 6e-3 for a in arch_vals | ms_vals):
                    continue
                bad.append(f"{os.path.basename(p)}: label {tok}%")
    check(not bad, f"every plotted literal is in the text or the archive ({bad})")


# ------------------------------------------------------------- G. the package
def check_cover_letter(ms: str) -> None:
    """The cover letter is a deliverable too, and it drifts.

    An audit found it claiming three non-pre-registered amendments where the
    manuscript says four, and asserting a purity benchmark the manuscript
    declines to claim. Nothing checked it, so these checks exist now."""
    head("H. The cover letter")
    path = f"{HOME}/CoverLetter_GigaScience.md"
    if not os.path.exists(path):
        check(False, "the cover letter exists")
        return
    cl = open(path, encoding="utf-8").read()
    check("\u2014" not in cl, f"no em dash in the cover letter "
                               f"({cl.count(chr(8212))} found)")
    dbl_cl = re.findall(r"\b(\w+)\s+\1\b", cl, re.I)
    check(not dbl_cl, f"no doubled word in the cover letter ({dbl_cl})")
    plural = re.findall(r"\b(?:we|our|us)\b", cl, re.I)
    check(not plural, f"first person singular in the cover letter ({plural})")
    brit = [w for w in ("organis", "colour", "behaviour", "centre", "recognise",
                        "summarise", "utilis", "tumour", "labelled", "modelling",
                        "characteris", "normalis", "analyse", "analysing", "favour", "fibre",
                        "programme", "practise", "licence")
            if w in cl.lower()]
    check(not brit, f"American spelling in the cover letter ({brit})")
    check("Sincerely," in cl and "Yours sincerely" not in cl,
          "the cover letter closes with the American form")
    ms_title = re.search(r"^# (.+)$", ms, re.M).group(1).strip()
    check(ms_title in cl, "the cover letter carries the manuscript's exact title")
    # The letter states a word count. It went stale twice, so it is computed
    # here from the manuscript rather than trusted, on the same definition:
    # Background through the end of Additional files, table rows excluded,
    # figure legends included.
    body = ms.split("## Background", 1)[1].split("## References", 1)[0]
    body = "\n".join(l for l in body.split("\n") if not l.strip().startswith("|"))
    txt_ = re.sub(r"`[^`]*`", " ", body)
    txt_ = re.sub(r"[*_#|>-]", " ", txt_)
    real_wc = len([w for w in txt_.split() if any(c.isalnum() for c in w)])
    m_wc = re.search(r"about ([\d,]+) words", cl)
    stated = int(m_wc.group(1).replace(",", "")) if m_wc else 0
    check(m_wc is not None and abs(stated - real_wc) <= 400,
          f"the letter's word count ({stated:,}) is within 400 of the "
          f"manuscript's ({real_wc:,})")

    # counts the letter states must be the manuscript's counts
    # Derive the anchor wording from the manuscript rather than fixing it here,
    # so that adding an anchor cannot leave the letter quietly behind.
    anchor_word = re.search(r"([\w ]+?) claims are bound one by one", ms)
    pairs = [("thirty-one", "amendments")]
    if anchor_word:
        pairs.append((anchor_word.group(1).strip(), "anchored claims"))
    for word, what in pairs:
        check(word in cl, f"the cover letter states {what} as '{word}'")
    m = re.search(r"(\w+) of them state in their own text that they are not "
                  r"pre-registered", cl)
    ms_m = re.search(r"(\w+) amendments in all are not pre-registered", ms)
    check(m is not None and ms_m is not None and m.group(1) == ms_m.group(1),
          "the cover letter and the manuscript agree on how many amendments "
          "are not pre-registered")
    m = re.search(r"injecting (\S+) known defects", cl)
    ms_m = re.search(r"injects (\S+) known defects", ms)
    check(m is not None and ms_m is not None and m.group(1) == ms_m.group(1),
          "the cover letter and the manuscript agree on the injected-defect count")
    # every DOI, repository URL and identifier must match the manuscript
    for tok in ("10.5281/zenodo.22658669",
                "https://github.com/drswnam69-debug/identity-retention",
                "0000-0001-5562-4775", "drswnam@catholic.ac.kr"):
        check(tok in cl and tok in ms or tok in cl and tok not in ms and
              tok.startswith("0000") is False or tok in cl,
              f"the cover letter carries {tok}")
    # claims the manuscript must not contradict
    check("not what adjusting for tumor purity does" not in cl,
          "the cover letter does not claim a purity benchmark the manuscript "
          "declines to claim")


def check_alt_text(ms: str) -> None:
    """Alt text must describe the figure it is numbered for.

    The alt text for Figures 8 and 9 was swapped: each described the other's
    panels. Nothing caught it, because both were well-written descriptions of a
    figure in this paper. The check compares distinctive words shared with each
    legend and requires the numbered pairing to be the best available match."""
    head("I. Figure alt text")
    path = f"{HOME}/submission_GigaScience/06_Figure_alt_text.md"
    if not os.path.exists(path):
        check(False, "the alt-text file exists")
        return
    alt = open(path, encoding="utf-8").read()
    legends, alts = {}, {}
    for m in re.finditer(r"\*\*Figure (\d+)\..*?\*\*(.*?)(?=\n\n)", ms, re.S):
        legends[int(m.group(1))] = m.group(2)
    for m in re.finditer(r"\*\*Figure (\d+)\.\*\*(.*?)(?=\n\n|\Z)", alt, re.S):
        alts[int(m.group(1))] = m.group(2)
    check(set(alts) == set(legends),
          f"alt text covers exactly the {len(legends)} numbered figures "
          f"(missing {sorted(set(legends) - set(alts))}, "
          f"extra {sorted(set(alts) - set(legends))})")
    stop = set("the a an and or of in to for with each its is are two three four "
               "five panels panel across against between per from that this on by "
               "as at one same their them it not no number values value shown "
               "with 50 percent".split())
    def words(t):
        return {w for w in re.findall(r"[a-z]{4,}", t.lower())} - stop
    for n in sorted(set(alts) & set(legends)):
        own = len(words(alts[n]) & words(legends[n]))
        scores = sorted(((len(words(alts[n]) & words(legends[k])), k)
                         for k in legends), reverse=True)
        best = scores[0]
        check(own > 0 and best[1] == n,
              f"Figure {n} alt text matches Figure {n}'s legend better than any "
              f"other ({own} shared terms; best is Figure {best[1]} with "
              f"{best[0]})")


def check_registered_rules(ms: str, proto: str) -> None:
    """Every registered decision rule must be answered in the manuscript.

    Four compliance failures reached the pre-submission audit: a registered
    verdict inverted in a heading, two pre-specified secondary analyses run and
    archived but never reported, and an outcome branch that required a
    statement in the abstract and did not get one. None of them was a wrong
    number, so nothing here could see them. This finds each place the protocol
    fixes an outcome in advance and requires that the amendment it belongs to
    is discussed in the manuscript, so a rule cannot be silently skipped.
    It cannot judge whether the verdict is right; it can insist the rule was
    not forgotten."""
    head("L. Registered decision rules are answered")
    # locate every rule statement and attribute it to its amendment
    rules = []
    cur = "the locked plan"
    for line in proto.split("\n"):
        m = re.match(r"^#{1,3}\s+§?(6[a-z]+)[\.\s]", line)
        if m:
            cur = "§" + m.group(1)
        if re.search(r"\*\*(Decision rule fixed in advance|Outcomes accepted "
                     r"in advance|Rule fixed in advance|Decision rule)", line):
            rules.append(cur)
    seen = {}
    for r in rules:
        seen[r] = seen.get(r, 0) + 1
    check(len(rules) >= 15,
          f"the protocol still states its decision rules ({len(rules)} found)")
    missing = [r for r in sorted(seen) if r.startswith("§") and r not in ms]
    check(not missing,
          f"every amendment that fixes an outcome is discussed in the "
          f"manuscript ({len(missing)} are not: {missing})")
    # a rule whose amendment says "in the abstract" must reach the abstract
    abstract = ms.split("## Abstract", 1)[1].split("\n## ", 1)[0].lower()
    for m in re.finditer(r"§?(6[a-z]+)", ""):
        pass
    blocks = re.split(r"^#{1,3}\s+§?(6[a-z]+)[\.\s]", proto, flags=re.M)
    for i in range(1, len(blocks) - 1, 2):
        sec, body = blocks[i], blocks[i + 1]
        if re.search(r"in the abstract", body, re.I):
            check("§" + sec in ms,
                  f"§{sec}, which fixes an outcome for the abstract, is in the "
                  f"manuscript")
    # the four amendments that route an outcome to the abstract, and the words
    # in the abstract that answer each. §6w's branch went unhonored for three
    # revisions because nothing checked the abstract itself.
    for sec, probe, what in (("6s", "premise", "the premise failure"),
                             ("6w", "replicate", "the non-replication"),
                             ("6ae", "random covariate", "the negative control")):
        if "§" + sec in proto:
            check(probe in abstract,
                  f"the abstract carries {what} that §{sec} routes to it")


def check_figure_glyphs() -> None:
    """Every character a figure script draws must exist in the figure font.

    Liberation Sans has no superscript minus and no superscript zero. Written
    as literal characters they render as empty boxes, and "9.0 x 10 to the
    minus ten" printed as "9.0 x 10[]1[]" in a submitted figure for three
    revisions without anything noticing, because nothing here had ever looked
    at a glyph. Mathtext is the fix and this is the guard."""
    head("K. Figures use only characters the font can draw")
    try:
        from matplotlib import font_manager
        from fontTools.ttLib import TTFont
    except Exception as e:                                  # pragma: no cover
        NOTE.append(f"    glyph check skipped ({e})")
        return
    path = font_manager.findfont("Liberation Sans")
    cov = set()
    for t in TTFont(path, fontNumber=0)["cmap"].tables:
        cov |= set(t.cmap.keys())
    bad_total = {}
    for f in sorted(glob.glob(f"{HOME}/make_figure*.py")):
        src = open(f, encoding="utf-8").read()
        # literal text the script draws, plus any \uXXXX escape it writes
        chars = set()
        for a, b in re.findall(r'"([^"\\]*)"|\'([^\'\\]*)\'', src):
            chars |= set(a or "") | set(b or "")
        for m in re.findall(r"\\u([0-9a-fA-F]{4})", src):
            chars.add(chr(int(m, 16)))
        bad = sorted(c for c in chars if ord(c) > 127 and ord(c) not in cov)
        if bad:
            bad_total[os.path.basename(f)] = [f"U+{ord(c):04X}" for c in bad]
    check(not bad_total,
          f"every character drawn by a figure script is in {os.path.basename(path)} "
          f"({len(bad_total)} script(s) use one that is not)")
    for k, v in bad_total.items():
        NOTE.append(f"    {k}: {', '.join(v)}")


def check_figure_deposit() -> None:
    """Every figure in the manuscript must be regenerable from the deposit.

    An audit found that only two of the ten figure scripts were deposited, and
    that one of them read its input from an upload directory outside the
    archived tree, so the submitted figures could not be rebuilt from the
    deposit alone. That contradicted the Availability section."""
    head("J. The figures are regenerable from the deposit")
    import glob
    scripts = sorted(os.path.basename(p) for p in glob.glob(f"{HOME}/make_figure*.py"))
    scripts.append("sync_figures.py")
    for d, label in ((f"{HOME}/release/identity-retention", "repository"),
                     (f"{HOME}/SF2build/SupplementaryFile2", "Supplementary File 2")):
        missing = [f for f in scripts if not os.path.exists(f"{d}/code/{f}")]
        check(not missing, f"the {label} carries every figure script "
                           f"({len(missing)} missing: {missing})")
    outside = []
    for p_ in glob.glob(f"{HOME}/make_figure*.py") + glob.glob(f"{CODE_DIR}/*.py"):
        for ln in open(p_, encoding="utf-8"):
            if ln.lstrip().startswith("#"):
                continue
            if "user-data/upl" + "oads" in ln:
                outside.append(f"{os.path.basename(p_)}: {ln.strip()[:70]}")
    check(not outside, f"no analysis or figure script reads outside the archive "
                       f"({len(outside)} do)")
    for o in outside:
        NOTE.append("    " + o)


def check_metadata(ms: str) -> None:
    """The submission form's title and abstract must be the manuscript's.

    This file is typed into the journal's own boxes. It was written by hand and
    went stale when the abstract was revised, which nothing caught, so it is now
    generated and verified."""
    head("G2. The submission metadata file matches the manuscript")
    path = f"{HOME}/submission_GigaScience/05_title_abstract_keywords.txt"
    if not os.path.exists(path):
        check(False, "05_title_abstract_keywords.txt exists")
        return
    meta = open(path, encoding="utf-8").read()
    title = re.search(r"^# (.+)$", ms, re.M).group(1).strip()
    check(title in meta, "the title matches the manuscript")
    ab = re.search(r"## Abstract\n\n(.*?)\n\n\*\*Keywords", ms, re.S).group(1)
    plain = re.sub(r"\*\*(.+?)\.\*\*", r"\1.", ab).replace("**", "").replace("*", "")
    miss = [p.strip()[:60] for p in plain.strip().split("\n\n")
            if p.strip() and p.strip() not in meta]
    check(not miss, f"every abstract paragraph matches the manuscript ({len(miss)} do not)")
    for m in miss:
        NOTE.append("    stale in 05_: " + m + " ...")
    kw = re.search(r"\*\*Keywords\.\*\*\s*(.+)", ms).group(1).strip()
    check(kw in meta, "the keyword list matches the manuscript")


def check_package(ms: str) -> None:
    head("G. The submission package")
    if not os.path.isdir(PACKAGE):
        check(False, f"package directory {PACKAGE} exists")
        return
    have = set(os.listdir(PACKAGE))
    need = ["00_CoverLetter.docx", "00_CoverLetter.pdf", "01_Manuscript.docx",
            "02_Manuscript_typeset.pdf", "03_SupplementaryFile1_Preregistration.zip",
            "04_SupplementaryFile2_code_and_data.zip", "05_title_abstract_keywords.txt",
            "06_Figure_alt_text.md", "07_Figure_legends.docx"]
    n_fig = len(re.findall(r"\*\*Figure (\d+)\.", ms))
    need += [f"Figure{i}.tif" for i in range(1, n_fig + 1)] + ["FigureS1.tif"]
    for f in need:
        check(f in have, f"package contains {f}")
    try:
        from PIL import Image
        for f in sorted(x for x in have if x.endswith(".tif")):
            p = os.path.join(PACKAGE, f)
            w, h = Image.open(p).size
            mm_h = h / 300 * 25.4
            check(w == 2055 and mm_h <= 225 and os.path.getsize(p) < 10e6,
                  f"{f}: {w}px wide, {mm_h:.0f} mm tall, {os.path.getsize(p)/1e6:.2f} MB")
    except ImportError:
        NOTE.append("    Pillow not available; figure dimensions not checked")
    # the DOCX must carry every figure, every table and every heading
    try:
        import zipfile
        z = zipfile.ZipFile(os.path.join(PACKAGE, "01_Manuscript.docx"))
        doc = z.read("word/document.xml").decode()
        txt = re.sub(r"<[^>]+>", "", doc)
        check(len([n for n in z.namelist() if n.startswith("word/media/")]) == n_fig,
              f"DOCX embeds {n_fig} figures")
        check(doc.count("<w:tbl>") == len(re.findall(r"\*\*Table (\d)\.", ms)),
              "DOCX carries every table")
        check('w:line="480"' in z.read("word/styles.xml").decode(), "DOCX is double spaced")
        for h in SECTIONS:
            check(h[3:] in txt, f"DOCX carries the {h[3:]} section")
    except Exception as e:
        check(False, f"DOCX readable ({e})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    ms = open(MS, encoding="utf-8").read()
    proto = open(PROTOCOL, encoding="utf-8").read()
    archive = Archive(load_archive())
    print(f"manuscript {os.path.basename(MS)}, {len(ms.split())} words; "
          f"{archive.n} archived values from {RESULTS}")
    check_style(ms)
    check_references(ms)
    check_anchored(ms, proto)
    check_numbers(ms, archive, a.verbose)
    check_internal(ms)
    check_protocol(ms, proto)
    check_protocol_numbers(proto)
    check_figures(ms, archive)
    check_package(ms)
    check_metadata(ms)
    check_cover_letter(ms)
    check_alt_text(ms)
    check_registered_rules(ms, proto)
    check_figure_glyphs()
    check_figure_deposit()
    if NOTE:
        print("\nDetail")
        print("------")
        for n in NOTE:
            print(n)
    print("\n" + ("ALL CHECKS PASS" if not FAIL else f"{len(FAIL)} CHECK(S) FAILED"))
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
