# Target journals — RSI manuscript (방향 2)

> **개정 2026-08-26 (2차).** 이 문서의 최우선 권고였던 “Tier 3 위로 가려면 탈분화
> 교란부터 해결하라”는 **완료되었고 GSE14520(213쌍, HBV, 3번째 플랫폼)에서 재현되었다**
> (PREREGISTRATION §6h/§6i). 그 결과 순위와 부채 목록이 바뀌었다. 개정판 전문은
> `claude/방향2_투고전략_v2.md`와 아티팩트 “어디에 보낼 것인가”에 있다.
>
> 아래 원문은 **기록으로 남긴다** — 어떤 판단이 어떤 시점의 데이터에 근거했는지가
> 이 프로젝트의 자산이기 때문이다.

---


Prepared 2026-08-26. Facts below come from journal, society or publisher pages,
or DOAJ. **Anything I could not confirm at a primary source is marked
"unverified"** — aggregator sites (manusights, journalmetrics, askbisht,
journalsearches) were treated as unusable. Wiley and Karger pages blocked this
session, so Liver International and Liver Cancer are largely unverified.

---

## 1. What we are actually selling

This changed with Phase E. It is **not** a negative paper.

**Assets**
- Two pre-registered hypotheses **supported**: H1 replicated across three
  cohorts (pooled ρ = 0.372, I² = 0%, n = 392) and H2 replicated with a large
  effect (Cliff's δ = +0.76; paired P = 8.1 × 10⁻⁸ over 52 matched pairs).
- A **SHA-256-locked protocol with timestamped amendments**. Essentially no
  reanalysis paper in this field has this. It is the single most differentiating
  feature and is currently under-used.
- The nulls are *disciplined*, not embarrassing: H4 reports a bounded null with
  a working model (BCLC C/D HR 4.56 in the same fit).

**Liabilities, in the order reviewers will raise them**
1. **No independent validation and no wet-lab work.** Fatal at *J Hepatol* and
   *Hepatology*; negotiable everywhere else.
2. **The dedifferentiation confound on H2.** mARC1/mARC2/POR are hepatocyte
   differentiation genes and HCC is dedifferentiated. A competent reviewer will
   ask this in the first round and we currently **cannot answer it** — the
   42-gene panel has no ALB/CYP2E1/HNF4A.
3. Zero newly generated data.
4. H4 powered on 23 events.

> **Strong recommendation: fix liability #2 before submitting anywhere above
> Tier 3.** Re-preparing GSE76427 with a genome-wide probe map and refitting H2
> with a differentiation covariate is roughly a day of work and is the
> difference between "interesting index" and "reviewer 2 kills it". It would
> also unlock H6 in both array cohorts.

---

## 2. The shortlist

| | Journal | IF (verified) | OA | APC | Short format | Fit |
|---|---|---|---|---|---|---|
| **1** | **Hepatology International** (APASL/Springer) | **7.8 (2025)** | hybrid | **free via subscription route** (OA £3,290) | Original Art. capped at 4,000 words *incl. refs* | best cost-adjusted; **7-day median first decision** |
| **2** | **JHEP Reports** (EASL/Elsevier) | **8.7 (2025 JCR)** | full OA | €2,000; **€1,000 EASL member** | **Brief Report 2,000 w** | best prestige-per-euro; MEDLINE-indexed |
| **3** | **Hepatology Communications** (AASLD) | 4.6 (JCR year unverified) | full OA | US$2,700 max | unverified | AASLD brand, positioned for work that misses *Hepatology* |
| **4** | **Clinical and Molecular Hepatology** (KASL) | **21.7**, #2 in field | full OA | US$1,500 → **US$2,000 from Feb 2026** | **Research Letter 1,500 w** | home society; full Article = likely desk reject |
| **5** | **Scientific Reports** | **4.9 (2025)** | full OA | US$2,850 | none | strongest formal pre-registration stance |
| **6** | **PLOS ONE** | — | full OA | US$2,477 | none | **explicitly accepts null and replication work** |
| **7** | **F1000Research** | — | OA + open review | **US$1,758** (Brief Report $1,268) | yes | explicitly welcomes negative/confirmatory results |
| **8** | **Liver International** (Wiley) | unverified | hybrid | **free via subscription route** | Research Letter (unverified) | plausible, everything unverified |

Also viable, lower priority: BMJ Open Gastroenterology (IF 2.6, 41-day median
decision, has a short-report type); BMC Medical Genomics (IF 2.6, **3-day**
median first decision, explicitly does not judge on perceived impact);
BMC Gastroenterology (IF 3.2, but US$3,390 — poor value); PeerJ (US$1,195 flat,
sound-science criteria, thin hepatology readership).

**MDPI (Cancers IF 4.8, IJMS IF 5.6, CHF 2,900 each):** fast and near-certain
acceptance, but two real costs — the APC, and that a paper concluding *"this is
not a prognostic biomarker"* sits badly in *Cancers*, where reviewers may push
to reframe H4 as "prognostic potential". That reframing would compromise the
paper. Only if speed dominates.

**Do not pursue:** *Journal of Hepatology* (IF 40.1) and *Hepatology* — certain
desk reject. *Liver Cancer* (Karger) — narrow HCC scope, the hypothesis its
readers care about (H3, field effect) is the one that failed, and the APC is the
highest on the list (up to US$3,920).

---

## 3. On pre-registration — an important negative finding

**No hepatology journal offers Registered Reports.** I could not find a
pre-registration policy at any of them. The Registered Report route at
*Scientific Reports*, *PLOS ONE* and *PeerJ* requires Stage 1 review **before
outcomes are known**, so it is **closed to this manuscript** — the analyses are
done.

What remains is to convert the protocol into visible credibility:
- a **Methods subsection** stating the lock date, the hash, and that every
  amendment (§6a–§6f) is timestamped before the corresponding computation;
- the **full PREREGISTRATION.md deposited as a supplementary file** — including
  §6f, where I misquoted H2's direction and corrected it against the binding §1
  table before reading the result. Leaving that visible is a feature: it
  demonstrates the log is real rather than reconstructed;
- one **cover-letter paragraph** leading with it.

At rigour-based journals this converts the nulls from a liability into the point
of the paper.

---

## 4. Recommended sequence

1. **Do the differentiation-adjusted H2 first.** Everything below is stronger
   after it, and if the effect survives, Tier 1 becomes genuinely reachable.
2. **JHEP Reports, Brief Report** (join EASL first — membership almost certainly
   costs less than the €1,000 APC differential). Lead on H1 + H2; report H3/H4/H6
   as pre-specified secondary outcomes.
3. **Hepatology International** (free via subscription route, 7-day decision) or
   **Hepatology Communications**. The 4,000-word-including-references cap at
   Hep Int forces compression that suits a mixed-result paper.
4. **Scientific Reports** or **PLOS ONE** — where the pre-registration argument
   does the most work.
5. **F1000Research** as the guaranteed-publication backstop.

**Parallel option worth considering:** a **CMH Research Letter (1,500 words)** on
H1 + H2 only. At IF 21.7 in the author's own society journal this is high value
if it lands, and a Research Letter carries a far lower desk-reject risk than a
full Article. It does mean splitting the story.

---

## 5. Before submitting — confirm these

- Whether the institution has a **read-and-publish agreement** (Springer Nature,
  Wiley, Elsevier) or MDPI IOAP membership. A Korean consortium deal could zero
  out several of these APCs and change the ranking outright.
- *Hepatology Communications*: article types, decision times, and the JCR year
  behind IF 4.6.
- *Liver International*: everything — Wiley was unreachable this session.
- *BMJ Open Gastroenterology* and *BMC Medical Genomics*: APCs.
- **Scope boundary with the other two manuscripts** (the R-based CYB5R3 paper
  and the SMG1/LIHC paper). PREREGISTRATION §6b already records the rule; it must
  be restated in the cover letter, and TCGA-LIHC must stay out of this paper —
  it is not outcome-blind.

---

## Sources

JHEP Reports / EASL impact factors: https://easl.eu/news/easl-journals-ifs-2025/ ·
https://www.jhep-reports.eu/content/authorinfo ·
CMH: https://www.e-cmh.org/ and https://e-cmh.org/authors/authors.php ·
AASLD: https://www.aasld.org/journals ·
Hepatology Communications APC (DOAJ): https://doaj.org/api/search/journals/issn%3A2471-254X ·
Hepatology International: https://link.springer.com/journal/12072 and
https://link.springer.com/journal/12072/submission-guidelines ·
Scientific Reports: https://www.nature.com/srep/journal-impact ·
https://www.nature.com/srep/journal-policies/registered-reports ·
PLOS ONE: https://journals.plos.org/plosone/s/criteria-for-publication ·
https://plos.org/publish/fees/ ·
F1000Research: https://f1000research.com/about ·
PeerJ: https://peerj.com/about/aims-and-scope/ ·
MDPI: https://www.mdpi.com/journal/cancers/apc · https://www.mdpi.com/journal/ijms/apc ·
Liver Cancer (DOAJ): https://doaj.org/api/search/journals/issn%3A1664-5553 ·
BMC: https://link.springer.com/journal/12920 · https://link.springer.com/journal/12876 ·
BMJ Open Gastroenterology: https://bmjopengastro.bmj.com/pages/authors/
