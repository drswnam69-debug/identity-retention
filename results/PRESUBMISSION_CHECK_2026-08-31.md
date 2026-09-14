# 투고 전 최종 검증 (2026-08-31)

인용문헌 **12–21번은 이번에 처음 독립 검증**했다(지난번은 1–11번만 했다).
동시에 편집 이력이 남긴 손상을 찾기 위해 원고 전체를 새로 통독시켰다.

---

## A. 초록이 철회한 생화학을 그대로 싣고 있었다 — 최우선 수정

**수정 전 초록 첫 문장:**
> Coenzyme Q is both the substrate of a protective reductive system and **the electron
> source of a pathogenic N-reductive system** ...

CoQ는 mARC의 전자 공급원이 아니다. 서론·Methods·§6l·커버레터가 모두
NADH → CYB5R3 → cytochrome b5 → mARC로 정정해 놓았는데, **정작 초록만 옛 전제를
그대로 주장**하고 있었다. 많은 독자가 초록만 읽는다는 점에서 가장 위험한 잔존 오류였다.

**수정 후:**
> Hepatic NADH supplies both a protective system that keeps the coenzyme Q pool reduced
> and a pathogenic mitochondrial amidoxime-reducing component (mARC) system that competes
> for the same electrons.

초록에 **POR 오류 자체도 한 문장 추가**했다 — 본문 Limitations가 "more consequential"이라
부르고 커버레터가 헤드라인으로 삼는 사안인데 초록에만 없었다.
> Recomputing without POR, locked into the drain on a mistaken premise, changed no verdict.

단어 수를 맞추기 위해 6곳을 다듬어 **249단어**(한계 250) 유지.

---

## B. 내가 지난번에 과잉교정한 것을 되돌렸다

인용 [20] Yan et al.을 "POR and, **to a lesser degree**, CYB5R1"로 고쳤었다.
논문 초록을 직접 확인하니 **두 효소를 대등하게** 다룬다("Genetic knockout of POR **and**
CYB5R1 decreases cellular hydrogen peroxide generation"; "POR/CYB5R1 oxidoreductase").
"CYB5R1의 활성이 낮다"는 서술은 논문이 아니라 **논평(H1 Connect)**에서 온 것이었다.
→ 순위 표현 삭제, "POR and CYB5R1"로 환원. 근거 없는 정량적 서열은 심사자가 바로 짚는다.

---

## C. 인용문헌 12–21 검증 결과

**서지 오류 0건.** 10편 전부 저자·제목·저널·연도·권·페이지가 실제 기록과 일치.
날조·잘못된 저널·잘못된 연도 없음. [13]의 `2025;83:1338–1352`도 정확.

**GEO 6건 전부 원 출판물·플랫폼·표본수가 GEO 기록과 일치:**

| 등록번호 | 플랫폼 (GEO) | n | 연결 PMID |
|---|---|---|---|
| GSE135251 | GPL18573 Illumina NextSeq 500 | 216 | 33268509 |
| GSE130970 | GPL16791 Illumina HiSeq 2500 | 78 | 31467298 |
| GSE167523 | GPL21290 Illumina HiSeq 3000 | 98 | 34105780 |
| GSE76427 | GPL10558 Illumina HumanHT-12 V4 | 167 (115 T + 52 NT) | 29117471 |
| GSE164760 | GPL13667 Affymetrix HG-U219 | 170 (53/29/74/8/6) | 33992698 |
| GSE14520 | GPL571 + **GPL3921** (HT_HG-U133A) | 488 총계, GPL3921 445 | 21159642 |

- **GSE167523 → Kozumi et al.** 연결 확인(GEO 기여자 Kodama·Murai·Takehara가 저자와 일치).
- **GSE14520 2개 플랫폼·488개·GPL3921 445개** 전부 확인. 445 = 225 T + 220 NT,
  488 − 445 = 43(GPL571).
- **213쌍은 GEO에 없는 유도값**이라는 지적 → 이미 환자 ID에서 직접 재유도해 두 번 확인했고
  (52쌍, 213쌍), 본문 Methods에 유도 방법을 명시해 두었으므로 그대로 유지.
- [13] 저자 이니셜 Kwan KKL → **Kwan KK**(PubMed 색인형)로 수정.

---

## D. 편집 이력이 남긴 손상 — 실제 결함 14건

| 항목 | 문제 | 조치 |
|---|---|---|
| Code availability | "eleven amendments" — §6l 추가 전 값 | twelve |
| 커버레터 | "**Two** concern the analysis log; the fourth…" — 산술이 깨짐(§6c·§6i·§6k = 3개) | Three |
| Methods | "both denominators are given in Table 2" — 표에는 중앙값만 있다 | "the composite denominator" |
| Limitations | "**A second** and more consequential error" — 첫 번째 오류가 없다(앞 두 항목은 플랫폼 한계) | "A more consequential error, recorded as §6l," |
| Limitations | AIFM2 문장이 자기모순 — "conservative rather than generous"가 이미 방향인데 "not comparable **in either direction**" | §6k가 실제로 기록한 논리로 복원 |
| Results 소제목 | "changes **nothing**" 과장 — 같은 문단이 51.1%→61.9% 이동을 보고 | "changes **no verdict**" |
| Discussion | "**most** of that fall tracks identity loss" — GSE76427(51.1%)과 DRAIN′(61.9%)이 반증 | "in the larger cohort, most of" |
| 서론 [13] | "**On the protective side** … 메발론산 경로가 간암을 **촉진**한다 [13]" — 인용이 정반대 편에 붙어 있었다 | 문장 분리 후 "That same capacity is not benign in tumors:" |
| 서론 [4–7] | [7]은 "종간 차이"로 가설을 **한정**하는 논문인데 "hepatoprotective [4–7]"로 뭉뚱그림 | "[4–6], though not identically across species [7]" |
| §6l Results | "both slightly more strongly … in GSE14520" 수식 오류 — GSE76427은 RSI′ +0.981 < +1.129로 **작다** | 코호트 한정 재작성 |
| Methods | "**pre-specified** sensitivity analysis" — D2/C2(원 잠금)와 DRAIN′(§6l 사후)을 같은 말로 부름 | "specified in §6l before it was computed" |
| Results | GSE164760 배제 규칙이 **세 번** 반복 | 1회로 축약 |
| Methods | MASH가 정의 전에 사용 | 첫 사용 시 전개 |
| Table 2 | `*… without* POR` — 이탤릭이 유전자명 앞에서 닫힘 | 범위 수정 |
| Title Page | Data availability가 본문 Code availability와 모순(코드가 이미 제출됐다고 주장) | 본문에 맞춤 |
| 전체 | 통계 기호 이탤릭 불일치 — Results에만 `P =` 14곳이 평문(편집 이음매 흔적) | `*P*`·`*z*`·`*n*`·`*I*²` 전면 통일 |

---

## E. 그림·테이블 재검증

- Figure 1: 유전자 15/4/3, §6l POR 주석, forest ρ = 0.372·I² = 0%·n = 216/78/98/392 — 전부 일치
- Figure 2: 중앙값 +1.13/+1.26, 양성 짝 88%/93%, Cliff δ +0.76/+0.86 — 원자료 재계산 일치
- Figure 3: 3패널 구성이 범례와 일치, 각주 VIF 1.12·1.14
- 그래픽초록: DRAIN 상자·캐스케이드 라벨 수정본 반영
- **Table 1 표본수를 phenotype 파일에서 재계산**: 167 / 445 / 170 (53·29·74·8·6) — GEO와도 일치
- **짝 수 재유도**: 52, 213
- **Table 2 41개 수치 전부 원본 JSON과 일치**

---

## F. 최종 수치

- 초록 **249단어**(한계 250) · 본문+참고문헌 **3,878단어**(한계 4,000)
- 인용 1→21 완전 오름차순, 21편 전부 인용, 누락·유령 0
- 표 2 · 그림 3 · 그래픽초록 1 · 보충파일 1 (한계 6 이내)
- 조항 12개(§6a–§6l), 그중 4개가 저자 오류를 기록
- 마크다운 마커 불균형 0, 표 열 수 정상, we/our 0건, 영국식 철자 0건(참고문헌 제외)
- .docx 13쪽 렌더링 육안 확인
