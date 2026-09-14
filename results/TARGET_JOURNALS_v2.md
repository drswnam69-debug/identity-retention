# 투고 전략 v2 — RSI 원저 (방향 2)

개정 2026-08-26, PREREGISTRATION §6h/§6i 완료 후.
아티팩트: 어디에 보낼 것인가 — https://claude.ai/code/artifact/fa6d1c54-680f-4a30-a938-cbdcd6fdda63

**출처 규율.** 모든 수치는 저널·발행사·학회 자체 페이지 또는 DOAJ에서 확인했다.
집계 사이트(manusights, journalmetrics, askbisht, resurchify, scijournal, bioxbio 등)는
사용하지 않았다. 확인하지 못한 항목은 "미확인"으로 표시했고 추정하지 않았다.
Wiley(onlinelibrary.wiley.com), Elsevier(sciencedirect.com), Ovid(edmgr.ovid.com)는
두 번의 조사 모두 네트워크 정책에서 차단되었다.

---

## 1. 무엇이 달라졌는가

전 판본의 최우선 권고 — "Tier 3 위로 가려면 탈분화 교란부터 해결하라" — 가 완료되었다.

- H2가 두 코호트에서 지지 (Cliff's δ +0.76, +0.86; 52쌍, 213쌍; 병인과 플랫폼이 다름)
- 탈분화 보정이 두 코호트에서 재현: REDUCTION 71%/84% 대 DRAIN 41%/25%
- 과보정 D2에서 DRAIN은 두 코호트 모두 유의성 상실(P=0.22, 0.36), REDUCTION은 69%/81% 유지

**전략적 재정의.** 이 논문은 더 이상 "새로운 간 지표"가 아니다. 지금 있는 것은
*MASLD→HCC 진행에서 간의 전자 배분이 환원 예비 쪽으로 이동하고, 그 이동의 환원 팔은
탈분화로 설명되지 않는다*는 기전 주장이다. 제목·초록을 이 방향으로 다시 쓰는 것이
무료로 얻는 가장 큰 개선이다.

## 2. 새 부채 순서

1. 새 데이터 0, 습식 실험 0, 직교 검증 0 (여전히 1번, 그러나 A2로 절반 해소 가능)
2. **세포 조성·종양 순도 — 미해결.** D1 보정은 간세포 *정체성*을 보정한 것이지
   검체 내 간세포 *비율*을 보정한 것이 아니다. 탈분화와 정확히 같은 종류의 반론이다.
3. AIFM2 부재가 가장 날카로운 결과 쪽에 있다 (구조적으로 가설에 유리)
4. 6개 가설 중 3개 미지지
5. "Index"라는 명명이 바이오마커 기대를 부르고, 그 답이 H4 영가설

## 3. 순위

| | 저널 | IF (연도) | APC | 단형식 | 판단 |
|---|---|---|---|---|---|
| 1 | **JHEP Reports** | 8.7 (2025) | €2,000 · EASL 회원 €1,000 | Brief Report 2,000w | 1차 추천. Original 6,000w로 A1·A2 완료 후 |
| 2 | **Hepatology International** | 7.8 (2025) | **0** (구독 경로 확인됨) | Letter 500w | 비용 0. Original 4,000w(참고문헌 포함) |
| 3 | **Clin Mol Hepatol** | 21.7 | US$2,000 | Research Letter 1,500w · 그림 또는 표 1개 · 참고문헌 15개 | 고위험 고수익. 본편과 동시 진행 불가 |
| 4 | **Antioxidants (MDPI)** | 8.2 (2025) | CHF 2,900 | Brief Report | **유일하게 사전등록 조항 있음.** 1차 결정 18.7일 |
| 5 | Hepatology Communications | 4.6 대 6 (불일치, 미확인) | 최대 US$2,700 | 미확인 | 브라우저로 직접 확인 필요 |
| 6 | Liver International | 미확인 | 미확인 | 미확인 | Wiley 차단, 전부 미확인 |

**제외:** J Hepatol·Hepatology(확실한 desk reject) · Redox Biology·FRBM(최근 게재
논문에 순수 계산 연구 전무) · Cell Reports Medicine(기전적 진전 부족 위험) ·
Liver Cancer(범위 좁음, 실패한 H3가 그 독자의 관심사) · Cancers(H4를 "예후 가능성"으로
재구성하라는 압력 → 정직성 훼손).

**사전등록 음성 소견 유지:** 간학회 저널 중 Registered Report 운영 없음. Sci Rep·PLOS
ONE·PeerJ의 RR 트랙은 Stage 1 사전 심사를 요구하므로 이 원고에는 닫혀 있다.

## 4. 개선 — 노력 대비 효과 순

### A1 · 세포 조성·종양 순도 보정 (공개 데이터, 2–3일, 효과 큼)
다음 리뷰어 2의 질문. 두 경로: (1) 디컨볼루션(CIBERSORTx·xCell·EPIC)으로 세포 분율을
§6g와 같은 방식으로 공변량에 투입, (2) 공개 HCC 단일세포 아틀라스(예: GSE149614)에서
CYB5R3·AIFM2·NQO1과 MTARC1·POR이 간세포/악성세포에서 우세함을 제시. (2)가 더 설득력
있고 더 싸다. **§6j로 먼저 등록하고 계산할 것.**

### A2 · CPTAC 간암 프로테오믹스 직교 검증 (공개 데이터, 3–5일, 효과 가장 큼)
"습식 실험 없다"는 1번 부채를 피펫 없이 절반 해소하는 유일한 수단. 종양·인접 짝지은
단백질 정량이 공개되어 있다. CYB5R3 단백질 상승 + mARC1 단백질 하강이 확인되면 다른
분자층·다른 코호트·다른 측정 기술에서의 재현이 된다.
**주의:** ① §6b의 원고 간 범위 규칙에 CPTAC이 걸리는지 먼저 판단 (TCGA-LIHC 배제 근거가
여기에 적용되는지는 별개). ② 결과를 보기 전에 방향과 판정 규칙을 등록.

### A3 · AIFM2 포함 세 번째 짝지은 코호트 (공개 데이터, 2–4일, 방어력 상승 큼)
GSE25097·GSE36376·GSE64041 등 (플랫폼 커버리지 확인 필요). 4유전자 REDUCTION으로
재현하면 "구조적으로 가설에 유리" 한계가 사라진다.

### B · 프레이밍 (비용 0, 반나절)
- 제목·초록을 환원 팔로 재작성. 복합 지표는 사전등록된 도구로 보고하되 논문의 이름이
  되지 않게 한다.
- 음성 결과를 초록에 설계로 올린다 ("6개 중 3개 미지지, 명세대로 보고").
- Methods에 사전등록 소절: 잠금일, SHA-256, §6a–§6i가 각 계산 이전에 시간 기록.
- `PREREGISTRATION.md` 전문을 보충자료로. **§6f(H2 방향 오기)와 §6i(HNF1A 누락)를
  지우지 말 것** — 그 두 정정이 로그가 사후 재구성이 아니라는 증거다.
- 코드·파생 데이터를 Zenodo에 올려 DOI 확보.

### C · 더 높은 곳을 노릴 때만
- C1: H5 전사체 수준 정량 (recount3 `SRP217231`)
- C2: GSE164760 전장 재준비 (H6 해금) — 우선순위 낮음
- C3: CYB5R3·mARC1 면역조직화학 (IRB, 수개월). *J Hepatol*급을 여는 유일한 길이지만
  **이 원고를 붙잡고 기다리지 말 것** — 다음 논문의 씨앗으로.

## 5. 제출 전 확인

- [ ] 기관 read-and-publish 계약 (Springer Nature·Wiley·Elsevier 한국 컨소시엄) — 순위가
      통째로 바뀔 수 있다
- [ ] EASL 연회비 (€1,000 절감보다 싼가) — easl.eu 가입 페이지가 응답하지 않았다
- [ ] Liver International·Hepatology Communications 기고 형식·단어 수·APC를 브라우저로 직접
- [ ] §6b 범위 규칙과 CPTAC의 관계를 커버레터에 명시. TCGA-LIHC는 계속 제외(결과 비맹검)
- [ ] A1·A2를 §6j·§6k로 등록한 뒤 계산 — 여기서 규율을 깨면 최대 자산이 손상된다

## 출처

easl.eu/news/easl-journals-ifs-2025/ · jhep-reports.eu/content/authorinfo ·
link.springer.com/journal/12072 및 /submission-guidelines 및 /how-to-publish-with-us ·
e-cmh.org/authors/authors.php · mdpi.com/journal/antioxidants 및 /instructions 및 /apc ·
wolterskluwer.com/en/solutions/ovid/hepatology-communications-16964 ·
doaj.org (2471-254X, 2213-2317, 2212-8778, 2054-4774) ·
link.springer.com/journal/12920/how-to-publish-with-us · bmjopengastro.bmj.com/pages/authors ·
journals.sagepub.com/author-instructions/are · shop.elsevier.com (Redox Biology, FRBM, Molecular Metabolism)
