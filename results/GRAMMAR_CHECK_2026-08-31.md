# 문법·철자 교정 (Grammarly + LanguageTool, 2026-08-31)

## 도구

Grammarly MCP 커넥터는 존재하지 않는다(레지스트리 검색 결과 없음). 대신 두 가지를 썼다.

1. **Grammarly Pro** — 사용자의 로그인된 계정(drswnam69@gmail.com)에 Chrome으로 접속해
   원고 본문과 커버레터를 붙여넣고 제안을 하나씩 검토. 유전자 기호·통계 기호·인용부호는
   Grammarly가 오탐을 쏟아내므로, 붙여넣기 전에 그런 요소를 일반 단어로 치환한 산문
   버전을 만들어 넣었다.
2. **LanguageTool 6.8** — 컨테이너에 로컬 설치해 마크다운 원본에 직접 실행. 유전자 기호·
   데이터셋 번호·저자명 등 약 300개 용어 허용목록을 만들어 오탐을 걸렀다.

두 도구를 함께 쓴 이유: Grammarly는 문장 구조를 잘 잡고, LanguageTool은 마크업이 붙은
원본 파일에 직접 돌릴 수 있어 서식 문제(이탤릭 누락 등)를 잡는다.

## 점수

| | 처음 | 교정 후 |
|---|---|---|
| 본문 Writing quality | 85 | **89** |
| 본문 Grammar 지적 | 24 | **2** (둘 다 오탐) |
| 커버레터 Writing quality | — | **92** |
| 커버레터 Grammar 지적 | 2 | **0** |

## 실제로 고친 것 (24건)

**주어-동사 불일치 1건 (커버레터)**
- "each of the amendments that **governs** a computation" → 선행사가 복수이므로 틀렸다.
  단수형으로 재구성: "each amendment that governs a computation".

**동사 형태 1건**
- "the effect **reproduced** at nearly identical magnitude" → "the effect **was** reproduced".
  reproduce를 자동사로 쓰면 "번식하다"가 되어 심사자가 짚을 수 있다.

**도입구 뒤 쉼표 9건**
Under the over-adjusted D2 model_, / At the module level_, / In GSE14520_, / In GSE76427_, /
Under the joint model_, / In the larger cohort_, / In the smaller cohort_, /
Where AIFM2 is measured_, / in the same model_, / By the threshold fixed in the protocol_,

**등위접속사 앞 쉼표 5건** (독립절 두 개를 잇는 and 앞)
different verdicts_, and I report / clears 50%_, and its rise / in GSE14520)_, and DRAIN′ /
expression level)_, and the pooled correlation / No new data were generated_, and there is

**접속부사 쉼표 1건**
- "instead the error is disclosed" → "instead_, the error is disclosed"

**연속 쉼표(Oxford comma) 4건**
- 원고가 어떤 목록에는 쓰고("CYP2E1, CYP3A4, CYP1A2, and CYP2C9") 어떤 목록에는
  안 쓰고 있었다("MTARC1, MTARC2 and POR"). 미국식 표준에 맞춰 통일.

**서식 2건**
- `**P = 8.1 × 10⁻⁸**` — 굵게 표시 안에 있어서 지난번 이탤릭 통일 작업이 놓친 P 두 개.
  다른 곳과 같은 `*P* =` 형태로 정리.

**표현 1건**
- "a **rationale error** that does not drive any result" → Grammarly는 "rational error"로
  고치라고 했는데 그건 뜻이 완전히 달라진다(합리적인 오류). 의도는 "근거의 오류"이므로
  **"an error in the rationale, not one that drives any result"**로 재작성.

**세 문장 연속 "In ..." 시작 1건**
- 세 번째(대조되는 코호트) 문장을 "H2 was not supported in GSE164760..."로 바꿔 단조로움 해소.

**철자 1건**
- "partialing/partialling out"은 미국식 표기가 갈리는 단어라 아예 피했다 →
  "survived **adjustment for** each sample's mean expression".

## 거부한 제안

- **"in tumor" → "in the tumor" / "in tumors"** (3건). 병리·유전체 문헌에서 "in tumor"는
  조직 구획을 가리키는 표준 용법이고 원고 전체가 일관되게 쓴다. 바꾸면 오히려 어색해진다.
- **"a SHA-256" → "an SHA-256"**. SHA는 "샤"로 읽으므로 a가 맞다.
- **"the covariate"** (3건). LanguageTool이 covariate를 명사로 인정하지 않아 생긴 오탐.
- **"log₂-transformed"**. 아래첨자 때문에 파서가 깨진 것.
- **"analyses" → "analyzes"**. 명사 복수형 analyses와 동사 analyze는 둘 다 미국식으로 맞다.
- **제목 끝에 마침표 추가**. 제목에는 마침표를 찍지 않는다.
- **Clarity 47건 / Delivery 1건**. 대부분 "문장을 나눠라", "더 간결하게" 같은 일반 산문용
  제안이라 학술 원고에는 맞지 않아 적용하지 않았다.

## 최종 수치

- 초록 **250단어**(한계 250) · 본문+참고문헌 **3,882단어**(한계 4,000)
- 인용 1→21 완전 오름차순, 21편 전부 인용, 누락·유령 0
- 마크다운 마커 불균형 0, 표 열 수 정상
- .docx 13쪽 렌더링 육안 확인
