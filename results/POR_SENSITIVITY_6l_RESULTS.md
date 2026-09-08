# §6l — POR sensitivity analysis: results

**Run:** 2026-08-30 · **Plan:** PREREGISTRATION §6l, recorded before this was computed.
**Index lock hash unchanged:** `7a2bf934fe6e51184a573057c93e0cda0f3e79e2c538224efb35cf728b680046`

## 0. Integrity check

The locked index was re-derived from the expression matrices inside this script
and compared against `rsi.tsv` as written by `02_compute_rsi.py`:

| Cohort | Pearson r | max abs deviation |
|---|---|---|
| GSE76427 | 1.0000000000 | 0.00 |
| GSE14520 | 1.0000000000 | 0.00 |

The machinery reproduces the locked index exactly, so any difference below is
attributable to dropping *POR* and to nothing else. D1 is 22/22 and C1 is 21/21
(GSE76427) and 20/21 (GSE14520, *CLEC4G* absent), the same as in §6g/§6j.

## 1. What was computed

**DRAIN′ = mean_z(*MTARC1*, *MTARC2*)** — the N-reductive arm proper, without *POR*
**RSI′ = ½(mean_z SUPPLY + mean_z REDUCTION) − mean_z DRAIN′**

put through exactly the §6g and §6j machinery, unchanged.

## 2. Direction is unchanged in both cohorts

| Cohort | quantity | paired median | Wilcoxon P |
|---|---|---|---|
| GSE76427 | RSI (locked) | +1.129 | 8.1 × 10⁻⁸ |
| GSE76427 | **RSI′** | **+0.981** | **4.5 × 10⁻⁷** |
| GSE76427 | DRAIN (locked) | −0.595 | 1.3 × 10⁻⁷ |
| GSE76427 | **DRAIN′** | **−0.611** | **1.9 × 10⁻⁵** |
| GSE14520 | RSI (locked) | +1.262 | 4.1 × 10⁻³⁴ |
| GSE14520 | **RSI′** | **+1.292** | **2.0 × 10⁻³¹** |
| GSE14520 | DRAIN (locked) | −0.686 | 1.9 × 10⁻²⁴ |
| GSE14520 | **DRAIN′** | **−0.756** | **4.8 × 10⁻¹⁸** |

RSI′ still rises in tumor; DRAIN′ still falls. In GSE14520 both effects are
slightly **larger** without *POR*.

## 3. The module dissociation survives, unchanged

Retained fraction under the joint D1 + C1 model, on the same denominator the
manuscript uses for module retentions (the paired mean of the module shift):

| Cohort | REDUCTION | DRAIN (locked, with *POR*) | **DRAIN′ (without *POR*)** |
|---|---|---|---|
| **GSE14520** | 89.0% | 31.0% | **31.8%** |
| GSE76427 | 51.2% | 51.1% | 61.9% |

- **GSE14520: the dissociation is intact** — 89.0% versus 31.8%, against 89.0%
  versus 31.0% with *POR*. A difference of 0.8 percentage points. *POR* is not
  carrying this result.
- **GSE76427: the two arms remain not distinguishable**, and the knife-edge
  (51.20% vs 51.12%) actually tips the other way without *POR* (51.2% vs 61.9%).
  This reinforces, rather than undermines, the manuscript's refusal to claim a
  split in this cohort.

## 4. Composite retention

| Cohort | RSI (locked) joint | RSI′ joint |
|---|---|---|
| GSE76427 | +0.451, 39.9% | +0.465, 47.3% |
| GSE14520 | +0.814, 64.5% | +0.798, 61.8% |

Both cohorts move by a few percentage points and neither crosses a verdict
boundary except GSE76427, which remains below the pre-registered 50% threshold
under both definitions (39.9% and 47.3%).

## 5. Verdict under the rule fixed in §6l

**Branch 1 applies.** Direction and every qualitative verdict are unchanged:
RSI′ rises, DRAIN′ falls, the GSE14520 module dissociation holds, and GSE76427
remains uninformative for the split. Therefore:

- *POR*'s inclusion in DRAIN is reported as a **rationale error that does not
  drive any result**.
- The locked index stands. It is not redefined.
- DRAIN′ is reported in the manuscript as a pre-specified sensitivity analysis,
  and the error is disclosed in the Introduction, Methods and Limitations.

## 6. A second inconsistency found while doing this

The §6g/§6j pipeline computes **composite** retained fractions against the paired
**median** of ΔRSI, but **module** retained fractions against the paired **mean**
of the module shift. Both conventions are defensible; using two of them without
saying so is not. No verdict depends on the choice — on the median denominator
GSE14520's split reads 94.4% versus 36.2% instead of 89.0% versus 31.0%, and the
GSE76427 arms remain indistinguishable either way — but the manuscript must state
which denominator is used where. This is disclosed rather than harmonised, since
harmonising it now would mean changing numbers computed under the registered plan.
