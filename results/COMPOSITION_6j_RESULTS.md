# §6j — cell composition and tumor purity

Plan fixed in PREREGISTRATION §6j **before any composition gene was extracted
from any cohort**: the covariate by gene symbol, the tests in order, the joint
model named as decisive, and the decision rule with its three branches.

**The question.** §6g/§6h asked whether the tumor RSI rise is loss of
*hepatocyte identity*, and answered it in two cohorts. This asks the different
question a reviewer asks next: is it a change in **which cells are in the
sample**? D1 cannot settle that — D1 is built from hepatocyte-specific genes, so
a sample with fewer hepatocytes scores lower on D1 for a reason unrelated to the
state of the hepatocytes that remain.

---

## Summary

| | GSE76427 | GSE14520 | GSE164760 |
|---|---|---|---|
| design | 52 pairs | **213 pairs** | unpaired |
| C1 genes present | 21/21 | 20/21 (CLEC4G absent) | 21/21 |
| C1 direction | **lower in tumor**, P = 8.9 × 10⁻⁵ | **lower in tumor**, P = 4.1 × 10⁻⁴ | **lower in tumor**, P = 0.005 |
| unadjusted ΔRSI (median) | +1.129 | +1.261 | *no effect* (P = 0.5) |
| adjusted for **C1 only** | +1.136 (**101%**) | +1.462 (**116%**) | — |
| **joint, D1 + C1** | **+0.451 (40%)**, P = 0.014 | **+0.814 (65%)**, P = 1.8 × 10⁻²² | uninterpretable |
| REDUCTION, joint | +0.373 (**51.20%**) | +0.764 (**89%**) | — |
| DRAIN, joint | −0.364 (**51.12%**) | −0.248 (**31%**) | — |
| corr(ΔD1, ΔC1) · VIF | −0.332 · 1.12 | −0.351 · 1.14 | −0.225 · 1.05 |
| **§6j composite verdict** | **PARTLY explained** | **NOT explained** | uninterpretable |
| **§6j module verdict** | SURVIVES — **but knife-edge** | **SURVIVES** | — |

---

## What the data actually say

**1. The composition difference is real, and its direction was not predicted.**
C1 is **lower in tumor** in all three cohorts. §6j explicitly declined to
predict a direction, and this is why: the comparator is cirrhotic or inflamed
non-tumoral liver, which carries more immune and stromal content than the tumor
does. Predicting that after the fact would have been the manoeuvre the protocol
exists to prevent.

**2. Composition on its own explains none of the rise.** Adjusting for C1 alone
*raises* the estimate slightly in both paired cohorts (101% and 116% retained).
Taken by itself, the objection is empirically wrong: the RSI rise is not
infiltrating cells.

**3. But identity and composition together bite in the small cohort.** The joint
model retains **65% in GSE14520 (P = 1.8 × 10⁻²²)** and only **40% in GSE76427
(P = 0.014)**. The registered rule therefore returns *NOT explained* in the
213-pair cohort and *PARTLY explained* in the 52-pair one. Both are reported.
This is not collinearity — ΔD1 and ΔC1 are close to independent in every cohort
(VIF 1.05–1.14), which was checked precisely so that this could not be the
explanation.

**4. The module split holds decisively in GSE14520 and not at all clearly in
GSE76427.** In the 213-pair cohort REDUCTION retains **89%** against DRAIN's
**31%** — the split §6h reported, essentially intact after removing both
confounds. In the 52-pair cohort the two retentions are **51.20% and 51.12%**.

> ⚠ **The GSE76427 module verdict is a coin flip and must not be leaned on.**
> The §6j rule is a strict inequality between two retentions, so it returns
> "SURVIVES" on a margin of **0.0008** — four orders of magnitude inside the
> estimates' own confidence intervals. The script flags this automatically. The
> honest statement is that in GSE76427, once *both* confounds are removed, the
> reduction and drain arms are **no longer distinguishable from one another**.

**5. C2 changes nothing.** The deliberately over-adjusted covariate — C1 plus
KRT19, KRT7, SOX9 and EPCAM, which carry tumor biology as well as composition —
moves every estimate by less than 0.03. Whatever is driving these results, it is
not the cholangiocyte/progenitor compartment.

**6. GSE164760 is uninterpretable here, as it was in §6h.** Its unadjusted
tumor effect is +0.120, P = 0.5 — there is no effect for an adjustment to
explain, so every retention ratio is a ratio to approximately zero. §6j said in
advance that this cohort is supporting context and not a replication. Its
mechanical verdict is left standing in the JSON and is **not** counted.

---

## Where this leaves the paper

**The composition objection does not overturn the finding — but it does narrow
the claim, and the narrowing is real.**

- The headline result now rests **squarely on GSE14520**: 213 patient-matched
  pairs, HBV etiology, a third platform, and a reduction arm that retains 89% of
  its shift after adjusting for hepatocyte identity *and* cell composition
  simultaneously.
- **GSE76427 no longer supports the module split once both confounds are
  removed.** It still supports the composite rise (+0.451, P = 0.014), but the
  dissociation between the two arms disappears there. With 52 pairs this may
  simply be power — the confidence intervals are wide enough to contain
  GSE14520's split — but that is an explanation, not a result, and it is
  labeled as such.
- **The reframing proposed in improvement B is still right, but its evidence
  base is now one cohort, not two.** The manuscript must say "the reduction arm
  survives adjustment for identity and composition in the larger of two paired
  cohorts, and the two arms are not separable in the smaller" — not "the split
  replicates".

### What this does not license

- None of this is causal. Adjusting for a covariate is not an experiment.
- **AIFM2's absence from HG-U133A remains a limitation, but it runs the other
  way** (corrected in PREREGISTRATION §6k, 2026-08-27). GSE14520's REDUCTION is
  a 3-gene score, and AIFM2 is the *strongest* member of REDUCTION where it is
  measured (GSE76427 Δz = +1.292, P = 2.31 × 10⁻¹⁶; GSE164760 Δz = +0.531,
  P = 0.0437). The platform removes the module's most informative gene from the
  cohort now carrying the result, so the 89% REDUCTION retention is achieved
  without it. The limitation has grown in importance, not shrunk — but it makes
  this cohort's retention conservative, not generous.
- CLEC4G is absent from HG-U133A, so GSE14520's C1 is 20 genes rather than 21.
  Reported, not substituted.

### Next

1. **Improvement A3 is now the highest-value remaining analysis, not A2.** A
   third paired cohort that carries AIFM2 would simultaneously measure the
   module's strongest gene under a paired design and tell us whether GSE76427's
   null split is power or substance.
2. A2 (CPTAC protein-level validation) remains valuable but no longer answers
   the sharpest open question.
