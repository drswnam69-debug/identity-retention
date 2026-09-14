"""
rsi_config.py -- LOCKED definition of the hepatic Reductive-Supply Index (RSI).

Project : SMG1C - HMGCR - CYB5R3 - MTARC1 axis in MASLD and HCC
Analysis: Direction 2, human liver transcriptome original article
Locked  : 2026-08-25

DO NOT EDIT the module membership or the index formula after the lock date.
Any change invalidates the pre-specification. If a change is unavoidable,
create a new file (rsi_config_v2.py), record the reason, and report BOTH
versions in the manuscript.

Module rationale
----------------
S  (supply)     isoprenoid / CoQ biosynthetic capacity: how much CoQ10 the
                mevalonate pathway can produce.
R  (reduction)  CoQ-reducing capacity: whether that CoQ10 is kept in the
                reduced, antioxidant form (CoQH2).
D  (drain)      pathogenic consumption of the same electron pool through the
                N-reductive system (NADH -> CYB5R3 -> CYB5B -> mARC1).
SHARED          adapters used by both the protective and the pathogenic arm;
                direction cannot be assigned, so they are reported separately
                and are NOT part of the index.
E  (effector)   downstream ferroptosis effectors, reported separately.
U  (upstream)   SMG1 complex and NMD machinery, reported separately.
"""

import hashlib
import json

LOCK_DATE = "2026-08-25"

# --- index modules (these three define RSI) ---------------------------------

MODULE_SUPPLY = [
    "HMGCR", "HMGCS1", "MVK", "PMVK", "MVD", "IDI1", "FDPS",
    "PDSS1", "PDSS2", "COQ2", "COQ3", "COQ5", "COQ6", "COQ7", "COQ9",
]

MODULE_REDUCTION = [
    "CYB5R3", "CYB5R1", "AIFM2", "NQO1",
]

MODULE_DRAIN = [
    "MTARC1", "MTARC2", "POR",
]

# --- reported separately, NOT in the index ----------------------------------

MODULE_SHARED = ["CYB5A", "CYB5B"]

MODULE_EFFECTOR = ["GPX4", "SLC7A11", "SCD1", "ACSL4", "LPCAT3"]

MODULE_UPSTREAM = ["SMG1", "SMG8", "SMG9", "UPF1", "UPF2", "SREBF2"]

# --- sensitivity analysis S2: cholesterol-arm genes added to SUPPLY ---------
# If RSI behaves identically once these are included, the claim that the
# index is CoQ-specific rather than cholesterol-driven is WEAKENED.
SENSITIVITY_S2_CHOLESTEROL = ["SC5D", "DHCR7", "DHCR24", "MSMO1"]

# --- genes highlighted individually in the main figures ---------------------
FOCUS_GENES = ["HMGCR", "MVD", "CYB5R3", "MTARC1", "AIFM2", "GPX4", "SCD1"]

# --- pre-specified directional hypotheses -----------------------------------
HYPOTHESES = {
    "H1": "RSI increases monotonically across the MASLD spectrum "
          "(control < NAFL < MASH F0-1 < F2 < F3 < F4).",
    "H2": "RSI is higher in HCC tumor than in adjacent non-tumor liver.",
    "H3": "RSI in adjacent non-tumor NASH liver is already higher than in "
          "NASH liver from patients without HCC (field effect).",
    "H4": "Higher RSI predicts worse overall survival in HCC (HR > 1).",
    "H5": "NMD-target burden (a proxy for SMG1 activity) is inversely "
          "correlated with RSI (rho < 0).",
    "H6": "RSI is positively correlated with a ferroptosis-resistance "
          "signature (rho > 0).",
}

# --- index formula ----------------------------------------------------------
# RSI = 0.5 * (mean_z(SUPPLY) + mean_z(REDUCTION)) - mean_z(DRAIN)
#
# SUPPLY and REDUCTION are weighted 0.5 each because SUPPLY has 15 genes and
# REDUCTION has 4; a pooled mean would let SUPPLY dominate the index, which
# does not match the biological claim ("how much is made" x "how much is
# reduced").
WEIGHT_SUPPLY = 0.5
WEIGHT_REDUCTION = 0.5
WEIGHT_DRAIN = 1.0

ALL_MODULES = {
    "SUPPLY": MODULE_SUPPLY,
    "REDUCTION": MODULE_REDUCTION,
    "DRAIN": MODULE_DRAIN,
    "SHARED": MODULE_SHARED,
    "EFFECTOR": MODULE_EFFECTOR,
    "UPSTREAM": MODULE_UPSTREAM,
}


def lock_hash() -> str:
    """Stable hash of the locked definition. Record this in the manuscript."""
    payload = json.dumps(
        {
            "lock_date": LOCK_DATE,
            "supply": MODULE_SUPPLY,
            "reduction": MODULE_REDUCTION,
            "drain": MODULE_DRAIN,
            "weights": [WEIGHT_SUPPLY, WEIGHT_REDUCTION, WEIGHT_DRAIN],
            "hypotheses": HYPOTHESES,
        },
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode()).hexdigest()


if __name__ == "__main__":
    for name, genes in ALL_MODULES.items():
        print(f"{name:<10} n={len(genes):>2}  {', '.join(genes)}")
    print(f"\nS2 sensitivity add-on: {', '.join(SENSITIVITY_S2_CHOLESTEROL)}")
    print(f"\nLOCK DATE : {LOCK_DATE}")
    print(f"LOCK HASH : {lock_hash()}")
