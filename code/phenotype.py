"""phenotype.py -- map GEO phenotype text onto the ordered MASLD spectrum."""

from __future__ import annotations

import re
import pandas as pd

# Ordered stages used by hypothesis H1 (Jonckheere-Terpstra alternative).
STAGE_ORDER = ["control", "NAFL", "NASH_F0-F1", "NASH_F2", "NASH_F3",
               "NASH_F4"]

# Cohorts annotate severity differently. Each axis below is an ordered ladder
# of increasing severity; which one a cohort can support is a property of its
# GEO metadata, not a choice made after seeing results. See PREREGISTRATION.md,
# "Replication axes", fixed 2026-08-25 before any replication cohort was run.
AXES = {
    "group":    STAGE_ORDER,                        # GSE135251
    "fibrosis": ["F0", "F1", "F2", "F3", "F4"],     # GSE130970
    "subtype":  ["NAFL", "NASH"],                   # GSE167523
    # GSE164760 is a tissue-type series, not a severity gradient within one
    # tissue. The ladder is ordered by proximity to cancer and is used only
    # for display and an exploratory trend; the PHASE C hypotheses are
    # specific pairwise contrasts, not a trend across this ladder.
    "tissue":   ["Healthy liver", "NASH liver", "Cirrhotic liver",
                 "Non-tumoral NASH liver adjacent to HCC", "NASH-HCC tumor"],
    # GSE76427 is a mixed-etiology HCC series with paired adjacent liver. It
    # is NOT the GSE164760 NASH ladder and must not be folded onto it: its
    # adjacent tissue is not NASH-specific. Two levels only, tumor last.
    "tissue_hcc": ["Adjacent non-tumor liver", "HCC tumor"],
}

# Raw GEO tissue strings -> the "tissue_hcc" axis. Kept explicit so an
# unrecognized label fails loudly rather than being silently dropped.
TISSUE_HCC_ALIASES = {
    "adjacent non-tumor liver tissue": "Adjacent non-tumor liver",
    "primary hepatocellular carcinoma tumor": "HCC tumor",
}


def normalize_stage(value: str) -> str | None:
    """Normalize a free-text group label to one of STAGE_ORDER."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    v = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if v in {"control", "normal", "healthy", "healthy_control"}:
        return "control"
    if v in {"nafl", "steatosis", "nafld_nafl"}:
        return "NAFL"
    m = re.search(r"nash.*?f\s*_?(\d)(?:\s*_?f?\s*(\d))?", v)
    if m:
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) else lo
        if lo == 0 or (lo, hi) == (0, 1) or (lo, hi) == (1, 1):
            return "NASH_F0-F1" if hi <= 1 else f"NASH_F{hi}"
        return f"NASH_F{lo}"
    if v.startswith("nash"):
        return "NASH_F0-F1"
    return None


def stage_from_fibrosis(disease: str, fibrosis: str) -> str | None:
    """Fallback: derive the stage from separate disease / fibrosis columns."""
    d = str(disease).strip().lower()
    if d in {"control", "normal", "healthy"}:
        return "control"
    try:
        f = int(float(re.sub(r"[^0-9.]", "", str(fibrosis))))
    except (ValueError, TypeError):
        return None
    if "nafl" in d and "nash" not in d:
        return "NAFL"
    return "NASH_F0-F1" if f <= 1 else f"NASH_F{f}"


def stage_from_fibrosis_only(pheno: pd.DataFrame) -> pd.Series:
    """Ordered axis = fibrosis stage 0-4, for cohorts with no diagnosis label."""
    col = next((c for c in ("fibrosis_stage", "fibrosis", "stage")
                if c in pheno.columns), None)
    if col is None:
        raise ValueError("No fibrosis column; available: "
                         + ", ".join(pheno.columns))
    v = pd.to_numeric(pheno[col].astype(str).str.extract(r"(\d+)")[0],
                      errors="coerce")
    labels = [f"F{int(x)}" if pd.notna(x) and 0 <= x <= 4 else None for x in v]
    return pd.Categorical(labels, categories=AXES["fibrosis"], ordered=True)


def stage_from_subtype(pheno: pd.DataFrame) -> pd.Series:
    """Ordered axis = NAFL then NASH, for cohorts with no fibrosis stage."""
    col = next((c for c in ("disease_subtype", "subtype", "disease_state")
                if c in pheno.columns), None)
    if col is None:
        raise ValueError("No disease-subtype column; available: "
                         + ", ".join(pheno.columns))
    v = pheno[col].astype(str).str.strip().str.upper()
    labels = [("NASH" if x.startswith("NASH") else
               "NAFL" if x.startswith("NAFL") else None) for x in v]
    return pd.Categorical(labels, categories=AXES["subtype"], ordered=True)


def stage_from_tissue(pheno: pd.DataFrame) -> pd.Series:
    """Ordered axis = tissue type, by proximity to cancer (GSE164760)."""
    col = next((c for c in ("tissue", "source_tissue") if c in pheno.columns),
               None)
    if col is None:
        raise ValueError("No tissue column; available: "
                         + ", ".join(pheno.columns))
    v = pheno[col].astype(str).str.strip()
    known = set(AXES["tissue"])
    labels = [x if x in known else None for x in v]
    missing = sorted({x for x, lab in zip(v, labels) if lab is None})
    if missing:
        raise ValueError(f"Unrecognized tissue labels: {missing}")
    return pd.Categorical(labels, categories=AXES["tissue"], ordered=True)


def stage_from_tissue_hcc(pheno: pd.DataFrame) -> pd.Series:
    """Ordered axis = adjacent liver then tumor (GSE76427)."""
    col = next((c for c in ("tissue", "source_tissue") if c in pheno.columns),
               None)
    if col is None:
        raise ValueError("No tissue column; available: "
                         + ", ".join(pheno.columns))
    v = pheno[col].astype(str).str.strip()
    labels = [TISSUE_HCC_ALIASES.get(x.lower()) for x in v]
    missing = sorted({x for x, lab in zip(v, labels) if lab is None})
    if missing:
        raise ValueError(f"Unrecognized tissue labels: {missing}")
    return pd.Categorical(labels, categories=AXES["tissue_hcc"], ordered=True)


def build_stage_column(pheno: pd.DataFrame, axis: str = "auto"):
    """Return (categorical stage column, ordered category list, axis name)."""
    if axis == "tissue-hcc":
        return (stage_from_tissue_hcc(pheno), AXES["tissue_hcc"],
                "tissue_hcc")
    if axis == "tissue":
        return stage_from_tissue(pheno), AXES["tissue"], "tissue"
    if axis == "fibrosis":
        return stage_from_fibrosis_only(pheno), AXES["fibrosis"], "fibrosis"
    if axis == "subtype":
        return stage_from_subtype(pheno), AXES["subtype"], "subtype"
    return _build_group_axis(pheno), AXES["group"], "group"


def _build_group_axis(pheno: pd.DataFrame) -> pd.Series:
    """Best-effort stage assignment from whatever columns GEO supplied."""
    for col in ("group_in_paper", "group", "diagnosis", "disease_state"):
        if col in pheno.columns:
            s = pheno[col].map(normalize_stage)
            if s.notna().mean() > 0.8:
                return pd.Categorical(s, categories=STAGE_ORDER, ordered=True)

    dis = next((c for c in ("disease", "disease_state", "diagnosis")
                if c in pheno.columns), None)
    fib = next((c for c in ("fibrosis_stage", "fibrosis", "stage")
                if c in pheno.columns), None)
    if dis and fib:
        s = [stage_from_fibrosis(a, b)
             for a, b in zip(pheno[dis], pheno[fib])]
        return pd.Categorical(s, categories=STAGE_ORDER, ordered=True)

    raise ValueError(
        "Could not derive a stage column. Available columns:\n  "
        + "\n  ".join(pheno.columns)
        + "\nInspect the series matrix and extend build_stage_column()."
    )


def numeric(pheno: pd.DataFrame, *candidates: str) -> pd.Series | None:
    for c in candidates:
        if c in pheno.columns:
            return pd.to_numeric(
                pheno[c].astype(str).str.extract(r"([-\d.]+)")[0],
                errors="coerce")
    return None
