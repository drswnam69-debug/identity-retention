"""rsi.py -- compute the Reductive-Supply Index and its sub-scores."""

from __future__ import annotations

import numpy as np
import pandas as pd

import rsi_config as cfg


def module_zscores(expr: pd.DataFrame, genes: list[str]
                   ) -> tuple[pd.Series, list[str], list[str]]:
    """Mean within-cohort z-score across the genes of one module.

    `expr` is genes x samples on a log scale. Z-scoring is done per gene
    across the samples of THIS cohort only, which is what makes the index
    comparable between platforms.
    """
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if not present:
        raise ValueError(f"None of the module genes are present: {genes}")
    sub = expr.loc[present]
    sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
    z = z.dropna(how="all")
    return z.mean(axis=0), present, missing


def compute_rsi(expr: pd.DataFrame, include_cholesterol: bool = False
                ) -> tuple[pd.DataFrame, dict]:
    """Return a per-sample table of RSI and its components, plus a report."""
    supply = list(cfg.MODULE_SUPPLY)
    if include_cholesterol:
        supply = supply + list(cfg.SENSITIVITY_S2_CHOLESTEROL)

    z_s, p_s, m_s = module_zscores(expr, supply)
    z_r, p_r, m_r = module_zscores(expr, cfg.MODULE_REDUCTION)
    z_d, p_d, m_d = module_zscores(expr, cfg.MODULE_DRAIN)

    rsi = (cfg.WEIGHT_SUPPLY * z_s
           + cfg.WEIGHT_REDUCTION * z_r
           - cfg.WEIGHT_DRAIN * z_d)

    out = pd.DataFrame({"z_supply": z_s, "z_reduction": z_r,
                        "z_drain": z_d, "RSI": rsi})

    for label, genes in (("shared", cfg.MODULE_SHARED),
                         ("effector", cfg.MODULE_EFFECTOR),
                         ("upstream", cfg.MODULE_UPSTREAM)):
        try:
            z, _, _ = module_zscores(expr, genes)
            out[f"z_{label}"] = z
        except ValueError:
            out[f"z_{label}"] = np.nan

    for g in cfg.FOCUS_GENES + cfg.MODULE_UPSTREAM:
        if g in expr.index:
            v = expr.loc[g]
            v = v.iloc[0] if isinstance(v, pd.DataFrame) else v
            out[f"gene_{g}"] = (v - v.mean()) / (v.std(ddof=1) or np.nan)

    report = {
        "n_samples": int(expr.shape[1]),
        "n_genes_in_matrix": int(expr.shape[0]),
        "supply_present": p_s, "supply_missing": m_s,
        "reduction_present": p_r, "reduction_missing": m_r,
        "drain_present": p_d, "drain_missing": m_d,
        "include_cholesterol_sensitivity": include_cholesterol,
        "lock_date": cfg.LOCK_DATE,
        "lock_hash": cfg.lock_hash(),
    }
    return out, report


def cronbach_alpha(expr: pd.DataFrame, genes: list[str]) -> float:
    """Internal consistency of a module (higher = genes co-vary)."""
    present = [g for g in genes if g in expr.index]
    if len(present) < 2:
        return float("nan")
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(
        sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna(how="all")
    k = z.shape[0]
    item_var = z.var(axis=1, ddof=1).sum()
    total_var = z.sum(axis=0).var(ddof=1)
    if total_var == 0:
        return float("nan")
    return float(k / (k - 1) * (1 - item_var / total_var))


def pc1_variance_explained(expr: pd.DataFrame, genes: list[str]) -> float:
    """Fraction of variance carried by the first principal component."""
    present = [g for g in genes if g in expr.index]
    if len(present) < 2:
        return float("nan")
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(
        sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0).dropna(how="all")
    if z.shape[0] < 2:
        return float("nan")
    sv = np.linalg.svd(z.values - z.values.mean(axis=1, keepdims=True),
                       compute_uv=False)
    return float(sv[0] ** 2 / (sv ** 2).sum())
