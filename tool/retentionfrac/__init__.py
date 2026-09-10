"""retentionfrac: how much of a paired tumor-adjacent contrast is not
dedifferentiation.

One function. It refuses to return a number when the premise the measure depends
on does not hold in the matrix you gave it, because that is the failure mode the
study behind it ran into and could not have caught any other way.

    from retentionfrac import identity_retention
    r = identity_retention(expr, tumor, adjacent, signature, identity_genes)
    print(r.retention, r.ci95, r.premise_ok)

Nam SW. How much of a tumor expression signature is dedifferentiation?
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

__all__ = ["identity_retention", "RetentionResult"]
__version__ = "1.0.1"


@dataclass
class RetentionResult:
    retention: float | None
    ci95: tuple[float, float] | None
    unadjusted_shift: float
    intercept: float
    n_pairs: int
    premise_ok: bool
    premise_delta: float
    premise_p: float
    stable: bool
    n_signature_genes_found: int
    n_identity_genes_found: int
    notes: list[str] = field(default_factory=list)

    def __repr__(self) -> str:
        if self.retention is None:
            return f"<RetentionResult NOT REPORTABLE: {'; '.join(self.notes)}>"
        lo, hi = self.ci95
        return (f"<RetentionResult {self.retention:.3f} ({lo:.3f} to {hi:.3f}), "
                f"{self.n_pairs} pairs>")


def _zmean(mat, genes, index):
    g = [x for x in genes if x in index]
    if not g:
        return None, 0
    sub = np.asarray([mat[x] for x in g], float)
    sd = sub.std(axis=1, ddof=1)
    keep = sd > 0
    sub, g = sub[keep], [x for x, k in zip(g, keep) if k]
    if not len(g):
        return None, 0
    z = (sub - sub.mean(axis=1, keepdims=True)) / sub.std(axis=1, ddof=1, keepdims=True)
    return z.mean(axis=0), len(g)


def _wilcoxon_p(d):
    d = np.asarray(d, float); d = d[d != 0]
    n = len(d)
    if n < 10:
        return float("nan")
    r = np.argsort(np.argsort(np.abs(d))) + 1.0
    w = float(r[d > 0].sum())
    mu, sig = n * (n + 1) / 4.0, np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    from math import erfc, sqrt
    return float(erfc(abs(w - mu) / (sig * sqrt(2.0))))


def _ols_intercept(y, X):
    return float(np.linalg.solve(X.T @ X, X.T @ y)[0])


def identity_retention(expr, tumor_samples: Sequence, adjacent_samples: Sequence,
                       signature_genes: Sequence[str],
                       identity_genes: Sequence[str],
                       composition_genes: Sequence[str] | None = None,
                       n_boot: int = 2000, seed: int = 0,
                       min_genes: int = 10, min_identity_genes: int = 8):
    """Identity-retention fraction for one signature in one paired cohort.

    expr              : mapping gene -> sequence of values, or a pandas DataFrame
                        with genes as the index and samples as columns, already on
                        a log scale and normalized as you would analyze it.
    tumor_samples,
    adjacent_samples  : equal-length, patient-matched, in the same order.
    signature_genes   : the gene set to measure.
    identity_genes    : the differentiated-cell covariate for this tissue. It must
                        share no gene with the signature or with composition_genes,
                        and it must not contain redox enzymes if what you are
                        measuring is redox.
    composition_genes : optional non-parenchymal covariate. Adjusting for it is a
                        different operation and is not a substitute.

    Returns a RetentionResult. `retention` is None, with the reason in `notes`,
    when the premise fails or the ratio has no stable denominator.
    """
    try:
        idx = set(expr.index); cols = list(expr.columns)
        get = lambda g: expr.loc[g].to_numpy(float)          # noqa: E731
        pos = {c: i for i, c in enumerate(cols)}
    except AttributeError:
        idx = set(expr); get = lambda g: np.asarray(expr[g], float)  # noqa: E731
        pos = None
    if len(tumor_samples) != len(adjacent_samples):
        raise ValueError("tumor_samples and adjacent_samples must be paired and "
                         "in the same order")
    n = len(tumor_samples)
    notes: list[str] = []

    class _M:
        def __getitem__(self, g):
            return get(g)
    M = _M()
    ti = ([pos[s] for s in tumor_samples] if pos else list(tumor_samples))
    ai = ([pos[s] for s in adjacent_samples] if pos else list(adjacent_samples))

    sig, n_sig = _zmean(M, signature_genes, idx)
    ident, n_id = _zmean(M, identity_genes, idx)
    if sig is None or n_sig < min_genes:
        return RetentionResult(None, None, float("nan"), float("nan"), n, False,
                               float("nan"), float("nan"), False, n_sig, n_id,
                               [f"only {n_sig} signature genes on this platform, "
                                f"minimum {min_genes}"])
    if ident is None or n_id < min_identity_genes:
        return RetentionResult(None, None, float("nan"), float("nan"), n, False,
                               float("nan"), float("nan"), False, n_sig, n_id,
                               [f"only {n_id} identity genes on this platform, "
                                f"minimum {min_identity_genes}"])

    d_sig = sig[ti] - sig[ai]
    d_id = ident[ti] - ident[ai]
    premise_delta = float(d_id.mean())
    premise_p = _wilcoxon_p(d_id)
    premise_ok = bool(premise_delta < 0 and premise_p < 0.05)
    if not premise_ok:
        notes.append("PREMISE FAILS: the identity score is not lower in tumor in "
                     "this matrix, so there is no confound to remove and no "
                     "retention value from it is interpretable")

    cols_X = [np.ones(n), d_id]
    if composition_genes is not None:
        comp, n_comp = _zmean(M, composition_genes, idx)
        if comp is not None:
            cols_X.append(comp[ti] - comp[ai])
        else:
            notes.append("composition covariate had no genes on this platform "
                         "and was dropped")
    X = np.column_stack(cols_X)
    unadj = float(d_sig.mean())
    a0 = _ols_intercept(d_sig, X)

    # n_boot = 0 asks for the point estimate without an interval. It used to
    # reach np.percentile on an empty array and raise IndexError, which is a
    # poor answer to a reasonable request; the stability guard then falls back
    # to the sign test on the paired shift itself.
    if n_boot and n_boot > 0:
        rng = np.random.default_rng(seed)
        IDX = rng.integers(0, n, size=(n_boot, n))
        den = d_sig[IDX].mean(axis=1)
        d_lo, d_hi = np.percentile(den, [2.5, 97.5])
        stable = not (d_lo <= 0.0 <= d_hi)
    else:
        IDX = None
        se = float(np.std(d_sig, ddof=1) / np.sqrt(n)) if n > 1 else float("inf")
        stable = abs(float(d_sig.mean())) > 1.96 * se
        notes.append("no bootstrap was requested, so the interval is not "
                     "reported and stability is judged from the paired shift's "
                     "own standard error")
    if not stable:
        notes.append("the unadjusted shift is not reliably different from zero, "
                     "so the ratio has no stable interval; nothing to retain")
    ci = None
    if stable and IDX is not None:
        Xb = X[IDX]
        beta = np.linalg.solve(np.einsum("bij,bik->bjk", Xb, Xb),
                               np.einsum("bij,bi->bj", Xb, d_sig[IDX])[..., None])[..., 0]
        ratio = np.abs(beta[:, 0]) / np.abs(den)
        ci = tuple(float(v) for v in np.percentile(ratio, [2.5, 97.5]))

    ret = (abs(a0) / abs(unadj)) if (stable and premise_ok) else None
    if premise_ok and stable and ci and (ci[1] - ci[0]) > 0.5:
        notes.append("the interval is wider than 0.50; this cohort is too small "
                     "to place this signature relative to any threshold")
    return RetentionResult(ret, ci if (premise_ok and stable) else None,
                           unadj, a0, n, premise_ok, premise_delta, premise_p,
                           stable, n_sig, n_id, notes)
