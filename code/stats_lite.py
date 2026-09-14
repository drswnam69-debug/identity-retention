"""stats_lite.py -- the statistics this analysis needs, without SciPy.

SciPy is used when it is installed (for exact t and normal tails), but every
function falls back to a pure-numpy/math implementation so the pipeline runs
on a stock python3 with numpy and pandas only.
"""

from __future__ import annotations

import math
import numpy as np

try:                                  # optional, only for exact tails
    from scipy import stats as _sp
except Exception:                     # pragma: no cover
    _sp = None


# --------------------------------------------------------------- tails
def norm_sf(z: float) -> float:
    """Upper tail of the standard normal."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def norm_two_sided(z: float) -> float:
    return 2.0 * norm_sf(abs(z))


def _t_two_sided(t: float, df: float) -> float:
    if _sp is not None:
        return float(2 * _sp.t.sf(abs(t), df))
    # regularized incomplete beta, so the t tail is exact without SciPy
    x = df / (df + t * t)
    return float(_betainc(df / 2.0, 0.5, x))


def _betacf(a: float, b: float, x: float, itmax: int = 300,
            eps: float = 3e-14) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbeta) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbeta) * _betacf(b, a, 1.0 - x) / b


# --------------------------------------------------------------- ranks
def rankdata(a: np.ndarray) -> np.ndarray:
    """Average ranks, ties shared (equivalent to scipy 'average')."""
    a = np.asarray(a, float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), float)
    sa = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def _tie_groups(a: np.ndarray) -> np.ndarray:
    _, counts = np.unique(np.asarray(a, float), return_counts=True)
    return counts


# --------------------------------------------------------------- tests
def mannwhitney_u(x, y):
    """Two-sided Mann-Whitney U with tie correction. Returns (U, z, p)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan"), float("nan")
    allv = np.concatenate([x, y])
    r = rankdata(allv)
    u1 = r[:n1].sum() - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    t = _tie_groups(allv)
    n = n1 + n2
    tie_term = (t ** 3 - t).sum()
    var = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)))
    if var <= 0:
        return float(u1), float("nan"), float("nan")
    z = (u1 - mu) / math.sqrt(var)
    return float(u1), float(z), float(norm_two_sided(z))


def jonckheere_terpstra(groups, alternative: str = "increasing"):
    """Jonckheere-Terpstra test for an ordered alternative.

    `groups` is a sequence of 1-D arrays given in the hypothesized order.
    Uses the tie-corrected large-sample variance. Returns (JT, z, p).
    """
    gs = [np.asarray(g, float) for g in groups]
    gs = [g[np.isfinite(g)] for g in gs]
    gs = [g for g in gs if len(g)]
    k = len(gs)
    if k < 2:
        return float("nan"), float("nan"), float("nan")

    jt = 0.0
    for i in range(k - 1):
        for j in range(i + 1, k):
            a, b = gs[i][:, None], gs[j][None, :]
            jt += float((a < b).sum() + 0.5 * (a == b).sum())

    n = np.array([len(g) for g in gs], float)
    N = n.sum()
    allv = np.concatenate(gs)
    t = _tie_groups(allv).astype(float)

    mean = (N ** 2 - (n ** 2).sum()) / 4.0

    if len(t) == len(allv):                      # no ties: closed form
        var = (N ** 2 * (2 * N + 3) - (n ** 2 * (2 * n + 3)).sum()) / 72.0
    else:
        term1 = (N * (N - 1) * (2 * N + 5)
                 - (n * (n - 1) * (2 * n + 5)).sum()
                 - (t * (t - 1) * (2 * t + 5)).sum()) / 72.0
        term2 = ((n * (n - 1) * (n - 2)).sum() * (t * (t - 1) * (t - 2)).sum()
                 / (36.0 * N * (N - 1) * (N - 2)))
        term3 = ((n * (n - 1)).sum() * (t * (t - 1)).sum()
                 / (8.0 * N * (N - 1)))
        var = term1 + term2 + term3

    if var <= 0:
        return float(jt), float("nan"), float("nan")
    z = (jt - mean) / math.sqrt(var)
    if alternative == "increasing":
        p = norm_sf(z)
    elif alternative == "decreasing":
        p = norm_sf(-z)
    else:
        p = norm_two_sided(z)
    return float(jt), float(z), float(p)


def jt_permutation_p(groups, n_perm: int = 20000, seed: int = 0,
                     alternative: str = "increasing") -> float:
    """Exact-ish permutation p-value for JT; used to validate the z-based p."""
    gs = [np.asarray(g, float) for g in groups]
    sizes = [len(g) for g in gs]
    allv = np.concatenate(gs)
    obs, _, _ = jonckheere_terpstra(gs)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(allv)
        out, idx = [], 0
        for s in sizes:
            out.append(perm[idx:idx + s]); idx += s
        stat, _, _ = jonckheere_terpstra(out)
        if alternative == "increasing":
            count += stat >= obs
        elif alternative == "decreasing":
            count += stat <= obs
        else:
            count += abs(stat - _jt_mean(sizes)) >= abs(obs - _jt_mean(sizes))
    return (count + 1) / (n_perm + 1)


def _jt_mean(sizes):
    n = np.array(sizes, float)
    return (n.sum() ** 2 - (n ** 2).sum()) / 4.0


def spearman(x, y):
    """Spearman rho with a t-approximation p-value. Returns (rho, p)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")
    rx, ry = rankdata(x), rankdata(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    if abs(rho) >= 1.0:
        return rho, 0.0
    t = rho * math.sqrt((n - 2) / (1 - rho ** 2))
    return rho, float(_t_two_sided(t, n - 2))


def benjamini_hochberg(p):
    """BH-FDR adjusted p-values, order preserved."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def ols(y, X, names=None):
    """Least squares with heteroscedasticity-naive SEs. Returns a dict."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    n, p = X.shape
    if n <= p:
        raise ValueError(f"n={n} observations for p={p} parameters")
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = n - p
    sigma2 = float(resid @ resid) / df
    xtx_inv = np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(xtx_inv) * sigma2)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = beta / se
    pvals = [_t_two_sided(float(ti), df) if np.isfinite(ti) else np.nan
             for ti in t]
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {
        "names": list(names) if names is not None
        else [f"x{i}" for i in range(p)],
        "beta": beta, "se": se, "t": t, "p": np.array(pvals),
        "n": int(n), "df": int(df),
        "r2": float(1 - (resid @ resid) / ss_tot) if ss_tot > 0 else np.nan,
    }
