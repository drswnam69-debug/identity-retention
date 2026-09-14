"""survival.py -- Cox proportional hazards, Kaplan-Meier and log-rank,
without SciPy or lifelines.

The pipeline runs on a stock python3 with numpy and pandas; the prognostic
step must not break that. Everything here is validated against lifelines in
code/test_survival.py.
"""

from __future__ import annotations

import math
import numpy as np

from stats_lite import norm_two_sided, rankdata


def _order_by_time(time, event):
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    idx = np.argsort(t, kind="mergesort")
    return t[idx], e[idx], idx


def cox_ph(time, event, X, names=None, max_iter: int = 100,
           tol: float = 1e-9) -> dict:
    """Cox model with Breslow handling of ties, by Newton-Raphson.

    Returns coefficients on the log-hazard scale, standard errors from the
    inverse observed information, Wald p-values, and hazard ratios.
    """
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    X = np.atleast_2d(np.asarray(X, float))
    if X.shape[0] != len(t):
        X = X.T
    ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(X).all(axis=1) & (t > 0)
    t, e, X = t[ok], e[ok], X[ok]
    n, p = X.shape
    if n == 0 or p == 0:
        raise ValueError("no usable rows for the Cox model")

    order = np.argsort(-t, kind="mergesort")   # descending: risk set accrues
    t, e, X = t[order], e[order], X[order]

    beta = np.zeros(p)
    loglik = float("nan")
    for _ in range(max_iter):
        eta = np.clip(X @ beta, -60, 60)
        w = np.exp(eta)
        # cumulative sums over the risk set {j : t_j >= t_i}
        S0 = np.cumsum(w)
        S1 = np.cumsum(w[:, None] * X, axis=0)
        S2 = np.cumsum(w[:, None, None] * X[:, :, None] * X[:, None, :],
                       axis=0)
        # ties: every event at the same time shares the same risk set
        last = np.zeros(len(t), int)
        i = 0
        while i < len(t):
            j = i
            while j + 1 < len(t) and t[j + 1] == t[i]:
                j += 1
            last[i:j + 1] = j
            i = j + 1
        S0e, S1e, S2e = S0[last], S1[last], S2[last]

        d = e.astype(bool)
        if not d.any():
            raise ValueError("no events")
        loglik = float((eta[d] - np.log(S0e[d])).sum())
        Z = S1e[d] / S0e[d][:, None]
        grad = (X[d] - Z).sum(axis=0)
        info = np.zeros((p, p))
        for k in np.flatnonzero(d):
            zk = S1e[k] / S0e[k]
            info += S2e[k] / S0e[k] - np.outer(zk, zk)
        try:
            step = np.linalg.solve(info, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(info) @ grad
        beta_new = beta + step
        if not np.all(np.isfinite(beta_new)):
            break
        if np.max(np.abs(step)) < tol:
            beta = beta_new
            break
        beta = beta_new

    cov = np.linalg.pinv(info)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = beta / se
    return {
        "names": list(names) if names is not None
        else [f"x{i}" for i in range(p)],
        "coef": beta, "se": se, "z": z,
        "p": np.array([norm_two_sided(float(v)) if np.isfinite(v) else np.nan
                       for v in z]),
        "hr": np.exp(beta),
        "hr_lo": np.exp(beta - 1.959963985 * se),
        "hr_hi": np.exp(beta + 1.959963985 * se),
        "n": int(n), "n_events": int(e.sum()), "loglik": loglik,
    }


def kaplan_meier(time, event):
    """Returns (times, survival, at_risk) for a single group."""
    t, e, _ = _order_by_time(time, event)
    times, surv, at_risk = [], [], []
    s = 1.0
    n = len(t)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and t[j + 1] == t[i]:
            j += 1
        d = float(e[i:j + 1].sum())
        risk = n - i
        if d > 0:
            s *= 1.0 - d / risk
        times.append(float(t[i])); surv.append(s); at_risk.append(risk)
        i = j + 1
    return np.array(times), np.array(surv), np.array(at_risk)


def median_survival(time, event) -> float:
    times, surv, _ = kaplan_meier(time, event)
    below = np.flatnonzero(surv <= 0.5)
    return float(times[below[0]]) if len(below) else float("nan")


def logrank(time, event, group) -> dict:
    """k-sample log-rank test (chi-square on k-1 df, normal-approx p for k=2)."""
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    g = np.asarray(group)
    ok = np.isfinite(t) & np.isfinite(e)
    t, e, g = t[ok], e[ok], g[ok]
    levels = list(dict.fromkeys(g.tolist()))
    k = len(levels)
    obs = np.array([e[g == L].sum() for L in levels], float)
    exp = np.zeros(k)
    var = np.zeros(k)
    for ut in np.unique(t[e == 1]):
        at_risk = np.array([float((t[g == L] >= ut).sum()) for L in levels])
        N = at_risk.sum()
        d = float(((t == ut) & (e == 1)).sum())
        if N <= 1 or d == 0:
            continue
        exp += d * at_risk / N
        var += d * (at_risk / N) * (1 - at_risk / N) * (N - d) / (N - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = float(np.nansum((obs - exp) ** 2 / np.where(exp > 0, exp, np.nan)))
    if k == 2 and var[0] > 0:
        z = (obs[0] - exp[0]) / math.sqrt(var[0])
        p = norm_two_sided(z)
    else:
        p = _chi2_sf(chi2, k - 1)
    return {"levels": [str(L) for L in levels],
            "observed": obs.tolist(), "expected": np.round(exp, 3).tolist(),
            "chi2": round(chi2, 4), "df": k - 1, "p": p}


def _chi2_sf(x: float, df: int) -> float:
    """Upper tail of the chi-square distribution, via the regularized gamma."""
    if x <= 0:
        return 1.0
    a, xx = df / 2.0, x / 2.0
    if xx < a + 1:                        # series
        term = 1.0 / a
        total = term
        for k in range(1, 500):
            term *= xx / (a + k)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return 1.0 - total * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
    b, c = xx + 1.0 - a, 1e300           # continued fraction
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h


def schoenfeld_ph_test(time, event, x, beta: float) -> dict:
    """Correlation of scaled Schoenfeld residuals with rank(time).

    A significant correlation indicates the hazard ratio is not constant.
    """
    t = np.asarray(time, float)
    e = np.asarray(event, float).astype(bool)
    x = np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(x)
    t, e, x = t[ok], e[ok], x[ok]
    resid, times = [], []
    for ut in np.unique(t[e]):
        risk = t >= ut
        w = np.exp(beta * x[risk])
        xbar = float((w * x[risk]).sum() / w.sum())
        for xi in x[(t == ut) & e]:
            resid.append(xi - xbar)
            times.append(ut)
    if len(resid) < 5:
        return {"rho": float("nan"), "p": float("nan"), "n_events": len(resid)}
    r = rankdata(np.array(times))
    rho = float(np.corrcoef(rankdata(np.array(resid)), r)[0, 1])
    n = len(resid)
    t_stat = rho * math.sqrt((n - 2) / max(1e-12, 1 - rho ** 2))
    from stats_lite import _t_two_sided
    return {"rho": round(rho, 4), "p": float(_t_two_sided(t_stat, n - 2)),
            "n_events": n}
