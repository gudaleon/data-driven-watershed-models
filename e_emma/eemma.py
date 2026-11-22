"""
eemma.py -- Ensemble End-Member Mixing Analysis (EEMMA), Python port.

A faithful port of Kirchner's ``eemma.R`` (version 1.0, build 2023.03.27):

    Kirchner, J.W. (2023), Mixing models with multiple, overlapping, or
    incomplete end-members, quantified using time series of a single
    tracer, Geophysical Research Letters, 50, e2023GL104147.

    Kirchner, J.W. (2023), eemma.R, an R script for Ensemble End-Member
    Mixing Analysis. EnviDat. doi:10.16904/envidat.410

The original script is (C) 2023 ETH Zurich and James Kirchner, released
under the GNU General Public License 3 (GPL-3).  This port is a
derivative work and is likewise GPL-3.  See PORT_NOTES.md.

EEMMA estimates source (end-member) fractions in a mixture from the
time-series correlations among tracer time series of a *single* tracer,
rather than from the multi-tracer geometry of classical EMMA.  It
handles multiple, overlapping, or incomplete (partly unmeasured)
end-members, can estimate the tracer time series of an unmeasured
end-member, and can estimate a first-order "memory" (carry-over)
coefficient.  Mixing fractions are estimated either by ordinary least
squares (OLS) or by an errors-in-variables (EIV) procedure implementing
Fuller's (1987) bias correction, with analytic (delta-method / first
order second moment) propagation of the coefficient covariance matrix
through every back-transformation.

Public functions, mirroring the R script:

    eemma(y, X, ...)  -- ensemble end-member mixing analysis
    eiv(y, X, ...)    -- errors-in-variables regression (Fuller 1987)
    ols(y, X, ...)    -- ordinary least squares regression

Conventions follow the R original: inputs may be array-like or pandas
objects; column names are carried through to the outputs.  Outputs are
returned as dicts mirroring the R lists, with ``coeffs`` and
``timeseries`` as pandas DataFrames, ``f``/``se`` as pandas Series and
``f_cov``/``beta_cov`` as pandas DataFrames.
"""

import warnings

import numpy as np
import pandas as pd
from scipy import optimize, stats

__all__ = ["eemma", "eiv", "ols"]
__version__ = "1.0-port.1"


# ----------------------------------------------------------------------
# small helpers replicating R semantics
# ----------------------------------------------------------------------

def _match_arg(arg, choices):
    """Replicate R's match.arg: exact or unique-prefix match."""
    if arg is None:
        return choices[0]
    if arg in choices:
        return arg
    matches = [c for c in choices if c.startswith(arg)]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(f"argument must be one of {choices}, got {arg!r}")


def _as_matrix(a):
    """as.matrix() for vectors / DataFrames / Series / ndarrays."""
    if isinstance(a, pd.DataFrame):
        return a.to_numpy(dtype=float), list(a.columns)
    if isinstance(a, pd.Series):
        return a.to_numpy(dtype=float).reshape(-1, 1), [a.name]
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr, None


def _complete_cases(Z):
    return ~np.isnan(Z).any(axis=1)


def _cov_rows(Z):
    """R cov(): sample covariance (n-1 denominator), Z is n x p."""
    return np.cov(Z, rowvar=False, ddof=1)


def _sd_cols(X):
    """R apply(X, 2, sd)."""
    return np.std(X, axis=0, ddof=1)


def _lag(v):
    """dplyr::lag with default 1: first element NA."""
    out = np.full_like(v, np.nan, dtype=float)
    out[1:] = v[:-1]
    return out


# ----------------------------------------------------------------------
# find_root: smallest root of det(cov_mat - lambda * err_mat)
# (nested inside eiv() in the R script; module-level here)
# ----------------------------------------------------------------------

def _find_root(cov_mat, err_mat):
    def func(x):
        return np.linalg.det(cov_mat - x * err_mat)

    npts = 100
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.diag(err_mat) / np.diag(cov_mat)
        ulim = 1.0 / np.max(ratios[np.isfinite(ratios)]) \
            if np.any(np.isfinite(ratios)) else np.inf
    if ulim == 0:
        return 1.0

    x = np.linspace(0.0, ulim * 1.02, npts)
    f = np.array([func(xi) for xi in x])

    # zero-crossing root (Brent's method, as R's uniroot)
    if np.min(f) < 0:
        j = int(np.argmax(f < 0))          # first index below zero
        root = optimize.brentq(func, x[j - 1], x[j])
    else:
        root = ulim

    # asymptotic root: local minima below a threshold
    pos = f[f > 0]
    fthresh = max(0.01 * f[0], np.min(pos)) if pos.size else 0.01 * f[0]
    localmin = ulim
    for j in range(1, npts - 1):
        if f[j] <= fthresh and f[j] < f[j - 1] and f[j] < f[j + 1]:
            localmin = x[j]
            break

    return localmin if localmin < root else root


# ----------------------------------------------------------------------
# eiv: errors-in-variables regression (Fuller 1987, eqs. 2.5.13/2.5.14)
# ----------------------------------------------------------------------

def eiv(y, X, no_intercept=False, err_all=1.0, err_y=None, err_x=None,
        rxx=None, rxy=None, Szz=None, uncert="relative", alpha=1.0):
    """Errors-in-variables linear regression (Fuller's bias correction).

    Parameters mirror the R function of the same name in eemma.R.
    ``err_y`` defaults to ``err_all`` and ``err_x`` to ``err_all``
    repeated, exactly as in the R signature.  ``Szz``, if given, is the
    error covariance matrix of Z = [y, X] and supersedes
    err_y/err_x/rxx/rxy.
    """
    uncert = _match_arg(uncert, ["relative", "absolute"])

    y_mat, _ = _as_matrix(y)
    X_mat, xnames = _as_matrix(X)
    if y_mat.shape[0] != X_mat.shape[0]:
        raise ValueError("fatal error in eiv: X and y have different n")

    k = X_mat.shape[1]
    if xnames is None:
        xnames = [f"x{i + 1}" for i in range(k)]

    if err_y is None:
        err_y = err_all
    if err_x is None:
        err_x = np.full(k, float(err_all))
    else:
        err_x = np.atleast_1d(np.asarray(err_x, dtype=float))
    if rxx is None:
        rxx = np.eye(k)
    else:
        rxx = np.asarray(rxx, dtype=float)
    if rxy is None:
        rxy = np.zeros(k)
    else:
        rxy = np.atleast_1d(np.asarray(rxy, dtype=float)).copy()

    # retain only complete cases
    Z = np.hstack([y_mat, X_mat])
    cc = _complete_cases(Z)
    y_c = y_mat[cc]
    X_c = X_mat[cc]
    Z = Z[cc]
    n = int(cc.sum())
    if n < (k + 1):
        raise ValueError("Fatal error in eiv: too few complete cases -- "
                         "no degrees of freedom left!")

    if err_x.size == 1:
        err_x = np.full(X_c.shape[1], float(err_x[0]))
    elif err_x.size != X_c.shape[1]:
        raise ValueError("Fatal error in eiv: err.x must be a single value, "
                         "or a vector of the same length as NCOL(X)")

    # add dummy variable for intercept if desired
    if not no_intercept:
        k = k + 1
        X_c = np.hstack([X_c, np.ones((n, 1))])
        rxy = np.append(rxy, 0.0)
        rxx = np.block([[rxx, np.zeros((rxx.shape[0], 1))],
                        [np.zeros((1, rxx.shape[1] + 1))]])
        rxx[k - 1, k - 1] = 1.0
        err_x = np.append(err_x, 0.0)
        if Szz is not None:
            Szz = np.block([[Szz, np.zeros((Szz.shape[0], 1))],
                            [np.zeros((1, Szz.shape[1] + 1))]])
        xnames = list(xnames) + ["intercept"]
        Z = np.hstack([y_c, X_c])

    # covariances and mean crossproducts
    mZZ = _cov_rows(Z) * (n - 1) / n
    MXX = (X_c.T @ X_c) / n
    MXY = (X_c.T @ y_c) / n

    # error covariance matrix Szz
    if Szz is None:
        if np.linalg.det(rxx) <= 0:
            raise ValueError("Fatal error: error correlation matrix must "
                             "be positive definite!")
        if np.any(np.diag(rxx) != 1):
            raise ValueError("Diagonal elements of correlation matrix "
                             "must be 1!")
        if np.any(np.concatenate([[err_y], err_x]) < 0):
            raise ValueError("Fatal error: uncertainties err.y and err.x "
                             "cannot be negative!")
        err_z = np.concatenate([[err_y], err_x])
        rzz = np.block([[np.ones((1, 1)), rxy.reshape(1, -1)],
                        [rxy.reshape(-1, 1), rxx]])
        Szz = np.outer(err_z, err_z) * rzz
    Szz = np.asarray(Szz, dtype=float)

    Sxy = Szz[1:, 0].copy()          # first column, minus its first element
    Sxx = Szz[1:, 1:].copy()

    # constants in X distort the determinant search: drop them there
    mZZ_root, Szz_root = mZZ, Szz
    if np.any(np.diag(mZZ) == 0):
        keep = np.diag(mZZ) != 0
        mZZ_root = mZZ[np.ix_(keep, keep)]
        Szz_root = Szz[np.ix_(keep, keep)]

    lam = _find_root(mZZ_root, Szz_root)

    # Fuller's corrections to the crossproducts
    if uncert == "relative":                       # Fuller eq. 2.5.14
        adj_MXX = MXX - lam * (1.0 - alpha / n) * Sxx
        adj_MXY = MXY - lam * (1.0 - alpha / n) * Sxy.reshape(-1, 1)
    else:                                          # Fuller eq. 2.5.13
        adj_MXX = MXX - (lam - alpha / n) * Sxx
        adj_MXY = MXY - (lam - alpha / n) * Sxy.reshape(-1, 1)

    inv_adj_MXX = np.linalg.inv(adj_MXX)
    beta = (inv_adj_MXX @ adj_MXY).ravel()
    e = (y_c - X_c @ beta.reshape(-1, 1)).ravel()
    s2e = float(np.sum(e ** 2) / (n - k))
    beta_cov = s2e / n * inv_adj_MXX
    se = np.sqrt(np.diag(beta_cov))

    return {
        "beta": pd.Series(beta, index=xnames),
        "se": pd.Series(se, index=xnames),
        "s2e": s2e,
        "beta_cov": pd.DataFrame(beta_cov, index=xnames, columns=xnames),
        "e": e,
        "lambda": lam,
    }


# ----------------------------------------------------------------------
# ols: ordinary least squares regression
# ----------------------------------------------------------------------

def ols(y, X, no_intercept=False):
    """Ordinary least squares multiple regression (port of ols in eemma.R)."""
    y_mat, _ = _as_matrix(y)
    X_mat, xnames = _as_matrix(X)
    if y_mat.shape[0] != X_mat.shape[0]:
        raise ValueError("fatal error in ols: X and y have different n")

    k = X_mat.shape[1]
    if xnames is None:
        xnames = [f"x{i + 1}" for i in range(k)]

    Z = np.hstack([y_mat, X_mat])
    cc = _complete_cases(Z)
    y_c = y_mat[cc]
    X_c = X_mat[cc]
    n = int(cc.sum())
    if n < (k + 1):
        raise ValueError("Fatal error in ols: too few complete cases -- "
                         "no degrees of freedom left!")

    if not no_intercept:
        k = k + 1
        X_c = np.hstack([X_c, np.ones((n, 1))])
        xnames = list(xnames) + ["intercept"]

    MXX = (X_c.T @ X_c) / n
    MXY = (X_c.T @ y_c) / n

    inv_MXX = np.linalg.inv(MXX)
    beta = (inv_MXX @ MXY).ravel()
    e = (y_c - X_c @ beta.reshape(-1, 1)).ravel()
    s2e = float(np.sum(e ** 2) / (n - k))
    beta_cov = s2e / n * inv_MXX
    se = np.sqrt(np.diag(beta_cov))

    return {
        "beta": pd.Series(beta, index=xnames),
        "se": pd.Series(se, index=xnames),
        "s2e": s2e,
        "beta_cov": pd.DataFrame(beta_cov, index=xnames, columns=xnames),
        "e": e,
    }


# ----------------------------------------------------------------------
# eemma: ensemble end-member mixing analysis
# ----------------------------------------------------------------------

def eemma(y, X, complete=True, memory=False, no_intercept=None,
          method="OLS", err_all=1.0, err_y=None, err_x=None,
          sercorr_err_y=0.0, rxx=0.0, rxy=0.0, uncert="relative",
          alpha=1.0):
    """Ensemble end-member mixing analysis.

    Parameters mirror the R function of the same name in eemma.R:

    y : mixture tracer time series (length n).
    X : end-member tracer time series (n x k); column names are used as
        end-member names (default EM1, EM2, ...).
    complete : if True the end-member set is assumed complete and the
        mixing fractions are constrained to sum to 1; if False,
        unmeasured end-members are allowed and the unknown source
        ("unk"), its fraction, and its tracer time series are estimated.
    memory : if True, estimate a first-order memory coefficient
        ("lagmix"): the fraction of the mixture carried over from the
        previous time step.
    no_intercept : defaults to ``complete``, as in R.
    method : "OLS" (default) or "EIV".
    err_all / err_y / err_x : measurement/sampling error standard
        deviations (only used by EIV).
    sercorr_err_y : serial correlation of the *errors* in y.
    rxx / rxy : error correlation matrix among the x's / between y and
        each x; scalars are expanded as in R.
    uncert : "relative" (default) or "absolute".
    alpha : Fuller's method parameter (1-2 in eiv; positive in eemma).

    Returns a dict with keys: coeffs (DataFrame with f, se, t_statistic,
    p_value -- one-tailed, as in R), f, se, f_cov, RMSE, e, timeseries,
    lambda (EIV only, else NaN).
    """
    uncert = _match_arg(uncert, ["relative", "absolute"])
    method = _match_arg(method, ["OLS", "EIV"])
    if no_intercept is None:
        no_intercept = bool(complete)

    y_mat, _ = _as_matrix(y)
    X_mat, xnames_in = _as_matrix(X)
    if y_mat.shape[0] != X_mat.shape[0]:
        raise ValueError("fatal error in eemma: data series have different n")

    k = X_mat.shape[1]
    xnames = list(xnames_in) if xnames_in is not None \
        else [f"EM{i + 1}" for i in range(k)]

    if err_y is None:
        err_y = err_all
    if err_x is None:
        err_x = np.full(k, float(err_all))
    else:
        err_x = np.atleast_1d(np.asarray(err_x, dtype=float))
        if err_x.size == 1:
            err_x = np.full(k, float(err_x[0]))
        elif err_x.size != k:
            raise ValueError("Fatal error in eemma: err.x must be a single "
                             "value, or a vector of the same length as NCOL(X)")

    # correlation matrix of measurement errors in x's
    rxx_arr = np.asarray(rxx, dtype=float)
    if rxx_arr.size == 1:
        rxx_m = np.full((k, k), float(rxx_arr.ravel()[0]))
        np.fill_diagonal(rxx_m, 1.0)
    else:
        rxx_m = rxx_arr.copy()
    if rxx_m.shape != (k, k):
        raise ValueError("Fatal error in eemma: rxx has wrong size")
    if np.any(np.abs(rxx_m) > 1):
        raise ValueError("Fatal error in eemma: all values of rxx must be "
                         "between -1 and 1, inclusive")
    if not np.allclose(rxx_m, rxx_m.T):
        raise ValueError("Fatal error in eemma: rxx must be symmetric")
    if np.min(np.linalg.eigvalsh(rxx_m)) <= 0:
        raise ValueError("Fatal error in eemma: rxx must be positive definite")

    rxy_arr = np.atleast_1d(np.asarray(rxy, dtype=float))
    if rxy_arr.size == 1:
        rxy_v = np.full(k, float(rxy_arr[0]))
    else:
        rxy_v = rxy_arr.copy()
    if rxy_v.size != k:
        raise ValueError("Fatal error in eemma: rxy has wrong size")
    if np.any(np.abs(rxy_v) > 1):
        raise ValueError("Fatal error in eemma: all values of rxy must be "
                         "between -1 and 1, inclusive")
    if abs(sercorr_err_y) >= 1:
        raise ValueError("Fatal error in eemma: sercorr.err.y must be "
                         "between -1 and 1")
    if alpha <= 0:
        raise ValueError("Fatal error in eemma: alpha must be positive")

    # first-order memory term: lagged mixture as an extra end-member
    if memory:
        X_mat = np.hstack([X_mat, _lag(y_mat.ravel()).reshape(-1, 1)])
        xnames = list(xnames) + ["lagmix"]
        err_x = np.append(err_x, err_y)
        k = k + 1
        rxx_m = np.block([[rxx_m, (rxy_v * sercorr_err_y).reshape(-1, 1)],
                          [(rxy_v * sercorr_err_y).reshape(1, -1),
                           np.ones((1, 1))]])
        rxy_v = np.append(rxy_v, sercorr_err_y)

    # retain only complete cases
    Z = np.hstack([y_mat, X_mat])
    cc = _complete_cases(Z)
    y_c = y_mat[cc]
    X_c = X_mat[cc]
    n = int(cc.sum())
    yv = y_c.ravel()

    # error covariance matrix Szz for Z = [y, X]
    if np.linalg.det(rxx_m) <= 0:
        raise ValueError("Fatal error: error correlation matrix must be "
                         "positive definite!")
    if np.any(np.diag(rxx_m) != 1):
        raise ValueError("Diagonal elements of correlation matrix must be 1!")
    if np.any(np.concatenate([err_x, [err_y]]) < 0):
        raise ValueError("Fatal error: uncertainties err.y and err.x "
                         "cannot be negative!")
    err_z = np.concatenate([[err_y], err_x])
    rzz = np.block([[np.ones((1, 1)), rxy_v.reshape(1, -1)],
                    [rxy_v.reshape(-1, 1), rxx_m]])
    if np.min(np.linalg.eigvalsh(rzz)) <= 0:
        raise ValueError("Fatal error in eemma: rzz must be positive definite")
    Szz = np.outer(err_z, err_z) * rzz
    if np.min(np.linalg.eigvalsh(Szz)) <= 0:
        raise ValueError("Fatal error in eemma: Szz must be positive definite")

    # ------------------------------------------------------------------
    # complete == TRUE: difference transformation, regress, back-transform
    # ------------------------------------------------------------------
    if complete:
        if not no_intercept:
            warnings.warn("Warning: eemma option no.intercept==FALSE is "
                          "being ignored because it is inconsistent with "
                          "complete==TRUE")
        if X_c.shape[1] == 1:
            raise ValueError("Fatal error: complete==TRUE requires two or "
                             "more end-members, or one end-member if "
                             "memory==TRUE!")
        sds = _sd_cols(X_c)
        if np.sum(sds == 0) > 2:
            raise ValueError("Fatal error: if complete==TRUE, there cannot "
                             "be more than two constant end-members")

        # transformation: subtract the last X column from y and the rest
        tr = np.eye(k + 1)
        tr[k, :] = -1.0
        tr = tr[:, :k]                       # drop the last column
        Zc = np.hstack([y_c, X_c])
        trz = Zc @ tr
        trSzz = tr.T @ Szz @ tr

        if method == "EIV":
            q = eiv(trz[:, [0]], trz[:, 1:], no_intercept=True, Szz=trSzz,
                    uncert=uncert, alpha=alpha)
        else:
            q = ols(trz[:, [0]], trz[:, 1:], no_intercept=True)
        beta = q["beta"].to_numpy()
        beta_cov = q["beta_cov"].to_numpy()

        if memory:
            sumbeta = np.sum(beta)
            f = np.append(beta / sumbeta, 1.0 - sumbeta)
            # R: J <- matrix(-beta/sumbeta^2, nrow=k-1, ncol=k-1) fills
            # column-wise, so J[i, j] = -beta[i]/sumbeta^2 for every j
            J = np.tile((-beta / sumbeta ** 2).reshape(-1, 1), (1, k - 1))
            J = J + np.eye(k - 1) / sumbeta
            J = np.vstack([J, np.full((1, k - 1), -1.0)])
            f_cov = J @ beta_cov @ J.T
            se = np.sqrt(np.diag(f_cov))
        else:
            btr = tr[1:, 1:]
            f = btr @ beta
            f[k - 1] = 1.0 + f[k - 1]
            f_cov = btr @ beta_cov @ btr.T
            se = np.sqrt(np.diag(f_cov))

        timeseries = pd.DataFrame(np.hstack([y_c, X_c]),
                                  columns=["mixture"] + list(xnames))
        return _package(f, se, f_cov, q, n, k, xnames, timeseries, method)

    # ------------------------------------------------------------------
    # complete == FALSE
    # ------------------------------------------------------------------
    sds = _sd_cols(X_c)
    if no_intercept:
        if np.sum(sds == 0) > 1:
            raise ValueError("Fatal error: if complete==FALSE, and "
                             "no.intercept==TRUE, there can be a maximum "
                             "of one constant end-member")
    else:
        if np.sum(sds == 0) > 0:
            raise ValueError("Fatal error: if complete==FALSE and "
                             "no.intercept==FALSE, there can be no "
                             "constant end-members")

    ybar = float(np.mean(yv)) if no_intercept else 0.0
    y_ctr = (y_c - ybar)
    X_ctr = (X_c - ybar)

    if method == "EIV":
        q = eiv(y_ctr, X_ctr, no_intercept=no_intercept, Szz=Szz,
                uncert=uncert, alpha=alpha)
    else:
        q = ols(y_ctr, X_ctr, no_intercept=no_intercept)
    beta = q["beta"].to_numpy()
    beta_cov = q["beta_cov"].to_numpy()
    e = q["e"]

    if memory:
        if no_intercept:
            # beta[k-1] (last) is the memory coefficient
            f = beta / (1.0 - beta[k - 1])
            f[k - 1] = beta[k - 1]
            J = np.eye(k - 1) / (1.0 - beta[k - 1])
            J = np.hstack([J, (beta[:k - 1] / (1.0 - beta[k - 1]) ** 2)
                           .reshape(-1, 1)])
            J = np.vstack([J, np.append(np.zeros(k - 1), 1.0)])
            f_cov = J @ beta_cov @ J.T
            se = np.sqrt(np.diag(f_cov))
            timeseries = pd.DataFrame(np.hstack([y_c, X_c]),
                                      columns=["mixture"] + list(xnames))
            out_names = list(xnames)
        else:
            # beta[k-1] is the memory term, beta[k] is the intercept
            f = np.zeros(k + 2)
            f[:k - 1] = beta[:k - 1] / (1.0 - beta[k - 1])
            f[k - 1] = beta[k - 1]
            f[k] = beta[k]
            f[k + 1] = 1.0 - np.sum(beta[:k - 1]) / (1.0 - beta[k - 1])
            J = np.eye(k - 1) / (1.0 - beta[k - 1])
            J = np.hstack([J, (beta[:k - 1] / (1.0 - beta[k - 1]) ** 2)
                           .reshape(-1, 1)])
            J = np.vstack([J, np.append(np.zeros(k - 1), 1.0)])
            J = np.hstack([J, np.zeros((k, 1))])
            J = np.vstack([J, np.append(np.zeros(k), 1.0)])
            J = np.vstack([J, np.concatenate([
                np.full(k - 1, -1.0 / (1.0 - beta[k - 1])),
                [np.sum(beta[:k - 1]) / (1.0 - beta[k - 1]) ** 2, 0.0]])])
            f_cov = J @ beta_cov @ J.T
            se = np.sqrt(np.diag(f_cov))
            out_names = list(xnames) + ["intercept", "unk"]

            unk = ybar + (f[k] + e) / (f[k + 1] * (1.0 - f[k - 1]))
            du_dint = 1.0 / (f[k + 1] * (1.0 - f[k - 1]))
            du_dfu = (-(f[k] + e) / (f[k + 1] * (1.0 - f[k - 1])) ** 2
                      * (1.0 - f[k - 1]))
            du_dfm = ((f[k] + e) / (f[k + 1] * (1.0 - f[k - 1])) ** 2
                      * f[k + 1])
            unk_se = np.sqrt(
                (se[k - 1] * du_dfm) ** 2 + (se[k] * du_dint) ** 2
                + (se[k + 1] * du_dfu) ** 2
                + 2 * du_dfm * du_dint * f_cov[k - 1, k]
                + 2 * du_dfm * du_dfu * f_cov[k - 1, k + 1]
                + 2 * du_dint * du_dfu * f_cov[k, k + 1])
            timeseries = pd.DataFrame(
                np.column_stack([y_c, X_c, np.full(n, f[k] + ybar),
                                 unk, unk_se]),
                columns=["mixture"] + out_names + ["unk_std_err"])
    else:
        f = beta.copy()
        f_cov = beta_cov.copy()
        se = np.sqrt(np.diag(f_cov))
        out_names = list(xnames)
        if not no_intercept:
            out_names = list(xnames) + ["intercept", "unk"]
            f = np.append(f, 1.0 - np.sum(f[:k]))
            f_tr = np.vstack([np.eye(k + 1),
                              np.append(np.full(k, -1.0), 0.0)])
            f_cov = f_tr @ f_cov @ f_tr.T
            se = np.sqrt(np.diag(f_cov))
            unk = ybar + (f[k] + e) / f[k + 1]
            unk_se = np.sqrt(
                (se[k] / f[k + 1]) ** 2
                + (se[k + 1] * (f[k] + e) / f[k + 1] ** 2) ** 2
                - 2 * f_cov[k, k + 1] * (f[k] + e) / f[k + 1] ** 3)
            timeseries = pd.DataFrame(
                np.column_stack([y_c, X_c, np.full(n, f[k] + ybar),
                                 unk, unk_se]),
                columns=["mixture"] + out_names + ["unk_std_err"])
        else:
            timeseries = pd.DataFrame(np.hstack([y_c, X_c]),
                                      columns=["mixture"] + out_names)

    return _package(f, se, f_cov, q, n, k, out_names, timeseries, method)


def _package(f, se, f_cov, q, n, k, names, timeseries, method):
    """Shape the eemma return value (mirrors the R list)."""
    f = np.asarray(f, dtype=float)
    se = np.asarray(se, dtype=float)
    t_stat = np.abs(f / se)
    p_val = stats.t.sf(t_stat, df=n - k)   # one-tailed, as in R
    coeffs = pd.DataFrame(
        {"f": f, "se": se, "t_statistic": t_stat, "p_value": p_val},
        index=list(names))
    return {
        "coeffs": coeffs,
        "f": pd.Series(f, index=list(names)),
        "se": pd.Series(se, index=list(names)),
        "f_cov": pd.DataFrame(f_cov, index=list(names),
                              columns=list(names)),
        "RMSE": float(np.sqrt(q["s2e"])),
        "e": q["e"],
        "timeseries": timeseries,
        "lambda": q["lambda"] if method == "EIV" else np.nan,
    }
