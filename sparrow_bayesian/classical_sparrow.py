"""Classical (NLLS) SPARROW core, ported from USGS RSPARROW (public domain).

This module is a faithful Python port of the computational core of
RSPARROW 1.1.0 -- the parts that define the model -- isolated from the
R source in:

    R/estimateFeval.R      model evaluation: delivery, decay, residuals
    R/estimateFevalNoadj.R unconditioned variant of the same
    R/predict.R            prediction accumulation (ptnoder)
    R/estimateOptimize.R   NLLS estimation driver (NLMRT `nlfb`)
    R/createDataMatrix.R   data-matrix layout and delivery mean adjustment
    R/setNLLSWeights.R     regression weights (default = 1)
    src/tnoder.for         conditioned accumulation returning log residuals
    src/ptnoder.for        accumulation returning per-reach total loads

Everything else in RSPARROW (Shiny DSS, mapping, plotting, scenario UI,
report generation) is deliberately not ported; see PORT_NOTES.md.

Model equations, exactly as in estimateFeval.R / predict.R
-----------------------------------------------------------
For reach i with source inputs S[i, s], delivery variables D[i, d],
stream-decay variables T[i, k] and reservoir-decay variables R[i, m]:

  delivery factor for source s:
      F[i, s] = exp( sum_d design[s, d] * beta_dlv[d] * D_adj[i, d] )
      where D_adj is D with each column mean-centred over all reaches
      when mean adjustment is on (createDataMatrix.R, default in the
      RSPARROW tutorial models), and `design` is the 0/1 matrix from
      design_matrix.csv selecting which delivery variables act on
      which sources.

  delivered incremental source load:
      G[i] = sum_s beta_src[s] * S[i, s] * F[i, s]

  stream (reach) decay factor:   A[i] = prod_k exp(-T[i, k] * beta_str[k])
  reservoir decay factor:        B[i] = prod_m 1 / (1 + R[i, m] * beta_res[m])

  incremental delivered load:    inc[i]  = sqrt(A[i]) * B[i] * G[i]
      (the square root applies the reach's own decay to only half of
      the reach's travel, for load generated inside the reach)
  transport (carry) factor:      carry[i] = frac[i] * A[i] * B[i]

Network accumulation (tnoder.for / ptnoder.for): reaches are processed
in hydrologic (upstream-to-downstream) order.  A node accumulator
collects the load delivered to each node:

      load[i] = inc[i] + carry[i] * node[fnode[i]]
      node[tnode[i]] += iftran[i] * load[i]

Calibration (tnoder with ifadjust = 1): at a reach carrying an observed
load (depvar > 0) the residual is log(depvar) - log(load) -- with the
Fortran guard that a non-positive predicted load is replaced by 1.0
before taking the log -- and the observed load is then substituted for
the predicted load in everything accumulated downstream ("monitoring
adjustment" / conditioning).  Residuals are multiplied by sqrt(weight);
the default RSPARROW weights are all 1 (NLLS_weights = "default").

Estimation is nonlinear least squares on those log residuals, with
parameter bounds taken from parameters.csv.  RSPARROW drives this with
the NLMRT `nlfb` routine; here scipy.optimize.least_squares is used.

Parameter vector order follows RSPARROW's beta layout:
    [SOURCE..., DELIVF..., STRM..., RESV...]
with one data variable of the same name paired to each parameter.

Environment: numpy / scipy / pandas (mapvenv).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

__all__ = [
    "Parameter",
    "SparrowSpec",
    "SparrowNetwork",
    "ClassicalSparrow",
    "topological_order",
]


# ----------------------------------------------------------------------
# Parameter / model specification
# ----------------------------------------------------------------------
@dataclass
class Parameter:
    """One row of an RSPARROW parameters.csv file."""

    name: str
    parm_type: str  # SOURCE | DELIVF | STRM | RESV | OTHER
    init: float
    parm_min: float
    parm_max: float

    @property
    def is_constant(self) -> bool:
        """RSPARROW readParameters.R: parmConstant is 1 exactly when
        init == min == max and init != 0."""
        return self.init == self.parm_max == self.parm_min and self.init != 0.0

    @property
    def is_selected(self) -> bool:
        """RSPARROW selectParmValues.R `pselect`: a parameter enters the
        model if it is a nonzero constant, or it is estimable
        (parmMax > 0, or parmMax >= 0 for DELIVF) with min < max, and
        SOURCE parameters additionally require parmMin >= 0."""
        if self.is_constant:
            return True
        if not (
            self.parm_max > 0 or (self.parm_type == "DELIVF" and self.parm_max >= 0)
        ):
            return False
        if not self.parm_min < self.parm_max:
            return False
        if self.parm_type == "SOURCE" and self.parm_min < 0:
            return False
        return True


@dataclass
class SparrowSpec:
    """Model specification: which data variables play which roles.

    `design` is the (n_sources x n_delivery) 0/1 design matrix from
    design_matrix.csv; row order matches `sources`, column order
    matches `delivery_vars`.  If delivery_vars is empty, design may be
    None and every delivery factor is 1.
    """

    sources: list[str]
    delivery_vars: list[str] = field(default_factory=list)
    design: np.ndarray | None = None
    stream_decay_vars: list[str] = field(default_factory=list)
    reservoir_decay_vars: list[str] = field(default_factory=list)
    mean_adjust_delivery: bool = True

    def beta_names(self) -> list[str]:
        return (
            list(self.sources)
            + list(self.delivery_vars)
            + list(self.stream_decay_vars)
            + list(self.reservoir_decay_vars)
        )


# ----------------------------------------------------------------------
# Network topology
# ----------------------------------------------------------------------
def topological_order(fnode: np.ndarray, tnode: np.ndarray) -> np.ndarray:
    """Return reach indices in an upstream-to-downstream processing
    order (Kahn's algorithm on the reach graph).

    Reach j is upstream of reach i when tnode[j] == fnode[i].  The
    returned order guarantees every reach is processed after all
    reaches flowing into its from-node, which is exactly the property
    the Fortran accumulation relies on (RSPARROW obtains it from the
    HYDSEQ variable; the order produced here is a valid hydrologic
    order for the same network).
    """
    n = len(fnode)
    # reaches entering each node
    entering: dict[int, list[int]] = {}
    for j in range(n):
        entering.setdefault(int(tnode[j]), []).append(j)
    indegree = np.zeros(n, dtype=int)
    for i in range(n):
        indegree[i] = len(entering.get(int(fnode[i]), []))
    # downstream adjacency: reach -> reaches whose fnode == its tnode
    leaving: dict[int, list[int]] = {}
    for i in range(n):
        leaving.setdefault(int(fnode[i]), []).append(i)
    queue = [i for i in range(n) if indegree[i] == 0]
    order = []
    head = 0
    while head < len(queue):
        j = queue[head]
        head += 1
        order.append(j)
        for i in leaving.get(int(tnode[j]), []):
            indegree[i] -= 1
            if indegree[i] == 0:
                queue.append(i)
    if len(order) != n:
        raise ValueError(
            "reach network contains a cycle or disconnected bookkeeping: "
            f"ordered {len(order)} of {n} reaches"
        )
    return np.asarray(order, dtype=int)


class SparrowNetwork:
    """Reach network with compressed node indices and Stan-ready arrays.

    Attributes
    ----------
    order : (n_reach,) int
        Topological (upstream-to-downstream) processing order, as
        positions in the input frame.
    fnode_idx, tnode_idx : (n_reach,) int
        Node indices compressed to 0..n_node-1, in input row order.
    downstream_index : (n_reach,) int
        For each reach *in topological order*, the position (also in
        topological order) of the first reach whose fnode equals this
        reach's tnode, or -1 for terminal reaches.  This is the array
        a Stan implementation uses to push load downstream reach by
        reach; see PORT_NOTES.md.
    n_downstream : (n_reach,) int
        Number of reaches leaving this reach's tnode (>1 flags a
        flow split, handled in RSPARROW through `frac`).
    """

    def __init__(self, fnode: np.ndarray, tnode: np.ndarray):
        fnode = np.asarray(fnode, dtype=np.int64)
        tnode = np.asarray(tnode, dtype=np.int64)
        nodes, self.fnode_idx = np.unique(fnode, return_inverse=True)
        # compress both endpoint arrays over one shared node set
        all_nodes = np.unique(np.concatenate([fnode, tnode]))
        self.fnode_idx = np.searchsorted(all_nodes, fnode)
        self.tnode_idx = np.searchsorted(all_nodes, tnode)
        self.node_ids = all_nodes
        self.n_node = len(all_nodes)
        self.order = topological_order(fnode, tnode)
        # downstream reach lookup in topological-order positions
        pos_of_reach = np.empty(len(fnode), dtype=int)
        pos_of_reach[self.order] = np.arange(len(fnode))
        first_leaving: dict[int, int] = {}
        count_leaving: dict[int, int] = {}
        for i in range(len(fnode)):
            fn = int(fnode[i])
            count_leaving[fn] = count_leaving.get(fn, 0) + 1
            if fn not in first_leaving:
                first_leaving[fn] = i
        self.downstream_index = np.full(len(fnode), -1, dtype=int)
        self.n_downstream = np.zeros(len(fnode), dtype=int)
        for pos, i in enumerate(self.order):
            j = first_leaving.get(int(tnode[i]), -1)
            self.downstream_index[pos] = pos_of_reach[j] if j >= 0 else -1
            self.n_downstream[pos] = count_leaving.get(int(tnode[i]), 0)


# ----------------------------------------------------------------------
# The model
# ----------------------------------------------------------------------
class ClassicalSparrow:
    """SPARROW model over a reach table (pandas DataFrame).

    Required columns: fnode, tnode, frac, iftran, plus every variable
    named in the spec, plus the observed-load column (default
    "depvar"; values <= 0 or NaN mark non-monitoring reaches).

    Rows are reordered internally into hydrologic order: by the
    `hydseq` column when present (as RSPARROW does), otherwise by the
    topological order computed from fnode/tnode.  All per-reach
    outputs are returned in the *original* row order.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        spec: SparrowSpec,
        depvar_col: str = "depvar",
        weight: np.ndarray | None = None,
    ):
        self.spec = spec
        df = data.reset_index(drop=True).copy()
        self.n_reach = len(df)
        if "hydseq" in df.columns:
            sort_idx = np.argsort(df["hydseq"].to_numpy(), kind="stable")
        else:
            sort_idx = topological_order(
                df["fnode"].to_numpy(), df["tnode"].to_numpy()
            )
        self._sort = sort_idx
        self._unsort = np.empty(self.n_reach, dtype=int)
        self._unsort[sort_idx] = np.arange(self.n_reach)
        sdf = df.iloc[sort_idx].reset_index(drop=True)
        self.data = sdf

        self.network = SparrowNetwork(
            sdf["fnode"].to_numpy(), sdf["tnode"].to_numpy()
        )
        # fnode/tnode in compressed node space, hydrologic order
        self._fnode = self.network.fnode_idx
        self._tnode = self.network.tnode_idx
        self._frac = sdf["frac"].to_numpy(dtype=float)
        self._iftran = sdf["iftran"].to_numpy(dtype=float)
        dep = sdf[depvar_col].to_numpy(dtype=float) if depvar_col in sdf else np.zeros(self.n_reach)
        self._depvar = np.where(np.isnan(dep), 0.0, dep)
        self.is_site = self._depvar > 0
        self.n_site = int(self.is_site.sum())
        if weight is None:
            weight = np.ones(self.n_site)
        self._sqrt_weight = np.sqrt(np.asarray(weight, dtype=float))

        def _cols(names: list[str]) -> np.ndarray:
            if not names:
                return np.empty((self.n_reach, 0))
            return sdf[names].to_numpy(dtype=float)

        self._src = _cols(spec.sources)
        self._dlv_raw = _cols(spec.delivery_vars)
        if spec.mean_adjust_delivery and self._dlv_raw.shape[1] > 0:
            self._dlv = self._dlv_raw - self._dlv_raw.mean(axis=0)
        else:
            self._dlv = self._dlv_raw
        self._dec = _cols(spec.stream_decay_vars)
        self._res = _cols(spec.reservoir_decay_vars)
        if spec.delivery_vars:
            design = (
                np.asarray(spec.design, dtype=float)
                if spec.design is not None
                else np.ones((len(spec.sources), len(spec.delivery_vars)))
            )
            if design.shape != (len(spec.sources), len(spec.delivery_vars)):
                raise ValueError("design matrix shape does not match spec")
            self._design = design
        else:
            self._design = np.empty((len(spec.sources), 0))

        # beta layout offsets
        ns = len(spec.sources)
        nd = len(spec.delivery_vars)
        nk = len(spec.stream_decay_vars)
        self._sl_src = slice(0, ns)
        self._sl_dlv = slice(ns, ns + nd)
        self._sl_dec = slice(ns + nd, ns + nd + nk)
        self._sl_res = slice(ns + nd + nk, None)
        self.n_beta = ns + nd + nk + len(spec.reservoir_decay_vars)

    # -- model pieces -------------------------------------------------
    def delivery_factors(self, beta: np.ndarray) -> np.ndarray:
        """(n_reach x n_sources) source delivery factors F[i, s]."""
        if self._dlv.shape[1] == 0:
            return np.ones((self.n_reach, len(self.spec.sources)))
        lin = (self._dlv * beta[self._sl_dlv]) @ self._design.T
        return np.exp(lin)

    def decay_factors(self, beta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Stream decay factor A and reservoir decay factor B per reach."""
        if self._dec.shape[1] > 0:
            a = np.exp(-(self._dec * beta[self._sl_dec]).sum(axis=1))
        else:
            a = np.ones(self.n_reach)
        if self._res.shape[1] > 0:
            b = 1.0 / (1.0 + (self._res * beta[self._sl_res]).sum(axis=1))
        else:
            b = np.ones(self.n_reach)
        return a, b

    def incremental_loads(
        self, beta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (inc, carry, delivered_source_load G) in hydrologic order.

        inc[i]   = sqrt(A[i]) * B[i] * G[i]
        carry[i] = frac[i] * A[i] * B[i]
        """
        beta = np.asarray(beta, dtype=float)
        f = self.delivery_factors(beta)
        g = (f * self._src * beta[self._sl_src]).sum(axis=1)
        a, b = self.decay_factors(beta)
        inc = np.sqrt(a) * b * g
        carry = self._frac * a * b
        return (
            np.where(np.isnan(inc), 0.0, inc),
            np.where(np.isnan(carry), 0.0, carry),
            g,
        )

    # -- accumulation (tnoder / ptnoder) -------------------------------
    def _accumulate(
        self, inc: np.ndarray, carry: np.ndarray, condition: bool
    ) -> tuple[np.ndarray, np.ndarray]:
        """Core loop, a line-for-line port of tnoder.for / ptnoder.for.

        Returns (total_load, log_residuals_at_sites) in hydrologic
        order.  With condition=True the observed load replaces the
        predicted load at monitoring reaches before it is accumulated
        downstream, and total_load at a monitoring reach is the
        observed load (ptnoder with ifadjust=1); residuals are then
        log(depvar) - log(predicted load before substitution).
        """
        node = np.zeros(self.network.n_node)
        total = np.empty(self.n_reach)
        resid = np.empty(self.n_site)
        i_obs = 0
        depvar = self._depvar
        iftran = self._iftran
        for i in range(self.n_reach):
            rchld = inc[i] + carry[i] * node[self._fnode[i]]
            if depvar[i] > 0:
                resid[i_obs] = np.log(depvar[i]) - np.log(
                    rchld if rchld > 0 else 1.0
                )
                i_obs += 1
                if condition:
                    rchld = depvar[i]
            total[i] = rchld
            node[self._tnode[i]] += iftran[i] * rchld
        return total, resid

    def predict_loads(self, beta: np.ndarray, condition: bool = False) -> np.ndarray:
        """Total load at every reach, in the original row order.

        condition=False reproduces RSPARROW's `pload_total`
        (ptnoder, ifadjust=0); condition=True reproduces
        `mpload_total` (ifadjust=1).
        """
        inc, carry, _ = self.incremental_loads(beta)
        total, _ = self._accumulate(inc, carry, condition)
        return total[self._unsort]

    def residuals(self, beta: np.ndarray) -> np.ndarray:
        """Weighted log residuals at monitoring reaches, in hydrologic
        order -- the vector RSPARROW's NLLS minimises (estimateFeval)."""
        inc, carry, _ = self.incremental_loads(beta)
        _, resid = self._accumulate(inc, carry, condition=True)
        return self._sqrt_weight * resid

    def sse(self, beta: np.ndarray) -> float:
        r = self.residuals(beta)
        return float(r @ r)

    # -- estimation ----------------------------------------------------
    def fit(
        self,
        parameters: list[Parameter],
        x_scale: str | float = "jac",
        **least_squares_kwargs,
    ) -> dict:
        """NLLS estimation over the non-constant selected parameters.

        `parameters` is the full parameter list in beta order (sources,
        delivery, stream decay, reservoir).  Constant parameters are
        held at their initial values; non-selected parameters must not
        be passed (build the spec from selected parameters only).
        Bounds follow estimateOptimize.R: a bound that does not bracket
        the initial value is relaxed (min -> init/10 if min > init;
        max -> init*10 if max < init).
        """
        if len(parameters) != self.n_beta:
            raise ValueError(
                f"expected {self.n_beta} parameters in beta order, "
                f"got {len(parameters)}"
            )
        beta_full = np.array([p.init for p in parameters], dtype=float)
        free = np.array([not p.is_constant for p in parameters])
        free_idx = np.where(free)[0]
        x0 = beta_full[free_idx]
        lo = np.empty(len(free_idx))
        hi = np.empty(len(free_idx))
        for j, i in enumerate(free_idx):
            p = parameters[i]
            lo[j] = p.parm_min if p.parm_min <= p.init else p.init / 10.0
            hi[j] = p.parm_max if p.parm_max >= p.init else p.init * 10.0

        def fun(x: np.ndarray) -> np.ndarray:
            b = beta_full.copy()
            b[free_idx] = x
            return self.residuals(b)

        out = least_squares(
            fun, x0, bounds=(lo, hi), x_scale=x_scale, **least_squares_kwargs
        )
        beta_hat = beta_full.copy()
        beta_hat[free_idx] = out.x
        resid = out.fun
        dof = max(self.n_site - int(free.sum()), 1)
        mse = float(resid @ resid) / dof
        # covariance from the Jacobian at the solution, as in
        # estimateNLLSmetrics.R: cov = MSE * (J'J)^-1
        se = np.full(self.n_beta, np.nan)
        try:
            cov = mse * np.linalg.inv(out.jac.T @ out.jac)
            se[free_idx] = np.sqrt(np.diag(cov))
        except np.linalg.LinAlgError:
            pass
        return {
            "beta": beta_hat,
            "se": se,
            "free": free,
            "sse": float(resid @ resid),
            "mse": mse,
            "rmse": float(np.sqrt(mse)),
            "n_site": self.n_site,
            "n_parm": int(free.sum()),
            "residuals": resid,
            "scipy_result": out,
        }
