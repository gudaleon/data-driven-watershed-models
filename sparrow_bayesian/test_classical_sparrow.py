"""Hand-computed verification of classical_sparrow on a synthetic network.

The 8-reach network, every parameter, and every expected load below
are worked out by hand in this file with plain Python arithmetic --
no function from classical_sparrow is used to produce an expected
value.  Run:

    ~/workspace/mapvenv/bin/python test_classical_sparrow.py

Network (fnode -> tnode), flow splits at node 5 (frac 0.75 / 0.25):

    R1: 1 -> 3      R2: 2 -> 3      R3: 3 -> 5   (monitoring site)
    R4: 4 -> 5      R5: 5 -> 6      R7: 5 -> 8
    R6: 6 -> 7      R8: 8 -> 7      (7 = outlet node)

Model: sources S1 (delivery-adjusted) and S2 (not adjusted);
one delivery variable D acting on S1 only; one stream-decay variable
T; one reservoir-decay variable Rv.  Delivery mean adjustment is OFF
for this hand check so the arithmetic stays transparent.
"""

import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from classical_sparrow import (  # noqa: E402
    ClassicalSparrow,
    Parameter,
    SparrowSpec,
)

# ----------------------------------------------------------------------
# Hand-set inputs
# ----------------------------------------------------------------------
# beta order: [S1, S2, DLV, STRM, RESV]
BETA = np.array([2.0, 0.5, 0.3, 0.4, 1.5])

# per reach: fnode, tnode, frac, iftran, S1, S2, D, T, Rv, depvar
ROWS = [
    # R1
    dict(fnode=1, tnode=3, frac=1.00, iftran=1, S1=10.0, S2=100.0, D=1.0, T=0.5, Rv=0.0, depvar=0.0),
    # R2
    dict(fnode=2, tnode=3, frac=1.00, iftran=1, S1=20.0, S2=0.0, D=2.0, T=0.25, Rv=0.0, depvar=0.0),
    # R3 (site, observed load 500)
    dict(fnode=3, tnode=5, frac=1.00, iftran=1, S1=5.0, S2=50.0, D=0.0, T=1.0, Rv=0.2, depvar=500.0),
    # R4
    dict(fnode=4, tnode=5, frac=1.00, iftran=1, S1=0.0, S2=200.0, D=1.0, T=0.0, Rv=0.0, depvar=0.0),
    # R5
    dict(fnode=5, tnode=6, frac=0.75, iftran=1, S1=8.0, S2=0.0, D=0.5, T=0.5, Rv=0.0, depvar=0.0),
    # R6
    dict(fnode=6, tnode=7, frac=1.00, iftran=1, S1=4.0, S2=20.0, D=1.5, T=0.75, Rv=0.1, depvar=0.0),
    # R7
    dict(fnode=5, tnode=8, frac=0.25, iftran=1, S1=12.0, S2=10.0, D=0.0, T=0.5, Rv=0.0, depvar=0.0),
    # R8
    dict(fnode=8, tnode=7, frac=1.00, iftran=1, S1=6.0, S2=0.0, D=1.0, T=0.25, Rv=0.0, depvar=0.0),
]
df = pd.DataFrame(ROWS)

spec = SparrowSpec(
    sources=["S1", "S2"],
    delivery_vars=["D"],
    design=np.array([[1.0], [0.0]]),  # delivery acts on S1 only
    stream_decay_vars=["T"],
    reservoir_decay_vars=["Rv"],
    mean_adjust_delivery=False,
)
model = ClassicalSparrow(df, spec)

# ----------------------------------------------------------------------
# Hand computation, reach by reach
#   G    = 2.0 * S1 * exp(0.3 * D) + 0.5 * S2
#   A    = exp(-0.4 * T)          (stream decay)
#   B    = 1 / (1 + 1.5 * Rv)     (reservoir decay)
#   inc  = sqrt(A) * B * G
#   load = inc + frac * A * B * (load arriving at fnode)
# ----------------------------------------------------------------------
def hand(r):
    G = 2.0 * r["S1"] * math.exp(0.3 * r["D"]) + 0.5 * r["S2"]
    A = math.exp(-0.4 * r["T"])
    B = 1.0 / (1.0 + 1.5 * r["Rv"])
    return math.sqrt(A) * B * G, r["frac"] * A * B

inc1, car1 = hand(ROWS[0]); load1 = inc1
inc2, car2 = hand(ROWS[1]); load2 = inc2
node3 = load1 + load2
inc3, car3 = hand(ROWS[2]); load3_pred = inc3 + car3 * node3
inc4, car4 = hand(ROWS[3]); load4 = inc4
node5_uncond = load3_pred + load4          # without conditioning
node5_cond = 500.0 + load4                 # R3 observed load substituted
inc5, car5 = hand(ROWS[4]); load5 = inc5 + car5 * node5_cond
inc7, car7 = hand(ROWS[6]); load7 = inc7 + car7 * node5_cond
node6 = load5
node8 = load7
inc6, car6 = hand(ROWS[5]); load6 = inc6 + car6 * node6
inc8, car8 = hand(ROWS[7]); load8 = inc8 + car8 * node8

expected_uncond = np.array([load1, load2, load3_pred, load4,
                            inc5 + car5 * node5_uncond, load6,
                            inc7 + car7 * node5_uncond, load8])
# note: load6/load8 above used conditioned node5 via load5/load7;
# rebuild the unconditioned chain explicitly:
load5u = inc5 + car5 * node5_uncond
load7u = inc7 + car7 * node5_uncond
expected_uncond = np.array([load1, load2, load3_pred, load4,
                            load5u, inc6 + car6 * load5u,
                            load7u, inc8 + car8 * load7u])
expected_cond = np.array([load1, load2, 500.0, load4, load5, load6, load7, load8])
expected_resid = np.array([math.log(500.0) - math.log(load3_pred)])

got_uncond = model.predict_loads(BETA, condition=False)
got_cond = model.predict_loads(BETA, condition=True)
got_resid = model.residuals(BETA)

for name, got, exp in [
    ("unconditioned loads", got_uncond, expected_uncond),
    ("conditioned loads", got_cond, expected_cond),
    ("log residual at R3", got_resid, expected_resid),
]:
    err = np.max(np.abs(got - exp) / np.maximum(np.abs(exp), 1e-12))
    print(f"{name}: max relative error vs hand computation = {err:.3e}")
    assert err < 1e-10, (name, got, exp)

# ----------------------------------------------------------------------
# NLLS recovery: observed loads generated from BETA at all 8 reaches
# (the hand-computed unconditioned chain), then refit.  With only the
# 3 sites R3/R6/R8 observed, 5 free parameters are underdetermined
# (a manifold of exact fits), so the recovery gate uses all reaches.
# ----------------------------------------------------------------------
df2 = df.copy()
df2["depvar"] = expected_uncond
model2 = ClassicalSparrow(df2, spec)
params = [
    Parameter("S1", "SOURCE", 1.0, 0.0, 100.0),
    Parameter("S2", "SOURCE", 0.1, 0.0, 100.0),
    Parameter("D", "DELIVF", 0.0, -10.0, 10.0),
    Parameter("T", "STRM", 0.1, 0.0, 100.0),
    Parameter("Rv", "RESV", 0.5, 0.0, 100.0),
]
fit = model2.fit(params)
print("recovery fit: SSE =", fit["sse"])
print("recovery fit: beta =", np.round(fit["beta"], 8), " (true:", BETA, ")")
assert fit["sse"] < 1e-16
assert np.max(np.abs(fit["beta"] - BETA)) < 1e-5

print("ALL HAND CHECKS PASSED")
