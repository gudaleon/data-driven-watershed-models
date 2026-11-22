"""Verify the NEWS2 phosphorus port against the reference output bundled
with the source repository (output/c00_NEWS2Output.csv, produced by the
original GNE code from inputs/preproctbls).

Compares DIP, DOP and PP yields (kg km-2 yr-1) and loads (Mg yr-1), and
TP as the sum of the three forms, basin by basin.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from news2 import run

SRC = Path.home() / "workspace/fei_lake_ontario/codebase/globalnews2gne"
res = run(SRC / "inputs/preproctbls").set_index("BasinID")
ref = pd.read_csv(SRC / "output/c00_NEWS2Output.csv").set_index("basinid")
ref = ref.loc[res.index]

pairs = [("DIP_Yld", "c0Yld_DIP"), ("DIP_Ld", "c0Ld_DIP"),
         ("DOP_Yld", "c0Yld_DOP"), ("DOP_Ld", "c0Ld_DOP"),
         ("PP_Yld", "c0Yld_PP"), ("PP_Ld", "c0Ld_PP")]
ref["TP_Yld_ref"] = ref["c0Yld_DIP"] + ref["c0Yld_DOP"] + ref["c0Yld_PP"]
ref["TP_Ld_ref"] = ref["c0Ld_DIP"] + ref["c0Ld_DOP"] + ref["c0Ld_PP"]
pairs += [("TP_Yld", "TP_Yld_ref"), ("TP_Ld", "TP_Ld_ref")]

print(f"basins compared: {len(res)}")
for mine, theirs in pairs:
    a, b = res[mine].to_numpy(float), ref[theirs].to_numpy(float)
    d = np.abs(a - b)
    rel = np.where(np.abs(b) > 0, d / np.abs(b), d)
    print(f"{mine:8s} max_abs_diff={d.max():.6g}  basins_exact={int((d == 0).sum())}"
          f"  basins_rel<=1e-9={int((rel <= 1e-9).sum())}  max_rel={rel.max():.3g}"
          f"  sum_mine={a.sum():,.2f}  sum_ref={b.sum():,.2f}")
