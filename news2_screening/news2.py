"""news2.py -- Modern Python 3 port of the Global NEWS 2 phosphorus model.

Port of Mayorga's Global NEWS 2 model, GNE implementation
(source: globalnews2gne/globalnews.py, gnecode/gncfg.py, vars.cfg,
constants.cfg; NEWS2_GNEManual.pdf, May 2009). GPL-2, as the source.

Scope of this port: the phosphorus forms only.

* DIP -- dissolved inorganic phosphorus (dissolved sub-model)
* DOP -- dissolved organic phosphorus (dissolved sub-model)
* PP  -- particulate phosphorus (particulate sub-model, via TSS)
* TP  -- DIP + DOP + PP (yields and loads)

The legacy GNE plumbing (vars.cfg / constants.cfg parsing, the IN/OUT
global dictionaries, dbf.py) is replaced by plain CSV readers
(pandas). Every equation and every constant is kept faithful to the
source; source file names are cited in comments next to each block.

Units (as in the source, vars.cfg and Mayorga et al. 2010):
  yields  kg km-2 yr-1
  loads   Mg yr-1 (ton/yr)
  runoff input in mm yr-1, converted to m yr-1 inside the model
  basin area in km2

Usage:
    python news2.py --input-dir inputs/preproctbls --output out.csv
or import and call ``run(input_dir)`` / ``compute_phosphorus(df)``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Calibrated constants -- source: constants.cfg, sections [CAL.DIP] / [CAL.DOP]
# (constants.cfg: "Calibrated and default constants for dissolved sub-models";
# see Mayorga et al. 2010, Supplementary Table A2.)
#   a, b : watershed export function f_F(Rnat) parameters
#   c    : point-source export fraction FEpnt_F (dummy/integration constant
#          for DIP, a calibrated fraction for DOP -- constants.cfg header)
#   e    : watershed export fraction multiplier, FEws_F = e * f_F(Rnat)
#   EC   : diffuse export coefficient (kg km-2 yr-1)
# ---------------------------------------------------------------------------
CAL = {
    "DIP": {"a": 0.85, "b": 2.0, "c": 1.0, "e": 0.29, "EC": 26.0},   # constants.cfg [CAL.DIP]
    "DOP": {"a": 0.95, "c": 0.010, "e": 0.010, "EC": 15.0},          # constants.cfg [CAL.DOP]
}

# Minimum runoff threshold, 3 mm/yr -- source: globalnews.py, R_MIN = 0.003
# (defined above model_dissolved()).
R_MIN = 0.003

# Maximum pre-dam TSS yield, ton km-2 yr-1 -- source: globalnews.py,
# Yld_TSS_pred_MAX = 5000 (above model_particulate()).
YLD_TSS_PRED_MAX = 5000.0

# Lithology chemical-weathering coefficients by dominant lithology class
# (Amiotte-Suchet et al. 2003) -- source: globalnews.py, Calc_Yld_TSS_pred().
LITHO_CHEM_COEF = {1: 0.0, 2: 1.0379, 3: 0.9037, 4: 0.3125,
                   5: 0.0, 6: -2.0678, 7: 0.6694, 8: 0.0}

# Reservoir-retention caps -- source: globalnews.py, Calc_FEriv():
# D_DIP is capped at 0.85 (DIN at 0.965 does not apply to phosphorus).
D_DIP_MAX = 0.85

# Input tables and columns -- source: vars.cfg, sections [IN.TBLS] / [IN.VARS].
# Keys are the model's internal names; values are (file, column) in the
# model's own basin-aggregated input format.
INPUT_COLUMNS = {
    "BasinID":    ("STN30v6ngNEWS.csv", "BASINID"),
    "A":          ("STN30v6ngNEWS.csv", "area"),
    "agric":      ("c00_diffsrc_landuse_ng.csv", "c0agrpct"),
    "wetlndrice": ("c00_diffsrc_landuse_ng.csv", "c0wrpct"),
    "marggrass":  ("c00_diffsrc_landuse_ng.csv", "c0mrgrpct"),
    "WSdif_fe_P": ("c00_diffsrc_landuse_ng.csv", "c0TPfe"),
    "WSdif_ma_P": ("c00_diffsrc_landuse_ng.csv", "c0TPma"),
    "WSdif_ex_P": ("c00_diffsrc_landuse_ng.csv", "c0TPex"),
    "Rnat":       ("c00_hydro_pntsrc_other_ng.csv", "c0R"),
    "frnrprecp":  ("c00_hydro_pntsrc_other_ng.csv", "c0frnrprecp"),
    "RSpntExc_P": ("c00_hydro_pntsrc_other_ng.csv", "c0TPsewhum"),
    "RSpntDet_P": ("c00_hydro_pntsrc_other_ng.csv", "c0TPsewdet"),
    "FQrem":      ("c00_disch_fqrem_reservret_ng.csv", "c0FQremSclGEMS"),
    "D_DIP":      ("c00_disch_fqrem_reservret_ng.csv", "c0Ddip0to1"),
    "D_TSS":      ("c00_disch_fqrem_reservret_ng.csv", "c0Dsed0to1"),
    "frnrslope":  ("LithologySlope_ng.csv", "FournierSlope"),
    "LiClass":    ("LithologySlope_ng.csv", "LithClass"),
}


def load_inputs(input_dir: str | Path) -> pd.DataFrame:
    """Read the NEWS2 basin-aggregated input tables (model's own format).

    Replaces the GNE loader (globalnews.py load_var_arrays() /
    gnecode/gncore.py) with pandas CSV reads. Rows are aligned on
    BasinID, in the basin order of the basins table (vars.cfg: the
    "basins" table, STN30v6ngNEWS.csv).
    """
    input_dir = Path(input_dir)
    frames: dict[str, pd.DataFrame] = {}
    out: pd.DataFrame | None = None
    for var, (fname, col) in INPUT_COLUMNS.items():
        if fname not in frames:
            frames[fname] = pd.read_csv(input_dir / fname)
        tbl = frames[fname]
        idcol = "BASINID" if "BASINID" in tbl.columns else "basinid"
        if var == "BasinID":
            part = tbl[[idcol]].rename(columns={idcol: "BasinID"})
        else:
            part = tbl[[idcol, col]].rename(columns={idcol: "BasinID", col: var})
        out = part if out is None else out.merge(part, on="BasinID", how="left")
    assert out is not None
    return out


# ---------------------------------------------------------------------------
# Dissolved sub-model -- source: globalnews.py, model_dissolved() and its
# helpers Calc_PointSources(), Calc_WSExportFunc(), Calc_WSdifExplicit(),
# Calc_RSdifExportCoeff(), Calc_FEriv().
# ---------------------------------------------------------------------------

def _dissolved(df: pd.DataFrame, form: str) -> dict[str, np.ndarray]:
    cal = CAL[form]
    A = df["A"].to_numpy(float)
    Rnat = df["Rnat"].to_numpy(float) / 1000.0  # mm/yr -> m/yr (model_dissolved)

    # --- Point (sewage) sources: Calc_PointSources() -----------------------
    # FEpnt_F = c for all non-DIN forms; RSpnt_E = human excreta + detergent P.
    FEpnt = cal["c"]
    RSpnt_E = df["RSpntExc_P"].to_numpy(float) + df["RSpntDet_P"].to_numpy(float)
    RSpnt = FEpnt * RSpnt_E

    # --- Watershed export function: Calc_WSExportFunc() --------------------
    # FEws_F = e_F * f_F(Rnat); f_F is set to 0 in arid basins (Rnat < R_MIN).
    if form == "DIP":
        # f = 1 / (1 + (Rnat/a)^-b); near-zero runoff is floored first so the
        # negative exponent cannot blow up (source comment, Calc_WSExportFunc).
        Rgt0 = Rnat.copy()
        Rgt0[Rnat < 0.00001] = 0.0000001
        f = 1.0 / (1.0 + (Rgt0 / cal["a"]) ** -cal["b"])
    else:  # DOP (and DIN/DON/DOC in the source): f = Rnat^a
        f = Rnat ** cal["a"]
    f = f.copy()
    f[Rnat < R_MIN] = 0.0
    FEws = cal["e"] * f
    FEws_nat = FEws.copy()  # only DIN differentiates FEws_nat (humid tropics)

    # --- Diffuse sources, explicit: Calc_WSdifExplicit() -------------------
    # WSdif_ant_E = fertilizer + manure - crop export, floored at 0.
    WSdif_ant = (df["WSdif_fe_P"].to_numpy(float)
                 + df["WSdif_ma_P"].to_numpy(float)
                 - df["WSdif_ex_P"].to_numpy(float))
    WSdif_ant = np.where(WSdif_ant < 0.0, 0.0, WSdif_ant)
    WSdif_nat = np.zeros_like(WSdif_ant)  # natural explicit diffuse P is 0
    RSdif_expl_ant = FEws * WSdif_ant
    RSdif_expl_nat = FEws_nat * WSdif_nat

    # --- Diffuse sources, export coefficient: Calc_RSdifExportCoeff() ------
    # RSdif_ec_F = f_F * EC_F, split agricultural / natural by agric fraction.
    RSdif_ec = f * cal["EC"]
    Agfr = df["agric"].to_numpy(float) / 100.0
    RSdif_ec_ant = Agfr * RSdif_ec
    RSdif_ec_nat = (1.0 - Agfr) * RSdif_ec

    RSdif_ant = RSdif_expl_ant + RSdif_ec_ant
    RSdif_nat = RSdif_expl_nat + RSdif_ec_nat

    # --- River export fraction: Calc_FEriv() -------------------------------
    # FEriv = (1 - L)(1 - D)(1 - FQrem). L_F = 0 for all non-DIN forms.
    # D_F: DIP uses the basin reservoir retention D_DIP capped at 0.85;
    # DOP (like DON/DOC) has D_F = 0.
    if form == "DIP":
        D = np.minimum(df["D_DIP"].to_numpy(float), D_DIP_MAX)
    else:
        D = 0.0
    FEriv = (1.0 - 0.0) * (1.0 - D) * (1.0 - df["FQrem"].to_numpy(float))

    # --- Yield and load: model_dissolved() ---------------------------------
    Yld = FEriv * (RSpnt + (RSdif_ant + RSdif_nat))     # kg km-2 yr-1
    Ld = Yld * A / 1000.0                                # Mg yr-1
    return {"Yld": Yld, "Ld": Ld, "FEws": FEws, "FEriv": np.broadcast_to(FEriv, A.shape).copy()
            if np.isscalar(FEriv) else FEriv, "RSpnt": np.full_like(A, FEpnt) * RSpnt_E,
            "RSdif_ant": RSdif_ant, "RSdif_nat": RSdif_nat}


# ---------------------------------------------------------------------------
# Particulate sub-model (PP via TSS) -- source: globalnews.py,
# model_particulate() and Calc_Yld_TSS_pred() (Beusen et al. 2005 TSS
# regression; Ludwig et al. %POC relation; PP eq. 4 of Mayorga et al. 2010).
# ---------------------------------------------------------------------------

def _tss_yield_pred(df: pd.DataFrame) -> np.ndarray:
    """Pre-dam TSS yield, ton km-2 yr-1: Calc_Yld_TSS_pred()."""
    clip = lambda s, lo, hi: s.to_numpy(float).clip(min=lo, max=hi)  # noqa: E731
    MargGrass = clip(df["marggrass"], 0, 84.6865)      # truncation ranges: source
    WetlnRice = clip(df["wetlndrice"], 0, 24.7206)
    FrnrPrecp = clip(df["frnrprecp"], 0.854985, 11.1682)
    FrnrSlope = clip(df["frnrslope"], 1.1301, 38.4531)
    LiClass = df["LiClass"].to_numpy(int)
    LithoChem = np.array([LITHO_CHEM_COEF[i] for i in LiClass])
    yld = np.exp(0.0301 * MargGrass + 0.1234 * WetlnRice + 0.2933 * FrnrPrecp
                 + 0.1027 * FrnrSlope + LithoChem + 0.5522)
    yld[LiClass == 8] = 0.0  # class 8 = ice / no lithology overlay (source)
    return yld


def _particulate_pp(df: pd.DataFrame) -> dict[str, np.ndarray]:
    A = df["A"].to_numpy(float)
    Rnat = df["Rnat"].to_numpy(float) / 1000.0  # m/yr (model_particulate)
    D_TSS = df["D_TSS"].to_numpy(float)

    Yld_TSS_pred = _tss_yield_pred(df)
    Yld_TSS_pred[Rnat < R_MIN] = 0.0
    Yld_TSS_pred[Yld_TSS_pred > YLD_TSS_PRED_MAX] = YLD_TSS_PRED_MAX
    Yld_TSS = (1.0 - D_TSS) * Yld_TSS_pred
    NoFlux = (Yld_TSS == 0) | (Rnat < R_MIN)
    Ld_TSS = Yld_TSS * A  # ton/yr; zeroed for NoFlux basins below, as in source

    with np.errstate(divide="ignore", invalid="ignore"):
        TSSc_pred = Yld_TSS_pred / Rnat            # mg/L (source unit note)
        logTSSc = np.log10(TSSc_pred)
        # %POC of TSS, Ludwig et al. polynomial in log10(TSSc): model_particulate()
        TSSpc_POC = (-0.160 * logTSSc ** 3 + 2.83 * logTSSc ** 2
                     - 13.6 * logTSSc + 20.3)
    TSSc_pred = np.where(NoFlux, 0.0, TSSc_pred)
    Ld_TSS = np.where(NoFlux, 0.0, Ld_TSS)
    TSSpc_POC = np.where(NoFlux, 0.0, TSSpc_POC)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # POC load under natural (pre-dam) conditions, kg/yr: model_particulate()
        Ld_POC_pred_kgyr = 1000.0 * (TSSpc_POC / 100.0) * Ld_TSS / (1.0 - D_TSS)
        # PP load, eq. 4: Ld_F_pred in ton/yr after the /1000 conversion.
        Ld_pred = np.exp(-3.09614 + 1.001769 * np.log(Ld_POC_pred_kgyr)) / 1000.0
        Ld = (1.0 - D_TSS) * Ld_pred
    Ld = np.where(NoFlux, 0.0, Ld)  # source zeroes Ld_F for invalid cases
    Ld = np.nan_to_num(Ld, nan=0.0, posinf=0.0, neginf=0.0)
    Yld = 1000.0 * Ld / A           # kg km-2 yr-1
    return {"Yld": Yld, "Ld": Ld, "Yld_TSS": Yld_TSS, "TSSpc_POC": TSSpc_POC}


# ---------------------------------------------------------------------------

def compute_phosphorus(df: pd.DataFrame) -> pd.DataFrame:
    """Compute DIP, DOP, PP and TP river export for each basin row of df."""
    out = pd.DataFrame({"BasinID": df["BasinID"].to_numpy(int),
                        "A_km2": df["A"].to_numpy(float)})
    for form in ("DIP", "DOP"):
        r = _dissolved(df, form)
        out[f"{form}_Yld"] = r["Yld"]
        out[f"{form}_Ld"] = r["Ld"]
        out[f"{form}_FEws"] = r["FEws"]
        out[f"{form}_FEriv"] = r["FEriv"]
    pp = _particulate_pp(df)
    out["PP_Yld"] = pp["Yld"]
    out["PP_Ld"] = pp["Ld"]
    out["Yld_TSS"] = pp["Yld_TSS"]
    # TP is the sum of the three phosphorus forms (cf. allformmodel.py, where
    # TPyld = DIPYld + DOPYld + PPYld for the ICEP calculation).
    out["TP_Yld"] = out["DIP_Yld"] + out["DOP_Yld"] + out["PP_Yld"]
    out["TP_Ld"] = out["DIP_Ld"] + out["DOP_Ld"] + out["PP_Ld"]
    return out


def run(input_dir: str | Path) -> pd.DataFrame:
    return compute_phosphorus(load_inputs(input_dir))


def main() -> None:
    ap = argparse.ArgumentParser(description="Global NEWS 2 phosphorus port (DIP, DOP, PP, TP)")
    ap.add_argument("--input-dir", required=True,
                    help="Directory holding the NEWS2 pre-processed input tables")
    ap.add_argument("--output", required=True, help="Output CSV path")
    args = ap.parse_args()
    res = run(args.input_dir)
    res.to_csv(args.output, index=False)
    print(f"Wrote {len(res)} basins to {args.output}")
    print("Global TP load (sum of basin loads): "
          f"{res['TP_Ld'].sum():,.1f} Mg/yr")


if __name__ == "__main__":
    main()
