"""
run_vignette.py -- run the EEMMA demonstration analyses from Kirchner's
vignette script (eemma_vignette_script.r) on the author's demo data
(vignette_data.txt), using the Python port in eemma.py.

Usage:  python run_vignette.py [path/to/vignette_data.txt]
"""

import sys

import pandas as pd

from eemma import eemma

pd.set_option("display.width", 120)


def main(path="vignette_data.txt"):
    dat = pd.read_csv(path, sep="\t")
    precip, soil, gw = dat["precip"], dat["soil"], dat["gw"]
    xylem1, xylem2, stream = dat["xylem1"], dat["xylem2"], dat["stream"]
    X3 = pd.DataFrame({"precip": precip, "soil": soil, "gw": gw})
    X2 = pd.DataFrame({"precip": precip, "soil": soil})
    X1 = pd.DataFrame({"precip": precip})

    cases = [
        ("xylem1 ~ precip+soil+gw (complete)", xylem1, X3, {}),
        ("xylem1 ~ precip+soil, complete=FALSE", xylem1, X2,
         {"complete": False}),
        ("xylem1 ~ precip+soil+gw, complete=FALSE", xylem1, X3,
         {"complete": False}),
        ("xylem1 ~ precip+soil (complete)", xylem1, X2, {}),
        ("xylem1 ~ precip+soil+gw, memory=TRUE", xylem1, X3,
         {"memory": True}),
        ("xylem2 ~ precip+soil+gw, memory=TRUE", xylem2, X3,
         {"memory": True}),
        ("xylem2 ~ precip+soil+gw, memory=FALSE", xylem2, X3, {}),
        ("stream ~ precip+soil+gw (complete)", stream, X3, {}),
        ("stream ~ precip+soil+gw, complete=FALSE", stream, X3,
         {"complete": False}),
        ("stream ~ precip+soil+gw, memory=TRUE", stream, X3,
         {"memory": True}),
        ("stream ~ precip only, complete=FALSE", stream, X1,
         {"complete": False}),
    ]
    for title, y, X, kw in cases:
        res = eemma(y, X, **kw)
        print(f"\n### {title}")
        print(res["coeffs"].to_string(float_format=lambda v: f"{v:.8g}"))
        print(f"RMSE = {res['RMSE']:.6g}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "vignette_data.txt")
