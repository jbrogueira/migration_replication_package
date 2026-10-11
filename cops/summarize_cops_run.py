"""Endpoint and welfare summary of a closed-budget preferential-tax run (tau_P = 0.13),
for the appendix tables tab:cops-transition and tab:cops-cev.

Usage: python summarize_cops_run.py --dir recalib_v6 --prefix v6
reads {prefix}_{cc}_taup_0p13.mat, pt_stocks_{prefix}_{cc}_0p13.mat and
cev_paths_{prefix}_{cc}_taup13.mat for cc in pt, es, de. End of the transition = slot 199
(index 198; slot 200 is the pinned terminal slice); announcement = slot 2 (index 1).
CEV rows: employed workers of the adopting country (row = country index).
"""
import argparse
from pathlib import Path
import scipy.io

CTRY = {"pt": 16, "es": 8, "de": 5}
END, ANN = 198, 1


def load(path):
    return scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)


def summary(d, prefix, cc):
    st = load(d / f"pt_stocks_{prefix}_{cc}_0p13.mat")["St"]
    dg = load(d / f"{prefix}_{cc}_taup_0p13.mat")["info"].dg
    cv = load(d / f"cev_paths_{prefix}_{cc}_taup13.mat")["out"]
    ic = CTRY[cc]
    lam = dg.lambda_hat_path[1:-1]
    def c(name):
        v = getattr(cv, f"cev_{name}_{ic}")
        return v[ANN], v[END]
    return {
        "h_stock": 100 * (st.h_ratio[END] - 1), "pool": 100 * st.pref_share_of_h[END],
        "l_stock": 100 * (st.l_ratio[END] - 1), "w_h": 100 * (dg.What_h_path[END] - 1),
        "w_l": 100 * (dg.What_l_path[END] - 1), "tau_h0": dg.tau_h_t_path[ANN],
        "tau_h1": dg.tau_h_t_path[END], "lam_dev": 100 * max(abs(lam.min() - 1), abs(lam.max() - 1)),
        "budget_resid": abs(dg.budget_resid).max(),
        "pref_eu15": c("eu15_h_pref"), "pref_nms8": c("nms8_h_pref"), "pref_oth": c("oth_h_pref"),
        "std_eu15": c("eu15_h_std"), "low_eu15": c("eu15_l"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--prefix", required=True)
    a = ap.parse_args()
    d = Path(__file__).parent / a.dir
    out = {cc: summary(d, a.prefix, cc) for cc in CTRY}
    for k in ("h_stock", "pool", "l_stock", "w_h", "w_l", "tau_h0", "tau_h1", "lam_dev", "budget_resid"):
        print(f"{k:13s}" + "".join(f"{out[cc][k]:>12.4g}" for cc in CTRY))
    for k in ("pref_eu15", "pref_nms8", "pref_oth", "std_eu15", "low_eu15"):
        print(f"{k:13s}" + "".join(f"{out[cc][k][0]:>8.3f} {out[cc][k][1]:>7.3f}" for cc in CTRY))


if __name__ == "__main__":
    main()
