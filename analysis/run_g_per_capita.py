# analysis/run_g_per_capita.py
# Small country with the preferential tax, government spending proportional to the
# resident population: the budget clears at TR = G0 · N, N = p00l + p10l + p00h + p10h
# (N = 1 at the baseline). λ0 clears the budget; everything else as in run_diff_dest.
import argparse
import os

import numpy as np
import pandas as pd

from model.equilibrium import pi_ge_eq, tau_star
from model.equilibrium_diff_dest import (
    err_diff_dest, equilibrium_row, required_spending, solve_diff_dest,
    solve_diff_dest_damped, sweep_diff_dest,
)
from analysis.build_equilibrium_tables import build_equilibrium_tables_diff_dest
from viz.plotting import plot_g_per_capita_laffer
from config import pref_tax as pref_tax_default, λ0_vals, G0, DD_TAG, τl

CSV_DIR = "output/tables/csv"
FIG_DIR = "output/figures"


def _solve_exact(pref_tax: float, g_per_capita: bool) -> np.ndarray:
    """SLSQP + polish, cross-validated against the damped/brentq solve."""
    x_p = solve_diff_dest(pref_tax, g_per_capita=g_per_capita)
    x_d = None
    for damping in (0.5, 0.15, 0.05):
        try:
            x_d = solve_diff_dest_damped(pref_tax, damping=damping, g_per_capita=g_per_capita)
            break
        except RuntimeError:
            continue
    assert x_d is not None, "damped solve failed at all damping values"
    assert np.max(np.abs(x_p - x_d)) < 1e-6, "SLSQP and damped solves disagree"
    return x_p


def main(pref_tax: float = pref_tax_default):
    os.makedirs(CSV_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    tag = f"pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}"

    # (1) zero-policy gate: at N = 1 the per-resident rule reproduces the baseline
    x_zero = solve_diff_dest(None, g_per_capita=True)
    assert np.max(np.abs(x_zero - pi_ge_eq)) < 1e-6, "zero-policy gate failed"
    row_b = equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4], fun=err_diff_dest(pi_ge_eq, out=1))
    base_path = f"{CSV_DIR}/diff_dest_baseline_{DD_TAG}_taul_{τl}_equilibrium.csv"
    pd.DataFrame([row_b]).to_csv(base_path, index=False)

    # (2) policy equilibria at τ_P: fixed G0 (as in the main text) and G0 · N
    rows = {}
    for name, gpc in (("fixedG", False), ("gpc", True)):
        x = _solve_exact(pref_tax, gpc)
        row = equilibrium_row(x[:4], x[4], τ10h_fixed=pref_tax,
                              fun=err_diff_dest(x, out=1, custom_τ10h=pref_tax, g_per_capita=gpc))
        row["N"] = float(np.sum(x[:4]))
        row["G"] = required_spending(x, gpc)
        assert abs(row["TR"] - row["G"]) < 1e-8, f"budget does not clear ({name})"
        rows[name] = row
        path = f"{CSV_DIR}/diff_dest_{name}_{tag}_equilibrium.csv"
        pd.DataFrame([row]).to_csv(path, index=False)

    # (3) table: baseline vs the G0 · N equilibrium
    _, _, tex_path = build_equilibrium_tables_diff_dest(
        baseline_eq_path=base_path,
        pref_eq_path=f"{CSV_DIR}/diff_dest_gpc_{tag}_equilibrium.csv",
        caption=("{Unilateral preferential tax in the small country with differentiated "
                 "destinations, government spending proportional to the resident population: "
                 r"baseline~(B) vs.\ policy~(P). Migration shares, wages, tax rates, "
                 "consumption, and CEV welfare effects.}"),
        label="tab:g-per-capita",
    )

    # (4) Laffer curves: revenue and required spending along the λ0 sweeps
    df_base = sweep_diff_dest(λ0_vals)
    df_pref = sweep_diff_dest(λ0_vals, custom_τ10h=pref_tax)
    plot_g_per_capita_laffer(
        df_base, df_pref, tau_star, G0,
        eq_fixed=(rows["fixedG"]["λ0"], rows["fixedG"]["N"]),
        eq_gpc=(rows["gpc"]["λ0"], rows["gpc"]["N"]),
        filename=os.path.join(FIG_DIR, f"g_per_capita_laffer_{tag}.png"),
    )

    return row_b, rows, tex_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pref-tax", type=float, default=pref_tax_default)
    args = ap.parse_args()
    row_b, rows, tex_path = main(args.pref_tax)
    cols = ["λ0", "N", "G", "TR", "p00l", "p10l", "p00h", "p10h",
            "w00l", "w10l", "w00h", "w10h", "tau00l", "tau10l", "tau00h", "tau10h",
            "inc00l", "inc10l", "inc00h", "inc10h"]
    print(pd.DataFrame({"baseline": row_b, **rows}).T.reindex(columns=cols).to_string(
        float_format=lambda v: f"{v:.5f}"))
    print(tex_path)
