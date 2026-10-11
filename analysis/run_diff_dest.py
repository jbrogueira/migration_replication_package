# analysis/run_diff_dest.py
# Differentiated-destinations small country: exact budget-clearing equilibria,
# Laffer sweeps, and the cross-check CSV (docs/model_reformulation_code_plan.md, Step 3).
import os

import numpy as np
import pandas as pd

from model.equilibrium import pi_ge_eq
from model.equilibrium_diff_dest import (
    err_diff_dest, err_fixed_foreign_netw, equilibrium_row, solve_diff_dest,
    solve_diff_dest_damped, sweep_diff_dest, solve_fixed_foreign_netw,
)
from model.parameters import κl, κh
from config import λ0_vals, pref_tax, G0, σ, τl, SIG_TAG, DD_TAG

CSV_DIR = "output/tables/csv"
TWO_COUNTRY_CSV = f"{CSV_DIR}/two_country_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv"


def _native_cev(row_p: dict, row_b: dict, c01_p: dict = None) -> dict:
    """
    Native CEV by skill: ((c00)^σ + (c01/κ)^σ)^(1/σ).
    The move option c01 is the consumption abroad: the foreign net wage (unchanged
    for the small country) plus country-0 land income, which moves with the policy;
    c01_p overrides it (two-country row, where the other country adjusts).
    """
    from analysis.build_equilibrium_tables import native_abroad_income
    out = {}
    for a, κ_a in (('l', κl), ('h', κh)):
        c01_b = row_b[f"inc10{a}"]
        c01 = c01_p[a] if c01_p is not None else native_abroad_income(row_b, row_p, a)
        num_b = (row_b[f"inc00{a}"]**σ + (c01_b / κ_a)**σ) ** (1.0 / σ)
        num_p = (row_p[f"inc00{a}"]**σ + (c01 / κ_a)**σ) ** (1.0 / σ)
        out[a] = num_p / num_b - 1.0
    return out


def _crosscheck_row(scenario: str, row_p: dict, row_b: dict, c01_p: dict = None) -> dict:
    cev = _native_cev(row_p, row_b, c01_p)
    return {
        "scenario": scenario,
        "d_pi10l": row_p["p10l"] - row_b["p10l"],
        "d_pi10h": row_p["p10h"] - row_b["p10h"],
        "pct_dw00h": row_p["w00h"] / row_b["w00h"] - 1.0,
        "pct_dc00h": row_p["inc00h"] / row_b["inc00h"] - 1.0,
        "cev_native_l": cev['l'],
        "cev_native_h": cev['h'],
        "TR": row_p["TR"],
    }


def main(pref_tax=pref_tax):
    os.makedirs(CSV_DIR, exist_ok=True)

    # (1) zero-policy consistency solve; exact baseline row
    # 1e-6 accommodates pi_ge_eq's own SLSQP residual at high σ; the σ=2 test
    # gate (tests/test_diff_dest.py) keeps 1e-8
    x_zero = solve_diff_dest(None)
    assert np.max(np.abs(x_zero - pi_ge_eq)) < 1e-6, "zero-policy gate failed"
    row_b = equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4], fun=err_diff_dest(pi_ge_eq, out=1))
    pd.DataFrame([row_b]).to_csv(
        f"{CSV_DIR}/diff_dest_baseline_{DD_TAG}_taul_{τl}_equilibrium.csv", index=False)

    # (2) policy equilibrium at τ_P, exact budget clearing; damped cross-validation
    # (smaller damping at high σ, where the stock map responds more strongly)
    x_p = solve_diff_dest(pref_tax)
    x_d = None
    for damping in (0.5, 0.15, 0.05):
        try:
            x_d = solve_diff_dest_damped(pref_tax, damping=damping)
            break
        except RuntimeError:
            continue
    assert x_d is not None, "damped solve failed at all damping values"
    assert np.max(np.abs(x_p - x_d)) < 1e-6, "SLSQP and damped solves disagree"
    row_p = equilibrium_row(x_p[:4], x_p[4], τ10h_fixed=pref_tax,
                            fun=err_diff_dest(x_p, out=1, custom_τ10h=pref_tax))
    pd.DataFrame([row_p]).to_csv(
        f"{CSV_DIR}/diff_dest_pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}_equilibrium.csv",
        index=False)

    # (3) Laffer sweeps with and without the preferential tax
    df_sweep0 = sweep_diff_dest(λ0_vals)
    df_sweep0.to_csv(f"{CSV_DIR}/diff_dest_{DD_TAG}_taul_{τl}.csv", index=False)
    df_sweep1 = sweep_diff_dest(λ0_vals, custom_τ10h=pref_tax)
    df_sweep1.to_csv(
        f"{CSV_DIR}/diff_dest_pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}.csv",
        index=False)

    # (4) cross-check CSV: differentiated destinations vs fixed-foreign-net-wages
    #     variant vs two-country PE
    rows = [_crosscheck_row("diff_dest", row_p, row_b)]

    x_f = solve_fixed_foreign_netw(pref_tax)
    row_f = equilibrium_row(x_f[:4], x_f[4], τ10h_fixed=pref_tax,
                            fun=err_fixed_foreign_netw(x_f, out=1, custom_τ10h=pref_tax))
    rows.append(_crosscheck_row("fixed_foreign_netw", row_f, row_b))

    if os.path.exists(TWO_COUNTRY_CSV):
        tc = pd.read_csv(TWO_COUNTRY_CSV).iloc[0].to_dict()
        tc["TR"] = tc["TR0"]
        rows.append(_crosscheck_row("two_country_pe", tc, row_b,
                                    c01_p={'l': tc["inc01l"], 'h': tc["inc01h"]}))
    else:
        print(f"two-country PE CSV not found, cross-check row skipped: {TWO_COUNTRY_CSV}")

    df_cross = pd.DataFrame(rows)
    df_cross.to_csv(
        f"{CSV_DIR}/diff_dest_crosscheck_pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}.csv",
        index=False)

    return row_b, row_p, df_sweep0, df_sweep1, df_cross


if __name__ == "__main__":
    main()
