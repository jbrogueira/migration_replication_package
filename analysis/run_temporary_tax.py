# analysis/run_temporary_tax.py
"""
Temporary tax robustness check (Referee 2, comment 2).

The baseline preferential tax (pref_tax=0.13) assumes the entire stock of
high-skill immigrants benefits permanently. Real policies are temporary
(e.g., 5 years for new residents). We rescale the tax cut by the share of
exempt years to get an effective lifetime rate.

Temporary: tau10h_eff = (T_exempt/T_total) * 0.13 + (1 - T_exempt/T_total) * tau10h_base,
with tau10h_base the high-skill immigrants' average tax rate at the calibrated baseline.
"""

import os
import numpy as np
from analysis.exact_refine import refine_crossing
from config import σ, τl, G0, pref_tax, SIG_TAG, DD_TAG
from model.equilibrium import pi_ge_eq, ge_pref_residual_unilateral, sweep_equilibrium_2country, tau10h_star
from analysis.run_pref_tax import main as run_pref_soe
from analysis.build_equilibrium_tables import (build_equilibrium_tables, build_equilibrium_tables1,
                                               build_equilibrium_tables_diff_dest)
from config import λ0_vals, bounds_lambda, bounds_lambda_8


# Baseline tax on high-skill immigrants (from the calibrated equilibrium)
TAU10H_BASELINE = float(tau10h_star)

def effective_pref_tax(tau_pref, tau_baseline, T_exempt, T_total):
    """Compute effective lifetime preferential tax rate."""
    share = T_exempt / T_total
    return share * tau_pref + (1 - share) * tau_baseline


def run_pref_large_temp(pref_tax_temp):
    """Two-country unilateral preferential tax with the temporary (effective) rate;
    λ0 clears country 0's budget, λ1 clears country 1's."""
    λ1_star = pi_ge_eq[4]
    x0_one = pi_ge_eq[:4]
    x0 = np.concatenate((x0_one, x0_one, [λ1_star]))

    err_ge = lambda params, out, custom_λ0, custom_λ1, *_, **__: ge_pref_residual_unilateral(
        params, out, custom_λ0, None, pref_tax_temp
    )
    df = sweep_equilibrium_2country(
        err_func=err_ge, λ0_vals=λ0_vals, λ1_fixed=None, bounds=bounds_lambda_8,
        x0=x0, free_lambda1=True, tau10h_override=pref_tax_temp)

    os.makedirs("output/tables/csv", exist_ok=True)
    cols = ["p00l", "p10l", "p00h", "p10h", "p11l", "p01l", "p11h", "p01h", "λ1"]
    warm = df[df["TR0"] >= G0].head(1)[cols].to_numpy()[0]
    eq = refine_crossing(
        df, "λ0", "TR0", G0,
        lambda l0: sweep_equilibrium_2country(err_ge, [l0], None, bounds_lambda_8, x0=warm,
                                              free_lambda1=True, tau10h_override=pref_tax_temp))
    tag = int(pref_tax_temp * 100)
    eq.to_csv(
        f"output/tables/csv/two_country_pref_{tag}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        index=False
    )
    return df


def main():
    # --- 5/25 years scenario (referee's example) ---
    tau_temp_5 = effective_pref_tax(pref_tax, TAU10H_BASELINE, T_exempt=5, T_total=25)
    print(f"Temporary tax (5/25 years): effective rate = {tau_temp_5:.3f}")

    # --- 10/25 years scenario (Portugal NHR = 10 years) ---
    tau_temp_10 = effective_pref_tax(pref_tax, TAU10H_BASELINE, T_exempt=10, T_total=25)
    print(f"Temporary tax (10/25 years): effective rate = {tau_temp_10:.3f}")

    # Small country with differentiated destinations (active θ): exact budget clearing
    import pandas as pd
    from model.equilibrium_diff_dest import solve_diff_dest, equilibrium_row
    for τ_t in (tau_temp_5, tau_temp_10):
        x = solve_diff_dest(τ_t)
        path = f"output/tables/csv/diff_dest_pref_tax_{int(τ_t * 100)}_{DD_TAG}_taul_{τl}_equilibrium.csv"
        pd.DataFrame([equilibrium_row(x[:4], x[4], τ10h_fixed=τ_t)]).to_csv(path, index=False)
        build_equilibrium_tables_diff_dest(pref_eq_path=path)

    # Run SOE (small country) for both scenarios
    print("\n=== SOE: 5/25 years ===")
    run_pref_soe(pref_tax=tau_temp_5)

    print("\n=== SOE: 10/25 years ===")
    run_pref_soe(pref_tax=tau_temp_10)

    # Run large country for both scenarios
    print("\n=== Large country: 5/25 years ===")
    run_pref_large_temp(tau_temp_5)

    print("\n=== Large country: 10/25 years ===")
    run_pref_large_temp(tau_temp_10)

    # Build comparison tables
    tag_5 = int(tau_temp_5 * 100)
    tag_10 = int(tau_temp_10 * 100)

    # SOE tables
    for tag, tau_t in [(tag_5, tau_temp_5), (tag_10, tau_temp_10)]:
        pref_path = f"output/tables/csv/soe_pref_tax_{tag}_{SIG_TAG}_taul_{τl}_equilibrium.csv"
        if os.path.exists(pref_path):
            build_equilibrium_tables(pref_eq_path=pref_path)

    # Large country tables
    for tag, tau_t in [(tag_5, tau_temp_5), (tag_10, tau_temp_10)]:
        pref_path = f"output/tables/csv/two_country_pref_{tag}_{SIG_TAG}_taul_{τl}_equilibrium.csv"
        if os.path.exists(pref_path):
            build_equilibrium_tables(pref_eq_path=pref_path)
            build_equilibrium_tables1(pref_eq_path=pref_path)

    print("\nDone. Results in output/tables/csv/ and output/tables/tex/")


if __name__ == "__main__":
    main()
