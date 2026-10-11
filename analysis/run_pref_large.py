import os
import numpy as np
from model.equilibrium import pi_ge_eq, sweep_equilibrium_2country, ge_pref_residual_unilateral
from analysis.exact_refine import refine_crossing
from config import λ0_vals, bounds_lambda_8, pref_tax, G0, σ, τl, SIG_TAG

def main():
    λ1_star = pi_ge_eq[4]
    x0_one = pi_ge_eq[:4]
    x0 = np.concatenate((x0_one, x0_one))

    # unilateral adoption: country 0 clears its budget with λ0; country 1 (no
    # preferential tax) clears its own with λ1. The sweep is the figure's; the
    # reported equilibrium is its exact budget-clearing point.
    err_ge = lambda params, out, custom_λ0, custom_λ1, *_, **__: ge_pref_residual_unilateral(
        params, out, custom_λ0, None, pref_tax
    )
    df_ge = sweep_equilibrium_2country(
        err_func=err_ge, λ0_vals=λ0_vals, λ1_fixed=None, bounds=bounds_lambda_8,
        x0=np.concatenate([x0, [λ1_star]]), free_lambda1=True, tau10h_override=pref_tax)
    df_ge.to_csv(f"output/tables/csv/two_country_pref_ge_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}.csv",
                 index=False)
    cols = ["p00l", "p10l", "p00h", "p10h", "p11l", "p01l", "p11h", "p01h", "λ1"]
    warm = df_ge[df_ge["TR0"] >= G0].head(1)[cols].to_numpy()[0]
    eq = refine_crossing(
        df_ge, "λ0", "TR0", G0,
        lambda l0: sweep_equilibrium_2country(err_ge, [l0], None, bounds_lambda_8, x0=warm,
                                              free_lambda1=True, tau10h_override=pref_tax))
    eq.to_csv(
        f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        index=False
    )
    return df_ge

if __name__ == "__main__":
    main()
