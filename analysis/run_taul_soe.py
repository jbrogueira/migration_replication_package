# run_taul_soe.py
import os
import pandas as pd

from config import τl as TAU_L_BASE, pref_tax, σ, G0, tau_grid as TAU_GRID, SIG_TAG
from model.equilibrium import pi_ge_eq, sweep_taul_soe_pref

OUTPUT = os.path.join(os.getcwd(), "output", "tables", "csv")
os.makedirs(OUTPUT, exist_ok=True)


def main():
    lam_star = float(pi_ge_eq[4])
    df_pref = sweep_taul_soe_pref(lam_star, TAU_GRID, τ10h_pref=pref_tax)

    # --- save full sweeps (optional) ---
    # df_pref.to_csv(
    #     os.path.join(OUTPUT, f"taul_case1_soe_pref_tax_{int(pref_tax*100)}.csv"),
    #     index=False,
    # )

    # exact: λ0 at λ*, τl clears the budget (θ = inf: immigrant indifference)
    import numpy as np
    from model.equilibrium_diff_dest import solve_diff_dest_taul, equilibrium_row
    try:
        p0, tau_l_star, τ0 = solve_diff_dest_taul(custom_τ10h=pref_tax, θ=np.inf)
    except (ValueError, AssertionError) as e:
        # ValueError: revenue minus G0 has no sign change on the τl bracket.
        # AssertionError: the stock solve fails at a bracket endpoint because no
        # allocation exists there (at ϵ = 20 the low-skill immigrant stock reaches the
        # solver bound below τl ≈ 0.10, with revenue above G0 at every solvable τl).
        print("DIAGNOSTIC FLAG: no budget-clearing progressivity τl in the perfect-substitutes "
              f"limit ({type(e).__name__}: {e}); the equilibrium CSV row is NOT an equilibrium")
        df_pref_eq = df_pref[df_pref["TR"] > G0].head(1).assign(budget_crossing=False)
    else:
        row = equilibrium_row(p0, lam_star, τ_override=τ0)
        assert abs(row["TR"] - G0) < 1e-8, "budget clearing failed in the τl exercise (θ = inf)"
        row["tau_l"] = tau_l_star
        df_pref_eq = pd.DataFrame([row]).assign(budget_crossing=True)
    df_pref_eq.to_csv(
        os.path.join(
            OUTPUT,
            f"case1_soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{TAU_L_BASE}_equilibrium.csv",
        ),
        index=False,
    )

    return df_pref


if __name__ == "__main__":
    main()
