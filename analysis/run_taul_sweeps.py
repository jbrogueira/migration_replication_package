# analysis/run_taul_sweeps.py
import os
from analysis.exact_refine import refine_crossing
from config import σ, τl as TAU_L_BASE, pref_tax, G0, tau_grid as TAU_GRID, SIG_TAG
from model.equilibrium import (
    pi_ge_eq,
    sweep_taul_large_case2,
    sweep_taul_large_case3,
)

OUTPUT = os.path.join(os.getcwd(), "output", "tables", "csv")
os.makedirs(OUTPUT, exist_ok=True)


def main():
    lam_star = float(pi_ge_eq[4])

    # Case 2: Country 0 pref (τ10h), Country 1 none
    df_case2 = sweep_taul_large_case2(
        lambda_star=lam_star,
        tau_grid=TAU_GRID,
        τ10h_pref=pref_tax,
    )

    # Case 3: Bilateral pref (τ10h = τ01h)
    df_case3 = sweep_taul_large_case3(
        lambda_star=lam_star,
        tau_grid=TAU_GRID,
        τ_pref=pref_tax,
    )

    # --- save full sweeps (optional) ---
    # path2 = os.path.join(
    #     OUTPUT, f"two_country_taul_case2_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{TAU_L_BASE}.csv"
    # )
    # path3 = os.path.join(
    #     OUTPUT, f"two_country_taul_case3_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{TAU_L_BASE}.csv"
    # )
    # df_case2.to_csv(path2, index=False)
    # df_case3.to_csv(path3, index=False)

    # === Select equilibrium rows: first with TR0 > G0 (error if none) ===
    if "TR0" not in df_case2.columns or "TR0" not in df_case3.columns:
        raise ValueError("Expected column 'TR0' in sweep outputs to select by revenue.")

    # exact budget clearing in country 0 by the progressivity parameter
    cols = ["p00l", "p10l", "p00h", "p10h", "p11l", "p01l", "p11h", "p01h"]
    warm2 = df_case2[df_case2["TR0"] > G0].head(1)[cols + ["tauL1"]].to_numpy()[0]
    df_case2_eq = refine_crossing(
        df_case2, "tauL0", "TR0", G0,
        lambda t: sweep_taul_large_case2(lam_star, [t], pref_tax, x0=warm2), select="gt")
    warm3 = df_case3[df_case3["TR0"] > G0].head(1)[cols].to_numpy()[0]
    df_case3_eq = refine_crossing(
        df_case3, "tauL0", "TR0", G0,
        lambda t: sweep_taul_large_case3(lam_star, [t], pref_tax, x0=warm3), select="gt")

    # Save equilibria (1-row)
    path2_eq = os.path.join(
        OUTPUT,
        f"two_country_taul_case2_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{TAU_L_BASE}_equilibrium.csv",
    )
    path3_eq = os.path.join(
        OUTPUT,
        f"two_country_taul_case3_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{TAU_L_BASE}_equilibrium.csv",
    )
    df_case2_eq.to_csv(path2_eq, index=False)
    df_case3_eq.to_csv(path3_eq, index=False)

    return df_case2, df_case3


if __name__ == "__main__":
    main()
