# analysis/run_taul_sweeps.py
import os
import numpy as np
import pandas as pd

from config import σ, τl as TAU_L_BASE, pref_tax, G0
from model.equilibrium import (
    pi_ge_eq,
    sweep_taul_large_case2,
    sweep_taul_large_case3,
)

OUTPUT = os.path.join(os.getcwd(), "output", "tables", "csv")
os.makedirs(OUTPUT, exist_ok=True)

def main(tau_grid=None):
    lam_star = float(pi_ge_eq[4])

    if tau_grid is None:
        tau_grid = np.linspace(0.15, 0.34, 50)

    # Case 2: Country 0 pref (τ10h), Country 1 none
    df_case2 = sweep_taul_large_case2(
        lambda_star=lam_star,
        tau_grid=tau_grid,
        τ10h_pref=pref_tax
    )

    # Case 3: Bilateral pref (τ10h = τ01h)
    df_case3 = sweep_taul_large_case3(
        lambda_star=lam_star,
        tau_grid=tau_grid,
        τ_pref=pref_tax
    )

    # Save full sweeps
    #path2 = os.path.join(
    #     OUTPUT, f"two_country_taul_case2_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{TAU_L_BASE}.csv"
    # )
    # path3 = os.path.join(
    #     OUTPUT, f"two_country_taul_case3_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{TAU_L_BASE}.csv"
    # )
    # df_case2.to_csv(path2, index=False)
    # df_case3.to_csv(path3, index=False)

    # === Select equilibrium rows: first with TR0 > G0 (error if none) ===
    if "TR0" not in df_case2.columns or "TR0" not in df_case3.columns:
        raise ValueError("Expected 'TR0' column in both sweep DataFrames.")

    mask2 = df_case2["TR0"] > G0
    if not mask2.any():
        raise ValueError(f"Case 2: no row with TR0 > {G0:.3f}.")
    idx2 = df_case2.index[mask2][0]
    df_case2_eq = df_case2.loc[[idx2]].copy()

    mask3 = df_case3["TR0"] > G0
    if not mask3.any():
        raise ValueError(f"Case 3: no row with TR0 > {G0:.3f}.")
    idx3 = df_case3.index[mask3][0]
    df_case3_eq = df_case3.loc[[idx3]].copy()

    # Save equilibria (1-row)
    path2_eq = os.path.join(
        OUTPUT, f"two_country_taul_case2_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{TAU_L_BASE}_equilibrium.csv"
    )
    path3_eq = os.path.join(
        OUTPUT, f"two_country_taul_case3_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{TAU_L_BASE}_equilibrium.csv"
    )
    df_case2_eq.to_csv(path2_eq, index=False)
    df_case3_eq.to_csv(path3_eq, index=False)

    return df_case2, df_case3

if __name__ == "__main__":
    main()
