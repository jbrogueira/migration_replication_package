# analysis/run_taul_soe.py
import os
import numpy as np
import pandas as pd

from config import τl as TAU_L_BASE, pref_tax, σ, G0
from model.equilibrium import (
    pi_ge_eq,
    #sweep_taul_multilateral,
    #sweep_taul_soe,
    sweep_taul_soe_pref,
)

OUTPUT = os.path.join(os.getcwd(), "output", "tables", "csv")
os.makedirs(OUTPUT, exist_ok=True)

def main(tau_grid=None):
    lam_star   = float(pi_ge_eq[4])
    tau_l_base = float(TAU_L_BASE)
    if tau_grid is None:
        tau_grid = np.linspace(0.05, 0.40, 40)

    #df_mlt  = sweep_taul_multilateral(lam_star, tau_grid)
    #df_soe  = sweep_taul_soe(lam_star, tau_grid)
    df_pref = sweep_taul_soe_pref(lam_star, tau_grid, τ10h_pref=pref_tax)

    # --- save full sweeps ---
    #df_mlt.to_csv(os.path.join(OUTPUT, "taul_mlt_sig.csv"), index=False)
    #df_soe.to_csv(os.path.join(OUTPUT, "taul_soe_sig.csv"), index=False)
    #df_pref.to_csv(os.path.join(OUTPUT, f"taul_case1_soe_pref_tax_{int(pref_tax*100)}.csv"), index=False)

 # --- 1-row policy snapshot for LaTeX builder (CASE 1; pick by TR threshold) ---

    if "TR" not in df_pref.columns:
        raise ValueError("df_pref must have a 'TR' column for selection by revenue.")
    
    mask = df_pref["TR"] > G0
    if not mask.any():
        raise ValueError(f"No row with TR > {G0:.3f} found in df_pref.")
    
    idx_eq = df_pref.index[mask][0]  # first row with TR > G0
    
    df_pref_eq = df_pref.loc[[idx_eq]].copy()
    df_pref_eq.to_csv(
        os.path.join(
            OUTPUT,
            f"case1_soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{TAU_L_BASE}_equilibrium.csv"
        ),
        index=False
    )


    #return df_mlt, df_soe, df_pref
    return df_pref

if __name__ == "__main__":
    main()
