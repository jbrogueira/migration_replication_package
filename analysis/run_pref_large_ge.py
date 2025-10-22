import os
import pandas as pd
from model.equilibrium import sweep_equilibrium_2country
from model.equilibrium import ge_pref_residual_global
from config import λ0_vals, bounds_lambda_8, pref_tax, G0, σ, τl
from model.parameters import κl, κh 

def main():
    """
    Run full GE large preferential τ10h sweep, allowing λ1 to float.
    """
    df_pref_large_ge = sweep_equilibrium_2country(
        err_func=ge_pref_residual_global,   # <-- use the GLOBAL residual
        λ0_vals=λ0_vals,
        λ1_fixed=None,
        bounds=bounds_lambda_8,
        x0=None,
        free_lambda1=True,
        tau_high_imm_override=None,         # leave None here
        custom_τ10h=pref_tax                # 0.13 applied to τ10h & τ01h
    )


    # Save and show
    os.makedirs("output/tables", exist_ok=True)
    #out_csv = f"output/tables/csv/global_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}.csv"
    # out_tex = f"output/tables/pref_{int(pref_tax*100)}_sig_{int(σ)}_large_gbc.tex"

    #df_pref_large_ge.to_csv(out_csv, index=False)

    pd.set_option('display.max_columns', None)
    print("Two-country GE preferential τ10h head:\n", df_pref_large_ge.head())
    
    # after you build df
    eq = df_pref_large_ge[df_pref_large_ge["TR0"] >= G0].head(1)
    eq.to_csv(f"output/tables/csv/global_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv", index=False)
    
    return df_pref_large_ge


if __name__ == "__main__":
    main()
