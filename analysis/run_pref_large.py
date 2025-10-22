import os
import numpy as np
import pandas as pd
from model.equilibrium import pi_ge_eq, err_pe_large_tau10h, sweep_equilibrium_2country
from config import λ0_vals, bounds_lambda_8, pref_tax, G0, pref_tax, σ, τl

def main():
    """
    Perform the two-country sweep with large preferential τ10h in Country 0.
    """
    # 1) Baseline λ1 from symmetric one-country equilibrium
    λ1_star = pi_ge_eq[4]

    # 2) Build initial guess for two-country populations
    #    (concatenate one-country shares twice)
    x0_one = pi_ge_eq[:4]
    x0 = np.concatenate((x0_one, x0_one))

    # 3) Wrap err_pe_large_tau10h in the sweep
    # custom_τ10h comes from pref_tax in config
    df_pref_large = sweep_equilibrium_2country(
        err_func=err_pe_large_tau10h,
        λ0_vals=λ0_vals,
        λ1_fixed=λ1_star,
        bounds=bounds_lambda_8,
        x0=x0,
        custom_τ10h=pref_tax     # unilateral; do NOT pass tau_high_imm_override here
    )


    # 4) Save and display results
    os.makedirs("output/tables", exist_ok=True)
    #out_csv = f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}.csv"
    # out_tex = f"output/tables/two_country_pref_{int(pref_tax*100)}_sig_{int(σ)}_large.tex"

    #df_pref_large.to_csv(out_csv, index=False)

    
    # after you build df
    eq = df_pref_large[df_pref_large["TR0"] >= G0].head(1)
    eq.to_csv(f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv", index=False)

    pd.set_option('display.max_columns', None)
    print("Two-country preferential τ10h head:\n", df_pref_large.head())
    return df_pref_large

if __name__ == "__main__":
    main()
