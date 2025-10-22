# analysis/run_two_country.py
import os
import pandas as pd
from model.equilibrium import pi_ge_eq, err_2country, sweep_equilibrium_2country
from config import λ0_vals, bounds_lambda_8, σ, τl

def main():
    """
    Perform the partial-equilibrium two-country sweep and save results.
    """
    # 1) Baseline λ1 from symmetric one-country equilibrium
    λ1_star = pi_ge_eq[4]

    # 2) Wrap err_2country in PE mode
    pe_residual = lambda params, out, custom_λ0, custom_λ1: err_2country(
        params=params,
        out=out,
        λ0=custom_λ0,
        λ1=custom_λ1,
        tau10h_override=None,
        mode='PE'
    )

    # 3) Run the sweep over λ0 holding λ1 fixed
    df_two_country_pe = sweep_equilibrium_2country(
        err_func=pe_residual,
        λ0_vals=λ0_vals,
        λ1_fixed=λ1_star,
        bounds=bounds_lambda_8
    )

    # 4) Save and display results
    CSV_DIR = os.path.join("output", "tables", "csv")
    os.makedirs(CSV_DIR, exist_ok=True)
    #out_csv = os.path.join(CSV_DIR, f"two_country_pe_sig_{int(σ)}_taul_{τl}.csv")

    #df_two_country_pe.to_csv(out_csv, index=False, encoding="utf-8")

    pd.set_option('display.max_columns', None)
    #print("Saved:", out_csv)
    print("Two-country PE head:\n", df_two_country_pe.head())
    return df_two_country_pe

if __name__ == "__main__":
    main()
