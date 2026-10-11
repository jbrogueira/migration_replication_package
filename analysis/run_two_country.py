import os
from model.equilibrium import pi_ge_eq, err_2country, sweep_equilibrium_2country
from config import λ0_vals, bounds_lambda_8, σ, τl, SIG_TAG

def main():
    λ1_star = pi_ge_eq[4]

    # Positional interface; swallow extras
    pe_residual = lambda params, out, custom_λ0, custom_λ1, *_, **__: err_2country(
        params=params, out=out, λ0=custom_λ0, λ1=custom_λ1, mode='PE'
    )

    df_two_country_pe = sweep_equilibrium_2country(
        err_func=pe_residual,
        λ0_vals=λ0_vals,
        λ1_fixed=λ1_star,
        bounds=bounds_lambda_8
    )

    CSV_DIR = os.path.join("output", "tables", "csv")
    os.makedirs(CSV_DIR, exist_ok=True)
    # out_csv = os.path.join(CSV_DIR, f"two_country_pe_{SIG_TAG}_taul_{τl}.csv")
    # df_two_country_pe.to_csv(out_csv, index=False, encoding="utf-8")
    return df_two_country_pe

if __name__ == "__main__":
    main()
