import os
from model.equilibrium import sweep_equilibrium, err_ge_eq_laffer
from config import λ0_vals, bounds_lambda, G0, σ, τl, SIG_TAG


def main():
    """
    Perform the multilateral (symmetric) sweep and save results.
    """
    df = sweep_equilibrium(err_ge_eq_laffer, λ0_vals, bounds_lambda, symmetric=True)

    # Save and display results
    os.makedirs("output/tables", exist_ok=True)
    #out_csv = f"output/tables/csv/baseline_{SIG_TAG}_taul_{τl}.csv"
    #df.to_csv(out_csv, index=False)
    # exact calibrated baseline (the grid row over-collects by up to one λ0 step)
    import pandas as pd
    from model.equilibrium import pi_ge_eq
    from model.equilibrium_diff_dest import equilibrium_row
    eq = pd.DataFrame([equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4])])
    eq.to_csv(f"output/tables/csv/baseline_{SIG_TAG}_taul_{τl}_equilibrium.csv", index=False)
    return df

if __name__ == "__main__":
    main()
