import os
from model.equilibrium import sweep_equilibrium, err_ge_soe_laffer, pi_ge_eq
from config import λ0_vals, bounds_lambda, σ, τl, SIG_TAG


def main():
    """
    Perform the unilateral (SOE) sweep and save results, starting
    from the symmetric-equilibrium shares to force divergence.
    """
    # initialize at symmetric GE equilibrium shares (first four elements)
    x0_soe = pi_ge_eq[:4]

    # pass the custom initial guess to the sweep
    df = sweep_equilibrium(
        err_ge_soe_laffer,
        λ0_vals,
        bounds_lambda,
        x0=x0_soe
    )

    os.makedirs("output/tables", exist_ok=True)
    csv_path = f"output/tables/csv/soe_{SIG_TAG}_taul_{τl}.csv"
    df.to_csv(csv_path, index=False)

    return df

if __name__ == "__main__":
    main()
