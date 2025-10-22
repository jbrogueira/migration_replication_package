import os
import pandas as pd
from model.equilibrium import sweep_equilibrium, err_ge_eq_laffer
from config import λ0_vals, bounds_lambda, G0, σ, τl 


def main():
    """
    Perform the multilateral (symmetric) sweep and save results.
    """
    df = sweep_equilibrium(err_ge_eq_laffer, λ0_vals, bounds_lambda)

    os.makedirs("output/tables", exist_ok=True)
    #out_csv = f"output/tables/csv/baseline_sig_{int(σ)}_taul_{τl}.csv"
    #out_tex = f"output/tables/tex/baseline.tex"
    #out_eq_csv = "output/tables/baseline_equilibrium.csv"

    #df.to_csv(out_csv, index=False)
    # df.round(2).to_latex(
    #     out_tex, index=False, longtable=True, float_format="%.2f", encoding="utf-8"
    # )

    # after you build df
    eq = df[df["TR"] >= G0].head(1)
    eq.to_csv(f"output/tables/csv/baseline_sig_{int(σ)}_taul_{τl}_equilibrium.csv", index=False)
    return df


if __name__ == "__main__":
    main()
