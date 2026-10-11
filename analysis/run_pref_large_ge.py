import os
from model.equilibrium import sweep_equilibrium_2country, ge_pref_residual_global
from analysis.exact_refine import refine_crossing
from config import λ0_vals, bounds_lambda_8, pref_tax, G0, σ, τl, SIG_TAG

def main():
    # Bind global pref tax; λ1 floats inside params[-1]; swallow extras
    ge_pref = lambda params, out, custom_λ0, custom_λ1, *_, **__: ge_pref_residual_global(
        params, out, custom_λ0, None, pref_tax
    )

    df_pref_large_ge = sweep_equilibrium_2country(
        err_func=ge_pref,
        λ0_vals=λ0_vals,
        λ1_fixed=None,
        bounds=bounds_lambda_8,
        x0=None,
        free_lambda1=True,
        # reporting-only
        tau_high_imm_override=pref_tax
    )

    os.makedirs("output/tables", exist_ok=True)
    out_csv = f"output/tables/csv/global_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}.csv"
    df_pref_large_ge.to_csv(out_csv, index=False)

    # exact budget clearing in country 0 (country 1's λ1 clears its own budget)
    cols = ["p00l", "p10l", "p00h", "p10h", "p11l", "p01l", "p11h", "p01h", "λ1"]
    warm = df_pref_large_ge[df_pref_large_ge["TR0"] >= G0].head(1)[cols].to_numpy()[0]
    eq = refine_crossing(
        df_pref_large_ge, "λ0", "TR0", G0,
        lambda l0: sweep_equilibrium_2country(ge_pref, [l0], None, bounds_lambda_8, x0=warm,
                                              free_lambda1=True, tau_high_imm_override=pref_tax))
    eq.to_csv(
        f"output/tables/csv/global_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        index=False
    )
    return df_pref_large_ge

if __name__ == "__main__":
    main()
