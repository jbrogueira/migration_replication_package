import os
import numpy as np
import pandas as pd
from scipy import optimize
from model.equilibrium import sweep_equilibrium, err_ge_pref_tax, pi_ge_eq
from config import λ0_vals, pref_tax, G0, σ, bounds_lambda, τl
from model.equations import wage, taul


def main(pref_tax=pref_tax, x0=pi_ge_eq[:4], bounds=bounds_lambda):
    """
    1) Sweep equilibrium fixing high-skill tax at pref_tax.
    2) Sweep over τ10h and compute total tax revenue.
    """
    # 1) Share sweep
    wrapped = lambda p, out, custom_λ0: err_ge_pref_tax(
        params=p,
        out=out,
        custom_λ0=custom_λ0,
        custom_τ10h=pref_tax
    )
    wrapped.tau10h_fixed = pref_tax  
    df_share = sweep_equilibrium(
        err_func=wrapped,
        λ_vals=λ0_vals,
        bounds=bounds,
        x0=x0
    )

    os.makedirs("output/tables", exist_ok=True)
    #csv_share = f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}.csv"
    # tex_share = f"output/tables/tex/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}.tex"

    #df_share.to_csv(csv_share, index=False)
    
    # Equilibrium values with preferential tax
    eq = df_share[df_share["TR"] >= G0].head(1)
    print(eq.shape)   
    eq.to_csv(f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv", index=False)

    # 2) Sweep over τ10h for Laffer curve
    λ_star = pi_ge_eq[4]
    tau10h_vals = np.linspace(0, 0.5, 25)
    results = []
    for τ10h in tau10h_vals:
        sol = optimize.minimize(
            lambda x: err_ge_pref_tax(
                params=x,
                out=1,
                custom_λ0=λ_star,
                custom_τ10h=τ10h
            ),
            x0, method='SLSQP', bounds=bounds, tol=1e-10
        )
        if not sol.success:
            continue

        p0 = sol.x
        w = np.array([
            wage(0, 0, 'l', p0),
            wage(1, 0, 'l', p0),
            wage(0, 0, 'h', p0),
            wage(1, 0, 'h', p0),
        ])
        τ_arr = np.array([taul(w[i], λ_star) for i in range(3)] + [τ10h])
        TR = np.sum(p0 * τ_arr * w)
        results.append({'tau10h': τ10h, 'TR': TR})

    df_tau = pd.DataFrame(results)

    #csv_tau = "output/tables/csv/soe_tau10h_sweep_sig_{int(σ)}_taul_{τl}.csv"
    #tex_tau = "output/tables/tau10h_sweep.tex"

    #df_tau.to_csv(csv_tau, index=False)
    #df_tau.round(2).to_latex(
    #    tex_tau, index=False, longtable=True, float_format="%.2f", encoding="utf-8"
    #)
    #print("Pref-tax tau10h head:\n", df_tau.head())

    return df_share, df_tau


if __name__ == "__main__":
    main()
