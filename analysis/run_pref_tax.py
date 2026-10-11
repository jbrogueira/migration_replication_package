import os
import numpy as np
import pandas as pd
from scipy import optimize
from model.equilibrium import sweep_equilibrium, err_ge_pref_tax, pi_ge_eq
from config import λ0_vals, pref_tax, G0, σ, bounds_lambda, τl, tau10h_vals, SIG_TAG
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
    csv_share = f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}.csv"
    df_share.to_csv(csv_share, index=False)
    

    # This part of the script just makes sure we look for the correct
    # equilibrium in cases where total revenues are not monotonic
    if df_share["TR"].is_monotonic_increasing:
        eq = df_share[df_share["TR"] >= G0].head(1)
    else:
        peak_idx = df_share["TR"].idxmax()
        tail = df_share.loc[peak_idx:]  
        eq = tail.loc[(tail["TR"] >= G0) & (tail["TR"].shift(-1) < G0)].head(1) 
        if eq.empty:
            ge = tail[tail["TR"] >= G0].tail(1)
            eq = ge if not ge.empty else df_share.iloc[0:0]


    # diagnostic flag: the selected row is a budget-clearing equilibrium only if the
    # revenue curve crosses G0 between it and its grid neighbour; otherwise the row
    # is the edge of the solvable λ0 range (e.g. the immigrant stock diverges before
    # revenue falls to G0) and is not an equilibrium
    crossing = False
    if not eq.empty:
        i = df_share.index.get_loc(eq.index[0])
        step = abs(λ0_vals[1] - λ0_vals[0])
        for j in (i - 1, i + 1):
            if 0 <= j < len(df_share):
                nb = df_share.iloc[j]
                adjacent = abs(nb["λ0"] - df_share.iloc[i]["λ0"]) < 1.5 * step
                if adjacent and np.sign(nb["TR"] - G0) != np.sign(df_share.iloc[i]["TR"] - G0):
                    crossing = True
    if crossing:
        # exact budget clearing (θ = inf: immigrant indifference), warm-started at the grid row
        from model.equilibrium_diff_dest import solve_diff_dest, equilibrium_row
        r0 = eq.iloc[0]
        x_init = np.array([r0["p00l"], r0["p10l"], r0["p00h"], r0["p10h"], r0["λ0"]])
        x = solve_diff_dest(pref_tax, x0_init=x_init, θ=np.inf)
        eq = pd.DataFrame([equilibrium_row(x[:4], x[4], τ10h_fixed=pref_tax)])
    eq = eq.assign(budget_crossing=crossing)
    if not crossing:
        print(f"DIAGNOSTIC FLAG: no budget-clearing crossing of TR = G0 on the solvable "
              f"λ0 range for the perfect-substitutes preferential tax (τ_P = {pref_tax}); "
              f"the equilibrium CSV row is NOT an equilibrium")

    eq.to_csv(f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv", index=False)

    # 2) Sweep over τ10h for Laffer curve
    λ_star = pi_ge_eq[4]
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

    #csv_tau = "output/tables/csv/soe_tau10h_sweep_{SIG_TAG}_taul_{τl}.csv"
    #df_tau.to_csv(csv_tau, index=False)
    return df_share, df_tau


if __name__ == "__main__":
    main()
