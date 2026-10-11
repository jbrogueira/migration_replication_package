# analysis/run_theta_sweep.py
# Small country with the preferential tax for several values of θ, the parameter of
# the choice among foreign destinations (appendix app:theta; overview table row).
# θ = σ: independent shocks across destinations; θ = inf: perfect substitutes.
# Every equilibrium clears the budget exactly (λ0 adjusts); foreign CEV is 0 for
# every θ (country 0 has measure zero in foreign choices).
import os

import numpy as np
import pandas as pd
from scipy import optimize

from model.equilibrium import pi_ge_eq, λ_star, tau_star
from model.equilibrium_diff_dest import (
    equilibrium_row, err_diff_dest, solve_diff_dest, sweep_diff_dest,
)
from analysis.run_diff_dest import _native_cev
from viz.plotting import plot_theta_laffer
from config import σ, θ as θ_base, λ0_vals, pref_tax, G0, τl, SIG_TAG

CSV_DIR = "output/tables/csv"
TEX_DIR = "output/tables/tex"
FIG_DIR = "output/figures"

THETA_GRID = (1, 2, 2.5, 3, 5, 10, 20, np.inf)   # 2.5 = baseline θ
THETA_FIG = (2.5, 5, 20, np.inf)


def _theta_str(θ) -> str:
    return r"$\infty$" if np.isinf(θ) else f"{θ:g}"


def _row(x, θ, row_b) -> dict:
    row = equilibrium_row(x[:4], x[4], τ10h_fixed=pref_tax,
                          fun=err_diff_dest(x, out=1, custom_τ10h=pref_tax, θ=θ))
    cev = _native_cev(row, row_b)
    return {"theta": θ, **row,
            "cev00l": cev['l'], "cev10l": 0.0, "cev00h": cev['h'], "cev10h": 0.0}


def solve_theta_path(thetas, row_b):
    """
    Budget-clearing equilibria along θ, warm-started (continuation). A θ with no
    budget-clearing solve is recorded with solved = False (e.g. at high ϵ in the
    perfect-substitutes limit the immigrant stock diverges before revenue falls to G0).
    """
    rows, sols, x = [], {}, None
    for θ in thetas:
        try:
            x = solve_diff_dest(pref_tax, x0_init=x, θ=θ)
        except AssertionError as e:
            print(f"DIAGNOSTIC FLAG: no budget-clearing equilibrium found at θ = {θ}: {e}")
            rows.append({"theta": θ, "solved": False})
            x = None
            continue
        sols[θ] = x
        rows.append({**_row(x, θ, row_b), "solved": True})
    return pd.DataFrame(rows), sols


def self_financing_theta(sols):
    """
    θ* with λ0(θ*) = λ*: the preferential tax clears the budget at the baseline tax
    level. Brackets on the finite grid; None if λ0 − λ* has no sign change there.
    """
    finite = sorted(θ for θ in sols if np.isfinite(θ))
    g = {θ: sols[θ][4] - λ_star for θ in finite}
    for lo, hi in zip(finite, finite[1:]):
        if np.sign(g[lo]) != np.sign(g[hi]):
            warm = {"x": sols[lo]}

            def _g(θ):
                warm["x"] = solve_diff_dest(pref_tax, x0_init=warm["x"], θ=θ)
                return warm["x"][4] - λ_star

            θ_star = optimize.brentq(_g, lo, hi, xtol=1e-8)
            return θ_star, solve_diff_dest(pref_tax, x0_init=warm["x"], θ=θ_star)
    return None, None


def cev_thresholds(df: pd.DataFrame, sols, row_b):
    """
    θ at which the native CEV of each skill turns from negative to positive (brentq
    between the bracketing grid points, warm-started). None if no sign change.
    """
    out = {}
    d = df[df["solved"] & np.isfinite(df["theta"])].sort_values("theta")
    for a in ("l", "h"):
        col = f"cev00{a}"
        out[a] = None
        for (_, lo), (_, hi) in zip(d.iterrows(), d.iloc[1:].iterrows()):
            if np.sign(lo[col]) != np.sign(hi[col]):
                warm = {"x": sols[lo["theta"]]}

                def _f(θ):
                    warm["x"] = solve_diff_dest(pref_tax, x0_init=warm["x"], θ=θ)
                    return _row(warm["x"], θ, row_b)[col]

                out[a] = optimize.brentq(_f, lo["theta"], hi["theta"], xtol=1e-8)
                break
    return out


def _tex_table(df: pd.DataFrame, θ_star, path: str):
    from analysis.build_equilibrium_tables import _pct_str as pct
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\begin{tabular}{ c | c c | c c c c }",
        r"\toprule",
        r" & & & \multicolumn{4}{c}{CEV} \\",
        r"$\theta$ & $\pi_{10}^h$ & $1-\lambda_0$ & $c_{00}^l$ & $c_{10}^l$ & $c_{00}^h$ & $c_{10}^h$ \\",
        r"\midrule",
        fr"B & {pi_ge_eq[3]:.3f} & {1 - λ_star:.3f} & & & & \\",
    ]
    for _, r in df.iterrows():
        if not r["solved"]:
            lines.append(fr"{_theta_str(r['theta'])} & \multicolumn{{6}}{{c}}{{no budget-clearing equilibrium}} \\")
            continue
        lines.append(fr"{_theta_str(r['theta'])} & {r['p10h']:.3f} & {1 - r['λ0']:.3f} & "
                     f"{pct(r['cev00l'])} & {pct(r['cev10l'])} & {pct(r['cev00h'])} & {pct(r['cev10h'])} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    note = (f"The preferential tax clears the budget at the baseline tax level at $\\theta = {θ_star:.2f}$."
            if θ_star is not None else
            "The preferential tax does not clear the budget at the baseline tax level for any finite $\\theta$ on the grid.")
    lines += [
        fr"\caption{{Unilateral preferential tax ($\tau^P = {pref_tax}$) in the small country, "
        fr"for values of $\theta$ ($\sigma = {σ:g}$). B: baseline. {note}}}",
        fr"\label{{tab:theta_sweep_{SIG_TAG}}}",
        r"\end{table}",
    ]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    for d in (CSV_DIR, TEX_DIR, FIG_DIR):
        os.makedirs(d, exist_ok=True)
    tag = f"pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}"

    row_b = equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4])
    thetas = [θ for θ in THETA_GRID if θ >= σ]
    df, sols = solve_theta_path(thetas, row_b)

    θ_star, x_star = self_financing_theta(sols)
    if θ_star is not None:
        df_star = pd.DataFrame([_row(x_star, θ_star, row_b)])
        df_star.to_csv(f"{CSV_DIR}/theta_star_{tag}.csv", index=False)
        print(f"self-financing θ* = {θ_star:.4f} (λ0 = λ* = {λ_star:.4f}); "
              f"native CEV l {df_star['cev00l'].iloc[0]:+.2%}, h {df_star['cev00h'].iloc[0]:+.2%}")
    else:
        print("no self-financing θ on the finite grid")

    thr = cev_thresholds(df, sols, row_b)
    pd.DataFrame([{"theta_star": θ_star, "theta_cev0_l": thr["l"], "theta_cev0_h": thr["h"]}]).to_csv(
        f"{CSV_DIR}/theta_thresholds_{tag}.csv", index=False)
    print("native CEV turns positive at θ =",
          {a: (None if v is None else round(v, 3)) for a, v in thr.items()})

    df.to_csv(f"{CSV_DIR}/theta_sweep_{tag}.csv", index=False)
    _tex_table(df, θ_star, f"{TEX_DIR}/theta_sweep_{tag}.tex")

    # Laffer curves with the preferential tax for a subset of θ
    df_base = sweep_diff_dest(λ0_vals, θ=θ_base if θ_base >= σ else σ)
    dfs = {}
    for θ in THETA_FIG:
        if θ < σ or θ not in sols:
            continue
        lab = (r"Preferential tax, $\theta = \infty$ (perfect substitutes)" if np.isinf(θ)
               else fr"Preferential tax, $\theta = {θ:g}$")
        dfs[lab] = sweep_diff_dest(λ0_vals, custom_τ10h=pref_tax, x0=sols[θ][:4], θ=θ)
    plot_theta_laffer(dfs, df_base, tau_star, G0,
                      filename=f"{FIG_DIR}/theta_laffer_{SIG_TAG}_taul_{τl}.png")

    return df, θ_star


if __name__ == "__main__":
    main()
