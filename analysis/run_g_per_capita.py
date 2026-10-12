# analysis/run_g_per_capita.py
# Preferential tax with government spending proportional to the resident population:
# each budget clears at TR_i = G_i · N_i, N_i = resident population of country i (N_i = 1
# at the baseline). The tax levels λ adjust; everything else as in the main experiments.
#   (1) small country, differentiated destinations (θ finite)
#   (2) small country, perfect-substitutes limit (θ = inf)
#   (3) large country, unilateral adoption (country 1's spending also scales with N_1)
#   (4) multilateral adoption (symmetric: N_0 = N_1 = 1, the equilibrium is unchanged)
import argparse
import os
from typing import Optional

import numpy as np
import pandas as pd

from model.equilibrium import (
    pi_ge_eq, tau_star, sweep_equilibrium_2country,
    ge_pref_residual_unilateral, ge_pref_residual_global,
)
from model.equilibrium_diff_dest import (
    err_diff_dest, equilibrium_row, required_spending, solve_diff_dest,
    solve_diff_dest_damped, sweep_diff_dest,
)
from analysis.exact_refine import refine_crossing
from analysis.build_equilibrium_tables import (
    build_equilibrium_tables, build_equilibrium_tables1, build_equilibrium_tables_diff_dest,
)
from viz.plotting import plot_g_per_capita_laffer
from config import (pref_tax as pref_tax_default, λ0_vals, G0, G1, bounds_lambda_8,
                    DD_TAG, SIG_TAG, τl)

CSV_DIR = "output/tables/csv"
FIG_DIR = "output/figures"
STOCKS0 = ['p00l', 'p10l', 'p00h', 'p10h']
STOCKS1 = ['p11l', 'p01l', 'p11h', 'p01h']
TABLE_TAIL = (r"baseline~(B) vs.\ policy~(P). Migration shares, wages, tax rates, "
              "consumption, and CEV welfare effects.}")


# --- (1)-(2) Small country -------------------------------------------------
def _solve_exact(pref_tax: float, g_per_capita: bool) -> np.ndarray:
    """SLSQP + polish, cross-validated against the damped/brentq solve."""
    x_p = solve_diff_dest(pref_tax, g_per_capita=g_per_capita)
    x_d = None
    for damping in (0.5, 0.15, 0.05):
        try:
            x_d = solve_diff_dest_damped(pref_tax, damping=damping, g_per_capita=g_per_capita)
            break
        except RuntimeError:
            continue
    assert x_d is not None, "damped solve failed at all damping values"
    assert np.max(np.abs(x_p - x_d)) < 1e-6, "SLSQP and damped solves disagree"
    return x_p


def _solve_perfect_substitutes(df_pref: pd.DataFrame, pref_tax: float,
                               g_per_capita: bool) -> Optional[np.ndarray]:
    """
    θ = inf: exact budget clearing from the unique crossing of the λ0 sweep.
    Returns None (with a diagnostic flag) when the sweep does not cross the budget
    exactly once between adjacent grid points, e.g. at ε = 20 with constant spending,
    where no budget-clearing equilibrium exists.
    """
    N = df_pref[STOCKS0].sum(axis=1)
    above = (df_pref['TR'] - (G0 * N if g_per_capita else G0) >= 0).to_numpy()
    λ = df_pref['λ0'].to_numpy()
    step = abs(λ0_vals[1] - λ0_vals[0])
    idx = np.where(above[1:] != above[:-1])[0]
    adjacent = [i for i in idx if abs(λ[i + 1] - λ[i]) < 1.5 * step]
    rule = "G0·N" if g_per_capita else "G0"
    if len(idx) != 1 or len(adjacent) != 1:
        print(f"DIAGNOSTIC FLAG: perfect substitutes, G = {rule}: {len(idx)} sign changes of "
              f"TR − G on the solvable λ0 grid ({len(adjacent)} between adjacent points); "
              f"no equilibrium reported")
        return None
    i = adjacent[0]
    r = df_pref.iloc[i]
    x = solve_diff_dest(pref_tax, x0_init=np.r_[r[STOCKS0].to_numpy(float), r['λ0']],
                        θ=np.inf, g_per_capita=g_per_capita)
    res = err_diff_dest(x, out=0, custom_τ10h=pref_tax, θ=np.inf, g_per_capita=g_per_capita)
    assert np.max(np.abs(res)) < 1e-10, "perfect-substitutes solve did not converge"
    assert min(λ[i], λ[i + 1]) <= x[4] <= max(λ[i], λ[i + 1]), \
        "perfect-substitutes solve left the bracketing grid interval"
    return x


def _small_row(x, pref_tax, g_per_capita, θ=None):
    kw = {} if θ is None else {"θ": θ}
    row = equilibrium_row(x[:4], x[4], τ10h_fixed=pref_tax,
                          fun=err_diff_dest(x, out=1, custom_τ10h=pref_tax,
                                            g_per_capita=g_per_capita, **kw))
    row["N"] = float(np.sum(x[:4]))
    row["G"] = required_spending(x, g_per_capita)
    assert abs(row["TR"] - row["G"]) < 1e-8, "budget does not clear"
    return row


def small_country(pref_tax: float, tag: str) -> dict:
    # zero-policy gate: at N = 1 the per-resident rule reproduces the baseline
    x_zero = solve_diff_dest(None, g_per_capita=True)
    assert np.max(np.abs(x_zero - pi_ge_eq)) < 1e-6, "zero-policy gate failed"
    row_b = equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4], fun=err_diff_dest(pi_ge_eq, out=1))
    base_path = f"{CSV_DIR}/diff_dest_baseline_{DD_TAG}_taul_{τl}_equilibrium.csv"
    pd.DataFrame([row_b]).to_csv(base_path, index=False)

    rows = {}
    # (1) differentiated destinations
    for name, gpc in (("fixedG", False), ("gpc", True)):
        rows[f"dd_{name}"] = _small_row(_solve_exact(pref_tax, gpc), pref_tax, gpc)
        pd.DataFrame([rows[f"dd_{name}"]]).to_csv(
            f"{CSV_DIR}/diff_dest_{name}_{tag}_equilibrium.csv", index=False)
    build_equilibrium_tables_diff_dest(
        baseline_eq_path=base_path,
        pref_eq_path=f"{CSV_DIR}/diff_dest_gpc_{tag}_equilibrium.csv",
        caption=("{Unilateral preferential tax in the small country with differentiated "
                 "destinations, government spending proportional to the resident population: "
                 + TABLE_TAIL),
        label="tab:g-per-capita",
    )
    df_base = sweep_diff_dest(λ0_vals)
    df_pref = sweep_diff_dest(λ0_vals, custom_τ10h=pref_tax)
    plot_g_per_capita_laffer(
        df_base, df_pref, tau_star, G0,
        eq_fixed=(rows["dd_fixedG"]["λ0"], rows["dd_fixedG"]["N"]),
        eq_gpc=(rows["dd_gpc"]["λ0"], rows["dd_gpc"]["N"]),
        filename=os.path.join(FIG_DIR, f"g_per_capita_laffer_{tag}.png"),
    )

    # (2) perfect-substitutes limit
    df_base_ps = sweep_diff_dest(λ0_vals, θ=np.inf)
    df_pref_ps = sweep_diff_dest(λ0_vals, custom_τ10h=pref_tax, θ=np.inf)
    for name, gpc in (("fixedG", False), ("gpc", True)):
        x = _solve_perfect_substitutes(df_pref_ps, pref_tax, gpc)
        rows[f"ps_{name}"] = None if x is None else _small_row(x, pref_tax, gpc, θ=np.inf)
    if rows["ps_gpc"] is not None:
        # file name keeps the SOE mirroring rule of build_equilibrium_tables (foreign side frozen)
        ps_path = (f"{CSV_DIR}/gpc_soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}"
                   f"_equilibrium.csv")
        pd.DataFrame([rows["ps_gpc"]]).to_csv(ps_path, index=False)
        build_equilibrium_tables(
            pref_eq_path=ps_path,
            caption=("{Unilateral preferential tax in the small country in the perfect-substitutes "
                     "limit, government spending proportional to the resident population: "
                     + TABLE_TAIL),
            label="tab:g-per-capita-ps",
        )
        eq_of = lambda r: None if r is None else (r["λ0"], r["N"])
        plot_g_per_capita_laffer(
            df_base_ps, df_pref_ps, tau_star, G0,
            eq_fixed=eq_of(rows["ps_fixedG"]), eq_gpc=eq_of(rows["ps_gpc"]),
            filename=os.path.join(
                FIG_DIR, f"g_per_capita_laffer_ps_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}.png"),
        )
    return rows


# --- (3)-(4) Two large countries --------------------------------------------
def _two_country(residual, pref_tax: float, g_per_capita: bool, tau_kw: dict, tau_pref):
    """λ0 sweep with λ1 clearing country 1's budget, then the exact crossing of
    country 0's budget (revenue per resident = G0 under the per-resident rule)."""
    err = lambda p, out, l0, l1, *_, **__: residual(p, out, l0, None, tau_pref,
                                                    g_per_capita=g_per_capita)
    x0 = np.concatenate([pi_ge_eq[:4], pi_ge_eq[:4], [pi_ge_eq[4]]])
    add_pc = lambda d: d.assign(TR0pc=d['TR0'] / d[STOCKS0].sum(axis=1))
    df = add_pc(sweep_equilibrium_2country(err, λ0_vals, None, bounds_lambda_8, x0=x0,
                                           free_lambda1=True, **tau_kw))
    trcol = 'TR0pc' if g_per_capita else 'TR0'
    warm = df[df[trcol] >= G0].head(1)[STOCKS0 + STOCKS1 + ['λ1']].to_numpy()[0]
    eq = refine_crossing(
        df, 'λ0', trcol, G0,
        lambda l0: add_pc(sweep_equilibrium_2country(err, [l0], None, bounds_lambda_8, x0=warm,
                                                     free_lambda1=True, **tau_kw)))
    r = eq.iloc[0]
    N0, N1 = r[STOCKS0].sum(), r[STOCKS1].sum()
    assert abs(r['TR0'] - G0 * (N0 if g_per_capita else 1.0)) < 1e-7, "country 0 budget"
    assert abs(r['TR1'] - G1 * (N1 if g_per_capita else 1.0)) < 1e-7, "country 1 budget"
    return df, eq.assign(N0=N0, N1=N1)


def large_countries(pref_tax: float) -> dict:
    out = {}

    # (3) unilateral adoption; the no-preferential-tax sweep is drawn for reference
    uni = dict(tau10h_override=pref_tax)
    _, out["lc_fixedG"] = _two_country(ge_pref_residual_unilateral, pref_tax, False, uni, pref_tax)
    df_pref, out["lc_gpc"] = _two_country(ge_pref_residual_unilateral, pref_tax, True, uni, pref_tax)
    df_none, _ = _two_country(ge_pref_residual_unilateral, pref_tax, True, {}, None)
    lc_path = f"{CSV_DIR}/gpc_two_country_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv"
    out["lc_gpc"].to_csv(lc_path, index=False)
    build_equilibrium_tables(
        pref_eq_path=lc_path,
        caption=("{Unilateral preferential tax in a large country, government spending "
                 "proportional to the resident population, country~0: " + TABLE_TAIL),
        label="tab:g-per-capita-lc",
    )
    build_equilibrium_tables1(
        pref_eq_path=lc_path,
        caption=("{Unilateral preferential tax in a large country, government spending "
                 "proportional to the resident population, country~1: " + TABLE_TAIL),
        label="tab:g-per-capita-lc-cty1",
    )
    to_small = lambda d: d.rename(columns={'TR0': 'TR'})
    r_f, r_g = out["lc_fixedG"].iloc[0], out["lc_gpc"].iloc[0]
    plot_g_per_capita_laffer(
        to_small(df_none), to_small(df_pref), tau_star, G0,
        eq_fixed=(r_f['λ0'], r_f['N0']), eq_gpc=(r_g['λ0'], r_g['N0']),
        filename=os.path.join(FIG_DIR, f"g_per_capita_laffer_lc_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}.png"),
    )

    # (4) multilateral adoption: symmetric, so N_0 = N_1 = 1 and the equilibrium is unchanged
    glob = dict(tau_high_imm_override=pref_tax)
    _, out["mlt_fixedG"] = _two_country(ge_pref_residual_global, pref_tax, False, glob, pref_tax)
    _, out["mlt_gpc"] = _two_country(ge_pref_residual_global, pref_tax, True, glob, pref_tax)
    cols = ['λ0', 'λ1'] + STOCKS0 + STOCKS1
    gap = np.max(np.abs(out["mlt_gpc"][cols].to_numpy() - out["mlt_fixedG"][cols].to_numpy()))
    assert gap < 1e-6, f"multilateral: per-resident rule moves the equilibrium ({gap:.1e})"
    out["mlt_gpc"].to_csv(
        f"{CSV_DIR}/gpc_global_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        index=False)
    return out


def main(pref_tax: float = pref_tax_default):
    os.makedirs(CSV_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    tag = f"pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}"
    # baseline row read by build_equilibrium_tables (written by run_baseline in run_all)
    if not os.path.exists(f"{CSV_DIR}/baseline_{SIG_TAG}_taul_{τl}_equilibrium.csv"):
        from analysis.run_baseline import main as run_baseline
        run_baseline()
    rows = small_country(pref_tax, tag)
    rows.update(large_countries(pref_tax))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pref-tax", type=float, default=pref_tax_default)
    args = ap.parse_args()
    rows = main(args.pref_tax)
    for k in ("dd_fixedG", "dd_gpc", "ps_fixedG", "ps_gpc"):
        r = rows[k]
        if r is None:
            print(f"{k:10s} no budget-clearing equilibrium")
            continue
        print(f"{k:10s} 1-λ0 = {1 - r['λ0']:.4f}  N = {r['N']:.4f}  G = {r['G']:.4f}")
    for k in ("lc_fixedG", "lc_gpc", "mlt_fixedG", "mlt_gpc"):
        r = rows[k].iloc[0]
        print(f"{k:10s} 1-λ0 = {1 - r['λ0']:.4f}  1-λ1 = {1 - r['λ1']:.4f}  "
              f"N0 = {r['N0']:.4f}  N1 = {r['N1']:.4f}")
