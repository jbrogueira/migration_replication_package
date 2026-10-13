"""Gates for government spending proportional to the resident population (TR = G0 · N)."""
import numpy as np

from config import G0, G1, pref_tax, λ0_vals, bounds_lambda_8
from model.equilibrium import pi_ge_eq
from model.equilibrium_diff_dest import (
    equilibrium_row, err_diff_dest, solve_diff_dest, solve_diff_dest_damped,
)

TAU_P = pref_tax


def test_zero_policy_reproduces_baseline():
    x = solve_diff_dest(None, g_per_capita=True)
    assert np.max(np.abs(x - pi_ge_eq)) < 1e-6


def test_default_rule_unchanged():
    x = np.concatenate([pi_ge_eq[:4], [0.9]])
    assert err_diff_dest(x, out=0, custom_τ10h=TAU_P) == \
        err_diff_dest(x, out=0, custom_τ10h=TAU_P, g_per_capita=False)


def test_budget_clears_at_g0_times_population():
    x = solve_diff_dest(TAU_P, g_per_capita=True)
    TR = equilibrium_row(x[:4], x[4], τ10h_fixed=TAU_P)["TR"]
    assert abs(TR - G0 * np.sum(x[:4])) < 1e-9


def test_slsqp_and_damped_agree():
    x_p = solve_diff_dest(TAU_P, g_per_capita=True)
    x_d = solve_diff_dest_damped(TAU_P, g_per_capita=True)
    assert np.max(np.abs(x_p - x_d)) < 1e-6


def test_two_country_zero_policy_baseline():
    from model.equilibrium import err_2country
    x = np.concatenate([pi_ge_eq[:4], pi_ge_eq[:4], [pi_ge_eq[4]]])
    res = err_2country(x, out=0, λ0=pi_ge_eq[4], g_per_capita=True)
    assert np.max(np.abs(res)) < 1e-6


def test_perfect_substitutes_budget_clears_at_g0_times_population():
    from analysis.run_g_per_capita import _solve_perfect_substitutes
    from model.equilibrium_diff_dest import sweep_diff_dest
    df = sweep_diff_dest(λ0_vals, custom_τ10h=TAU_P, θ=np.inf)
    x = _solve_perfect_substitutes(df, TAU_P, g_per_capita=True)
    TR = equilibrium_row(x[:4], x[4], τ10h_fixed=TAU_P)["TR"]
    assert np.sum(x[:4]) > 1.01
    assert abs(TR - G0 * np.sum(x[:4])) < 1e-9


def test_perfect_substitutes_no_crossing_is_flagged():
    import pandas as pd
    from analysis.run_g_per_capita import _solve_perfect_substitutes
    df = pd.DataFrame({"λ0": λ0_vals[:3], "TR": [1.0, 1.1, 1.2],
                       "p00l": 0.5, "p10l": 0.1, "p00h": 0.3, "p10h": 0.1})
    assert _solve_perfect_substitutes(df, TAU_P, g_per_capita=False) is None


def test_two_country_country1_spending_scales_with_population():
    """Unilateral adoption at the baseline λ0: country 1 loses residents (N1 < 1) and
    clears its budget at G1·N1."""
    from model.equilibrium import sweep_equilibrium_2country, ge_pref_residual_unilateral
    err = lambda p, out, l0, l1, *_, **__: ge_pref_residual_unilateral(
        p, out, l0, None, TAU_P, g_per_capita=True)
    x0 = np.concatenate([pi_ge_eq[:4], pi_ge_eq[:4], [pi_ge_eq[4]]])
    r = sweep_equilibrium_2country(err, [pi_ge_eq[4]], None, bounds_lambda_8, x0=x0,
                                   free_lambda1=True, tau10h_override=TAU_P).iloc[0]
    N1 = r[["p11l", "p01l", "p11h", "p01h"]].sum()
    assert N1 < 0.999
    assert abs(r["TR1"] - G1 * N1) < 1e-9


# --- Share s of spending proportional to the population: TR = G0·[(1 − s) + s·N] ---
S_RIVAL = 0.62


def test_share_one_and_zero_match_the_two_rules():
    x = np.concatenate([pi_ge_eq[:4] * 1.03, [0.9]])
    assert err_diff_dest(x, out=0, custom_τ10h=TAU_P, g_per_capita=1.0) == \
        err_diff_dest(x, out=0, custom_τ10h=TAU_P, g_per_capita=True)
    assert err_diff_dest(x, out=0, custom_τ10h=TAU_P, g_per_capita=0.0) == \
        err_diff_dest(x, out=0, custom_τ10h=TAU_P, g_per_capita=False)


def test_share_zero_policy_reproduces_baseline():
    x = solve_diff_dest(None, g_per_capita=S_RIVAL)
    assert np.max(np.abs(x - pi_ge_eq)) < 1e-6


def test_share_budget_clears_and_lies_between_the_two_rules():
    x = solve_diff_dest(TAU_P, g_per_capita=S_RIVAL)
    N = np.sum(x[:4])
    TR = equilibrium_row(x[:4], x[4], τ10h_fixed=TAU_P)["TR"]
    assert abs(TR - G0 * ((1 - S_RIVAL) + S_RIVAL * N)) < 1e-9
    lam_fixed = solve_diff_dest(TAU_P, g_per_capita=False)[4]
    lam_pc = solve_diff_dest(TAU_P, g_per_capita=True)[4]
    assert min(lam_fixed, lam_pc) - 1e-9 <= x[4] <= max(lam_fixed, lam_pc) + 1e-9


def test_share_slsqp_and_damped_agree():
    x_p = solve_diff_dest(TAU_P, g_per_capita=S_RIVAL)
    x_d = solve_diff_dest_damped(TAU_P, g_per_capita=S_RIVAL)
    assert np.max(np.abs(x_p - x_d)) < 1e-6
