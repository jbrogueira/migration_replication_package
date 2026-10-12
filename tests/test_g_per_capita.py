"""Gates for government spending proportional to the resident population (TR = G0 · N)."""
import numpy as np

from config import G0
from model.equilibrium import pi_ge_eq
from model.equilibrium_diff_dest import (
    equilibrium_row, err_diff_dest, solve_diff_dest, solve_diff_dest_damped,
)

TAU_P = 0.13


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
