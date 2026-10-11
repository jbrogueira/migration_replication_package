# tests/test_theta.py
# Gates for θ in the small-country solver (docs/nested_theta_code_plan.md §3, §6).
import numpy as np
import pytest

from config import σ, pref_tax
from model.equilibrium import pi_ge_eq, p_foreign, τ1
from model.equations import Pi
from model import equilibrium_diff_dest as dd

THETAS = sorted({σ, 2.0, 5.0, 20.0} - {x for x in (2.0, 5.0, 20.0) if x < σ})


def _old_density(a, p0, τ0):
    # pre-θ expression, V_move^σ / (V_stay^σ + D)
    return dd.move_value(a, p0, τ0) ** σ / (dd.dd_baseline.V_stay[a] ** σ + dd.dd_baseline.D[a])


def test_theta_equal_sigma_identity():
    x = dd.solve_diff_dest(pref_tax, θ=σ)
    τ0 = np.array([dd.equilibrium_row(x[:4], x[4], τ10h_fixed=pref_tax)[k]
                   for k in ("tau00l", "tau10l", "tau00h", "tau10h")])
    for p, τ in ((p_foreign, τ1), (x[:4], τ0)):
        for a in ('l', 'h'):
            q_new = dd.location_density(a, p, τ, θ=σ)
            assert abs(q_new - _old_density(a, p, τ)) < 1e-14


@pytest.mark.parametrize("θ", THETAS + [np.inf])
def test_zero_policy_invariant_in_theta(θ):
    x = dd.solve_diff_dest(None, θ=θ)
    assert np.max(np.abs(x - pi_ge_eq)) < 1e-8


@pytest.mark.parametrize("θ", THETAS)
def test_inflow_elasticity_equals_theta(θ):
    # elasticity of the inflow with respect to the value of moving to country 0, at
    # fixed wages; with respect to the net wage it is θ times the labor-income share
    # of that value (1 at γ = 0)
    h = 1e-6
    tau_h = τ1[3]

    def q_and_v(tau):
        t = τ1.copy()
        t[3] = tau
        return (dd.location_density('h', p_foreign, t, θ=θ),
                dd.move_value('h', p_foreign, t))

    (q_lo, v_lo), (q_hi, v_hi) = q_and_v(tau_h + h), q_and_v(tau_h - h)
    el = (np.log(q_hi) - np.log(q_lo)) / (np.log(v_hi) - np.log(v_lo))
    assert abs(el - θ) < 1e-5 * max(1.0, θ)


def test_native_staying_probability_independent_of_theta():
    # the native residuals do not involve θ; the foreign staying probability is θ-free
    for a in ('l', 'h'):
        P = Pi(0, 0, a, p_foreign, p_foreign, τ1, τ1)
        assert abs(dd.dd_baseline.P_stay[a] - P) < 1e-14


def test_limit_theta_to_infinity():
    # continuation in θ: the policy equilibrium approaches the θ = inf (indifference) solve
    try:
        x_inf = dd.solve_diff_dest(pref_tax, θ=np.inf)
    except AssertionError:
        pytest.skip("no budget-clearing equilibrium in the perfect-substitutes limit "
                    "for this specification")
    x = None
    gaps = []
    for θ in (σ, 5.0, 20.0, 50.0, 100.0, 200.0):
        if θ < σ:
            continue
        x = dd.solve_diff_dest(pref_tax, x0_init=x, θ=θ)
        gaps.append(np.max(np.abs(x[:4] - x_inf[:4])))
    assert all(g2 < g1 for g1, g2 in zip(gaps, gaps[1:]))
    # convergence is slow (at θ = 200 the gap is 0.24 of its θ = σ value at σ = 1,
    # 0.31 at σ = 2): the immigrant stock ratio is (V_move/v̄)^θ, so the value gap
    # closes only like 1/θ
    assert gaps[-1] < 0.5 * gaps[0]


def test_theta_inf_matches_soe_solver():
    # the θ = inf solve at the SOE grid equilibrium's λ0 reproduces its allocation
    from model.equilibrium import err_ge_pref_tax
    from scipy import optimize
    from config import bounds_lambda
    λ0 = 0.9
    sol = optimize.minimize(lambda v: err_ge_pref_tax(v, 1, λ0, pref_tax), p_foreign,
                            method='SLSQP', bounds=bounds_lambda,
                            tol=1e-14, options={"ftol": 1e-14, "maxiter": 1000})
    sw = dd.sweep_diff_dest([λ0], custom_τ10h=pref_tax, θ=np.inf)
    p_dd = sw[["p00l", "p10l", "p00h", "p10h"]].iloc[0].to_numpy()
    # relative: at high ϵ the stocks are large and weakly pinned by the wage condition
    assert np.max(np.abs(p_dd / sol.x - 1.0)) < 1e-4
