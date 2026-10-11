# tests/test_diff_dest.py
# Gates for the differentiated-destinations small-country solver
# (docs/model_reformulation_code_plan.md, Step 2).
import os

import numpy as np
import pytest

from config import σ, θ, pref_tax, G0, SIG_TAG, DD_TAG
from model.equilibrium import pi_ge_eq, err_ge_eq, p_foreign, λ_star, τ1
from model.equations import Pi
from model import equilibrium_diff_dest as dd

TWO_COUNTRY_CSV = f"output/tables/csv/two_country_pref_13_{SIG_TAG}_taul_0.216_equilibrium.csv"


@pytest.fixture(scope="session")
def policy_slsqp():
    return dd.solve_diff_dest(pref_tax)


@pytest.fixture(scope="session")
def policy_damped():
    return dd.solve_diff_dest_damped(pref_tax)


# --- Gate 1: residual identities at the baseline ---
def test_residual_identity_at_baseline():
    res_new = dd.err_diff_dest(pi_ge_eq, out=0)
    res_old = err_ge_eq(pi_ge_eq, out=0)
    # budget clearing, native low, native high coincide elementwise
    for i in (0, 1, 2):
        assert abs(res_new[i] - res_old[i]) < 1e-14


def test_location_density_equals_move_probability_at_baseline():
    for a in ('l', 'h'):
        q = dd.location_density(a, p_foreign, τ1)
        P_stay = Pi(0, 0, a, p_foreign, p_foreign, τ1, τ1)
        assert abs(q - (1.0 - P_stay)) < 1e-15


# --- Gate 2: zero-policy solve reproduces the baseline ---
def test_zero_policy_solve_warm_start():
    x = dd.solve_diff_dest(None)
    assert np.max(np.abs(x - pi_ge_eq)) < 1e-8


def test_zero_policy_solve_perturbed_start():
    x0 = np.concatenate([p_foreign * 1.05, [λ_star * 0.98]])
    x = dd.solve_diff_dest(None, x0_init=x0)
    assert np.max(np.abs(x - pi_ge_eq)) < 1e-6


# --- Gate 3: inflow elasticity θ at fixed baseline wages (exact only at γ = 0) ---
@pytest.mark.skipif(__import__("config").γ != 0, reason="net-wage elasticity is θ only at γ = 0; see test_theta.py")
def test_inflow_elasticity_theta():
    tau_h = τ1[3]
    h = 1e-6

    def q_h(tau):
        t = τ1.copy()
        t[3] = tau
        return dd.location_density('h', p_foreign, t)

    el = (np.log(q_h(tau_h - h * (1 - tau_h)))
          - np.log(q_h(tau_h + h * (1 - tau_h)))) / (2 * h)
    assert abs(el - θ) < 1e-6


# --- Gate 5: budget clearing and solver agreement at τ_P = 0.13 ---
def test_budget_clearing(policy_damped):
    row = dd.equilibrium_row(policy_damped[:4], policy_damped[4], τ10h_fixed=pref_tax)
    assert abs(row["TR"] - G0) < 1e-8


def test_solver_agreement(policy_slsqp, policy_damped):
    assert np.max(np.abs(policy_slsqp - policy_damped)) < 1e-6


# --- Gate 6 (soft): cross-checks ---
@pytest.mark.skipif(θ != σ, reason="binary-denominator cross-check has elasticity σ(1−s); comparable only at θ = σ")
def test_cross_check_fixed_foreign_netw(policy_slsqp):
    x_f = dd.solve_fixed_foreign_netw(pref_tax)
    resp_dd = policy_slsqp[3] - pi_ge_eq[3]   # high-skill immigrant stock response
    resp_f = x_f[3] - pi_ge_eq[3]
    # σ(1−s) vs σ: expected gap ~s_h = 0.155, allowed ~25%
    assert abs(resp_f / resp_dd - 1.0) < 0.25


# --- Gate 7: progressivity adjustment clears the budget at fixed λ* ---
def test_taul_adjustment_budget_clearing():
    from model.equilibrium import λ_star
    p0, tau_l_star, τ0 = dd.solve_diff_dest_taul(custom_τ10h=pref_tax)
    row = dd.equilibrium_row(p0, λ_star, τ_override=τ0)
    assert abs(row["TR"] - G0) < 1e-8
    # revenue lost to the preferential rate is recovered by raising progressivity
    assert tau_l_star > 0.216


# --- Gate 4: foreign CEV analytic, mirroring guard ---
def test_needs_policy_mirroring_guard():
    from analysis.build_equilibrium_tables import _needs_policy_mirroring
    assert _needs_policy_mirroring(f"diff_dest_pref_tax_13_{DD_TAG}_taul_0.216") is False
    assert _needs_policy_mirroring(f"soe_pref_tax_13_{SIG_TAG}_taul_0.216") is True


def test_foreign_cev_analytic_zero(policy_slsqp, tmp_path, monkeypatch):
    import pandas as pd
    import analysis.build_equilibrium_tables as bt
    from config import pref_tax as τP

    # exact B and P rows written to a temporary directory
    row_b = dd.equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4])
    row_p = dd.equilibrium_row(policy_slsqp[:4], policy_slsqp[4], τ10h_fixed=τP)
    b_path = tmp_path / f"diff_dest_baseline_{DD_TAG}_taul_0.216_equilibrium.csv"
    p_path = tmp_path / f"diff_dest_pref_tax_13_{DD_TAG}_taul_0.216_equilibrium.csv"
    pd.DataFrame([row_b]).to_csv(b_path, index=False)
    pd.DataFrame([row_p]).to_csv(p_path, index=False)

    # spy: the binary CES numerator must never be called for foreign origins (i=1)
    calls = []
    orig = bt._ces_num_for_i

    def spy(inc_map, a, i):
        calls.append(i)
        return orig(inc_map, a, i)

    monkeypatch.setattr(bt, "_ces_num_for_i", spy)

    _, _, tex_path = bt.build_equilibrium_tables_diff_dest(
        baseline_eq_path=str(b_path), pref_eq_path=str(p_path),
        out_dir=str(tmp_path), tex_dir=str(tmp_path),
    )
    assert calls and all(i == 0 for i in calls)

    tex = open(tex_path, encoding="utf-8").read()
    cev_line = next(ln for ln in tex.splitlines() if ln.startswith("{CEV}"))
    cells = [c.strip() for c in cev_line.rstrip("\\").split("&")[5:]]
    assert cells[1] == r"+0.0\%" and cells[3].startswith(r"+0.0\%")


@pytest.mark.skipif(not os.path.exists(TWO_COUNTRY_CSV),
                    reason="two-country PE equilibrium CSV not present")
def test_signs_match_two_country_pe(policy_slsqp):
    import pandas as pd
    tc = pd.read_csv(TWO_COUNTRY_CSV).iloc[0]
    resp_dd = policy_slsqp[:4] - pi_ge_eq[:4]
    resp_tc = np.array([tc["p00l"], tc["p10l"], tc["p00h"], tc["p10h"]]) - pi_ge_eq[:4]
    # compare signs only where both responses are away from zero
    mask = (np.abs(resp_dd) > 1e-4) & (np.abs(resp_tc) > 1e-4)
    assert mask.any()
    assert np.all(np.sign(resp_dd[mask]) == np.sign(resp_tc[mask]))
