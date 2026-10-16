# tests/test_calibration.py
# Gates for the calibration (docs/nested_theta_code_plan.md §2).
import numpy as np
import pytest

import config
from config import TARGETS, CALIBRATION, σ, ρ, ϵ, φ, τl, S, G0, γ
from model import calibration as cal
from model.equations import wage, Nia, Fi, Pi, taul
from model.equilibrium import pi_ge_eq, err_ge_eq

pytestmark = pytest.mark.skipif(config.LEGACY_CALIB, reason="legacy calibration active")


def _alloc():
    return np.array(CALIBRATION.alloc)


def test_wage_formulas_match_equations():
    p = _alloc()
    w_cal, F_cal, r_cal = cal.symmetric_wages(p, S, ρ, ϵ, φ, γ)
    w_eq = [wage(0, 0, 'l', p), wage(1, 0, 'l', p), wage(0, 0, 'h', p), wage(1, 0, 'h', p)]
    assert np.max(np.abs(w_cal - np.array(w_eq))) < 1e-14
    assert abs(F_cal - Fi(0, p)) < 1e-14
    from model.equations import ri
    assert abs(r_cal - ri(0, p)) < 1e-14


def test_calibrated_point_is_an_equilibrium():
    x = np.concatenate([_alloc(), [CALIBRATION.lam]])
    assert np.max(np.abs(err_ge_eq(x, out=0))) < 1e-13


def test_baseline_solve_returns_targets():
    x = np.concatenate([_alloc(), [CALIBRATION.lam]])
    assert np.max(np.abs(pi_ge_eq - x)) < 1e-10


def test_targets_hit():
    p = pi_ge_eq[:4]
    f = TARGETS.foreign_share
    assert abs(p[1] + p[3] - f) < 1e-10                           # immigrant share
    assert abs(p[3] - f * TARGETS.tertiary_foreign) < 1e-10       # high-skill immigrants
    Nl, Nh = Nia(p[0], p[1]), Nia(p[2], p[3])
    bill_l = p[0] * wage(0, 0, 'l', p) + p[1] * wage(1, 0, 'l', p)
    bill_h = p[2] * wage(0, 0, 'h', p) + p[3] * wage(1, 0, 'h', p)
    assert abs((bill_h / Nh) / (bill_l / Nl) - TARGETS.skill_premium) < 1e-10
    assert abs(G0 / Fi(0, p) - TARGETS.g_share) < 1e-10


def test_kappa_inverts_staying_probability():
    p = pi_ge_eq[:4]
    τ = np.array([taul(wage(j, 0, a, p), pi_ge_eq[4])
                  for j, a in ((0, 'l'), (1, 'l'), (0, 'h'), (1, 'h'))])
    assert abs(Pi(0, 0, 'l', p, p, τ, τ) - p[0] / config.origin0_low) < 1e-10
    assert abs(Pi(0, 0, 'h', p, p, τ, τ) - p[2] / config.origin0_high) < 1e-10


def test_theta_below_sigma_rejected():
    import subprocess, sys, os
    env = dict(os.environ, MIG_SIGMA="2", MIG_THETA="1")
    r = subprocess.run([sys.executable, "-c", "import config"], env=env,
                       capture_output=True, text=True)
    assert r.returncode != 0 and "θ ≥ σ" in r.stderr


def test_land_owned_by_natives_equals_land_stock():
    # land market clears: native endowments add up to the land stock (app:land)
    assert abs(config.ell_l * config.origin0_low + config.ell_h * config.origin0_high
               - config.L0) < 1e-12
    assert abs(config.ell_h / config.ell_l - TARGETS.skill_premium) < 1e-12


def test_phi_recalibrated_to_baseline_wage_ratio_when_eps_differs():
    # φ keeps the mean immigrant/native wage ratio of the baseline (φ = 0.72, ε = 7)
    import config
    from model.calibration import immigrant_wage_ratio, phi_for_wage_ratio
    p = config.TARGETS.allocation()[0]
    target = immigrant_wage_ratio(p, 7, 0.72)
    assert abs(target - 0.90) < 0.005
    phi20 = phi_for_wage_ratio(p, 20, target)
    assert abs(immigrant_wage_ratio(p, 20, phi20) - target) < 1e-12
    assert abs(round(phi20, 3) - 0.833) < 1e-12
