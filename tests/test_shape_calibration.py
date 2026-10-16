# tests/test_shape_calibration.py
# The calibrated σ, θ (model/shape_params.json) reproduce the two target estimates
# (app:shape-parameters) to the precision of their three-decimal rounding.
import os

import pytest

pytestmark = pytest.mark.skipif(
    any(os.environ.get(k) for k in ("MIG_SIGMA", "MIG_THETA")),
    reason="σ, θ overridden by the environment")


def test_theta_matches_take_up_corrected_estimate():
    import config
    from analysis.shape_moments import arc_eF, E_F_TARGET
    assert abs(arc_eF(config.θ) - E_F_TARGET) < 2e-3


def test_marginal_and_arc_elasticities_agree_at_baseline():
    # at the baseline the one percent and the policy-size elasticities agree to
    # three decimals (app:shape-parameters); at ε = 20 they differ by about 0.01
    import config
    from analysis.shape_moments import arc_eF, marginal_eF
    if config.SHAPE_KEY != "eps=7,gamma=0,rho=3,phi=0.72":
        pytest.skip("baseline specification only")
    assert abs(marginal_eF(config.θ) - arc_eF(config.θ)) < 2e-3


def test_sigma_matches_pooled_logit_coefficient():
    import config
    from analysis.shape_moments import pooled_alpha, ALPHA_TARGET
    assert abs(pooled_alpha(config.θ) - ALPHA_TARGET) < 3e-3
