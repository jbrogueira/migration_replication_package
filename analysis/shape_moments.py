# analysis/shape_moments.py
"""
Model counterparts of the two estimates that set the shape parameters (app:shape-parameters).

θ: the elasticity of the stock of high-skill immigrants in a small country with respect
to their net-of-tax rate when the preferential rate τ_P is introduced, as an arc
elasticity in logs, with wages adjusting and the tax level of residents (λ0) held at
its baseline value. Target: 2.392, Kleven et al. (2014), Table II, panel A4, column 3
(foreigners in their first three years, controlling for imperfect take-up).

σ: the logit coefficient implied by the response of natives at home (both skills
pooled) to a one percent rise in the net-of-tax rate of every resident of a small
country (λ0 × e^0.01), with wages adjusting: d ln(natives at home) / d ln λ0, divided
by one minus the pooled staying probability. Target: 1.323, Kleven, Landais and Saez
(2013), Table 2, panel A, column 1.

σ is read from config at import (MIG_SIGMA); θ is an argument, so one process solves
θ for a given σ. analysis/calibrate_shape.py iterates on σ across processes.

Usage, from the package root:  MIG_SIGMA=1.275 MIG_THETA=10 python -m analysis.shape_moments
"""
import json
import sys

import numpy as np
from scipy import optimize

from config import σ, pref_tax, origin0_low, origin0_high
from model.equations import wage, taul
from model.equilibrium import p_foreign, λ_star
from model.equilibrium_diff_dest import err_diff_dest

E_F_TARGET = 2.392       # Kleven et al. (2014), Table II, A4, col. 3
ALPHA_TARGET = 1.323     # Kleven, Landais and Saez (2013), Table 2, A, col. 1

_base = np.asarray(p_foreign, dtype=float)
_tau10h_base = taul(wage(1, 0, 'h', _base), λ_star)


def _solve(λ0, τ10h=None, θ=3.0):
    """Country-0 allocation at a fixed λ0 (no budget clearing), small country."""
    f = lambda v: err_diff_dest(v, out=0, custom_λ0=λ0, custom_τ10h=τ10h, θ=θ)
    sol = optimize.root(f, _base.copy(), method='hybr', tol=1e-12)
    if not sol.success or np.max(np.abs(sol.fun)) > 1e-10:
        raise RuntimeError(f"solve failed at θ = {θ}: {sol.message}")
    return sol.x


def arc_eF(θ: float) -> float:
    """Arc elasticity of high-skill immigrants over the policy cut, λ0 fixed."""
    q = _solve(λ_star, pref_tax, θ)
    return float(np.log(q[3] / _base[3]) / np.log((1.0 - pref_tax) / (1.0 - _tau10h_base)))


def marginal_eF(θ: float, h: float = 0.01) -> float:
    """Same elasticity for a one percent rise in the net-of-tax rate."""
    q = _solve(λ_star, 1.0 - (1.0 - _tau10h_base) * np.exp(h), θ)
    return float(np.log(q[3] / _base[3]) / h)


def pooled_alpha(θ: float, h: float = 0.01) -> float:
    """Logit coefficient implied by natives at home, wages adjusting."""
    p = _solve(λ_star * np.exp(h), None, θ)
    nat_b, nat_p = _base[0] + _base[2], p[0] + p[2]
    P_bar = nat_b / (origin0_low + origin0_high)
    return float((np.log(nat_p / nat_b) / h) / (1.0 - P_bar))


def theta_star(target: float = E_F_TARGET) -> float:
    """θ at which the arc elasticity equals the target (θ ≥ σ)."""
    return float(optimize.brentq(lambda t: arc_eF(t) - target, max(σ, 1.0) + 1e-6, 15.0,
                                 xtol=1e-6))


def main():
    θs = theta_star()
    out = {"sigma": σ, "theta_star": θs, "arc_eF": arc_eF(θs),
           "marginal_eF": marginal_eF(θs), "alpha": pooled_alpha(θs)}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
