"""
Small country with differentiated destinations (main text ssec:sc-diff, eq:q-nested,
eq:q-small-nested; appendix app:nested).

Country 0 is a measure-zero point on a measure-one continuum of symmetric foreign
countries. A worker draws a Frechet(σ) taste shock for the home country and, for the
foreign countries, a common taste for living abroad times country-specific
Poisson-Frechet(θ) shocks (θ ≥ σ). Foreign wages, taxes, the foreign stay value, and
the inclusive value D are unaffected by country-0 policy and are frozen at the
symmetric baseline (pi_ge_eq, solved in model.equilibrium).

Inflow from the (symmetric) rest of the world into country 0:
    q^a(0) = (1 − P^a_stay) (V_move / v̄^a)^θ,
with P^a_stay the foreign natives' staying probability (independent of θ and of
country-0 policy) and v̄^a the move value to any other foreign country. At θ = σ this
is V_move^σ / (V_stay^σ + D). θ = inf is the perfect-substitutes limit: the
immigrant-stock condition becomes indifference, V_move = v̄^a.

Unknowns: country-0 composition [p00l, p10l, p00h, p10h] and, when the budget
clears, λ0. Residual ordering in err_diff_dest: [budget clearing (only when λ0 is
free), native stay l, native stay h, immigrant stock l, immigrant stock h].

err_fixed_foreign_netw / solve_fixed_foreign_netw are cross-check-only (binary
denominators with foreign net wages frozen at baseline, inflow elasticity σ(1−s);
ported from code/model_update-large-economy_fixed wages.ipynb, cells 17-18).
"""
from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import pandas as pd
from scipy import optimize

from .equations import wage, taul, Pi, ri, land_income_origin, cap_tax_revenue
from .parameters import κl, κh
from .equilibrium import pi_ge_eq, p_foreign, λ_star, foreign_w, τ1
from config import σ, θ as θ_default, G0, tau_k, ell_l, ell_h, origin0_low, origin0_high, bnds, bounds_lambda


# --- Foreign-side objects frozen at the symmetric baseline ---
@dataclass(frozen=True)
class DiffDestBaseline:
    C: Dict[str, float]       # origin land income of a foreign-born worker (zero at γ = 0)
    V_stay: Dict[str, float]  # foreign native stay value
    v_bar: Dict[str, float]   # move value to any foreign destination (symmetric)
    D: Dict[str, float]       # inclusive value of the foreign continuum, v_bar^σ for every θ
    P_stay: Dict[str, float]  # foreign natives' staying probability


def _compute_dd_baseline() -> DiffDestBaseline:
    C, V_stay, v_bar, D, P_stay = {}, {}, {}, {}, {}
    for a, i_nat, i_imm, κ_a in (('l', 0, 1, κl), ('h', 2, 3, κh)):
        C[a] = land_income_origin(1, a, p_foreign)
        V_stay[a] = foreign_w[(0, a)] * (1.0 - τ1[i_nat]) + C[a]
        # measure-one continuum of identical destinations: the integral collapses
        v_bar[a] = (foreign_w[(1, a)] * (1.0 - τ1[i_imm]) + C[a]) / κ_a
        D[a] = v_bar[a] ** σ
        P_stay[a] = V_stay[a] ** σ / (V_stay[a] ** σ + D[a])
    return DiffDestBaseline(C=C, V_stay=V_stay, v_bar=v_bar, D=D, P_stay=P_stay)


dd_baseline = _compute_dd_baseline()


# --- Location density q^a(0) ---
def move_value(a: str, p0, τ0) -> float:
    """Value of moving to country 0 for a foreign-born worker: consumption (net wage plus
    origin land income) divided by κ, as in Pi."""
    κ_a = κl if a == 'l' else κh
    tau_imm = τ0[1] if a == 'l' else τ0[3]
    w_move_net = wage(1, 0, a, p0) * (1.0 - tau_imm)
    return max((w_move_net + dd_baseline.C[a]) / κ_a, 1e-12)


def location_density(a: str, p0, τ0, θ: float = θ_default) -> float:
    """
    q^a(0) = (1 − P^a_stay) (V_move / v̄^a)^θ, eq:q-small-nested (finite θ).
    P^a_stay and v̄^a are frozen at baseline: country 0 has measure zero abroad.
    """
    if np.isinf(θ):
        raise ValueError("θ = inf has no density; use the indifference condition (err_diff_dest)")
    ratio = move_value(a, p0, τ0) / dd_baseline.v_bar[a]
    return (1.0 - dd_baseline.P_stay[a]) * ratio ** θ


def _immigrant_residuals(p0, τ0, θ):
    """Immigrant-stock conditions [low, high]; indifference in the limit θ = inf."""
    if np.isinf(θ):
        return [move_value('l', p0, τ0) - dd_baseline.v_bar['l'],
                move_value('h', p0, τ0) - dd_baseline.v_bar['h']]
    return [p0[1] - location_density('l', p0, τ0, θ) * origin0_low,
            p0[3] - location_density('h', p0, τ0, θ) * origin0_high]


# --- Residual ---
def required_spending(p0, g_per_capita: float = False) -> float:
    """Government spending G0·[(1 − s) + s·N], N = resident population, s = g_per_capita
    the share of spending proportional to the population (False/0: constant G0;
    True/1: G0·N). Equals G0 at the baseline, where the resident population is one."""
    s = float(g_per_capita)
    return G0 * ((1.0 - s) + s * float(np.sum(p0[:4])))


def err_diff_dest(params, out: int = 1, custom_λ0: Optional[float] = None,
                  custom_τ10h: Optional[float] = None, θ: float = θ_default,
                  g_per_capita: float = False):
    p0 = np.asarray(params, dtype=float)[:4]
    λ0 = params[4] if custom_λ0 is None else custom_λ0

    w00l = wage(0, 0, 'l', p0); w10l = wage(1, 0, 'l', p0)
    w00h = wage(0, 0, 'h', p0); w10h = wage(1, 0, 'h', p0)
    τ00l = taul(w00l, λ0); τ10l = taul(w10l, λ0)
    τ00h = taul(w00h, λ0)
    τ10h = float(custom_τ10h) if custom_τ10h is not None else taul(w10h, λ0)
    τ0 = np.array([τ00l, τ10l, τ00h, τ10h])

    res = [
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * origin0_low,   # natives low
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * origin0_high,  # natives high
        *_immigrant_residuals(p0, τ0, θ),                             # immigrants l, h
    ]

    if custom_λ0 is None:
        rev_labor0 = (p0[0]*τ00l*w00l + p0[1]*τ10l*w10l +
                      p0[2]*τ00h*w00h + p0[3]*τ10h*w10h)
        rev_cap0 = cap_tax_revenue(p0, ri(0, p0), ri(1, p_foreign))
        res.insert(0, rev_labor0 + rev_cap0 - required_spending(p0, g_per_capita))  # budget clearing

    return np.sum(np.square(res)) if out == 1 else res


# --- Shared CSV-row helper (same column set as sweep_equilibrium) ---
def equilibrium_row(p0, λ0: float, τ10h_fixed: Optional[float] = None,
                    fun: float = np.nan, τ_override=None) -> dict:
    """τ_override = [τ00l, τ10l, τ00h, τ10h] replaces the λ0-schedule rates
    (used by the progressivity-adjustment exercise, where τl clears the budget)."""
    w00l = wage(0, 0, 'l', p0); w10l = wage(1, 0, 'l', p0)
    w00h = wage(0, 0, 'h', p0); w10h = wage(1, 0, 'h', p0)

    if τ_override is not None:
        tau00l, tau10l, tau00h, tau10h = (float(t) for t in τ_override)
    else:
        tau00l = taul(w00l, λ0); tau10l = taul(w10l, λ0)
        tau00h = taul(w00h, λ0)
        tau10h = float(τ10h_fixed) if τ10h_fixed is not None else taul(w10h, λ0)

    # origin-based land returns; the foreign return is frozen at baseline
    r0 = ri(0, p0)
    r1 = ri(1, p_foreign)

    land00l = r0 * ell_l * (1.0 - tau_k); land10l = r1 * ell_l * (1.0 - tau_k)
    land00h = r0 * ell_h * (1.0 - tau_k); land10h = r1 * ell_h * (1.0 - tau_k)

    inc00l = w00l*(1 - tau00l) + land00l
    inc10l = w10l*(1 - tau10l) + land10l
    inc00h = w00h*(1 - tau00h) + land00h
    inc10h = w10h*(1 - tau10h) + land10h

    rev_labor0 = (p0[0]*tau00l*w00l + p0[1]*tau10l*w10l +
                  p0[2]*tau00h*w00h + p0[3]*tau10h*w10h)
    rev_cap0 = cap_tax_revenue(p0, r0, r1)
    TR = rev_labor0 + rev_cap0

    return {
        "λ0": λ0,
        "p00l": p0[0], "p10l": p0[1], "p00h": p0[2], "p10h": p0[3],
        "w00l": w00l, "w10l": w10l, "w00h": w00h, "w10h": w10h,
        "tau00l": tau00l, "tau10l": tau10l, "tau00h": tau00h, "tau10h": tau10h,
        "inc00l": inc00l, "inc10l": inc10l, "inc00h": inc00h, "inc10h": inc10h,
        "TR": TR, "fun": fun,
    }


# --- Solvers ---
def solve_diff_dest(custom_τ10h: Optional[float] = None, x0_init=None,
                    θ: float = θ_default, g_per_capita: float = False) -> np.ndarray:
    """Joint budget-clearing solve; returns [p00l, p10l, p00h, p10h, λ0]."""
    x0i = (np.asarray(x0_init, dtype=float) if x0_init is not None
           else np.concatenate([p_foreign, [λ_star]]))
    sol = optimize.minimize(
        lambda v: err_diff_dest(v, out=1, custom_τ10h=custom_τ10h, θ=θ,
                                g_per_capita=g_per_capita),
        x0i, method='SLSQP', bounds=bnds, tol=1e-12,
        options={"ftol": 1e-12, "maxiter": 1000},
    )
    assert sol.success, f"diff-dest solve failed (θ={θ}): {sol.message}"
    assert sol.fun < 1e-8, f"diff-dest residual too large (θ={θ}): {sol.fun:.2e}"
    # polish on the residual vector (square system)
    pol = optimize.root(lambda v: err_diff_dest(v, out=0, custom_τ10h=custom_τ10h, θ=θ,
                                              g_per_capita=g_per_capita),
                        sol.x, method='hybr', tol=1e-14)
    res = lambda v: np.max(np.abs(err_diff_dest(v, out=0, custom_τ10h=custom_τ10h, θ=θ,
                                                g_per_capita=g_per_capita)))
    return pol.x if res(pol.x) < res(sol.x) else sol.x


def solve_diff_dest_damped(custom_τ10h: Optional[float] = None, damping: float = 0.5,
                           tol: float = 1e-12, maxit: int = 2000,
                           θ: float = θ_default, g_per_capita: float = False) -> np.ndarray:
    """
    Cross-validation solve: inner damped iteration on the stock map, outer brentq
    on budget clearing g(λ0) = TR − G, with G = G0 or G0 times the resident population.
    Returns [p00l, p10l, p00h, p10h, λ0].
    """
    if np.isinf(θ):
        raise ValueError("damped stock iteration needs a finite θ")
    warm = {"p": np.asarray(p_foreign, dtype=float).copy()}

    def _stocks_at(λ0: float) -> np.ndarray:
        p = warm["p"].copy()
        for _ in range(maxit):
            τ0 = np.array([
                taul(wage(0, 0, 'l', p), λ0),
                taul(wage(1, 0, 'l', p), λ0),
                taul(wage(0, 0, 'h', p), λ0),
                (float(custom_τ10h) if custom_τ10h is not None
                 else taul(wage(1, 0, 'h', p), λ0)),
            ])
            target = np.array([
                Pi(0, 0, 'l', p, p_foreign, τ0, τ1) * origin0_low,
                location_density('l', p, τ0, θ) * origin0_low,
                Pi(0, 0, 'h', p, p_foreign, τ0, τ1) * origin0_high,
                location_density('h', p, τ0, θ) * origin0_high,
            ])
            p_new = p + damping * (target - p)
            if np.max(np.abs(p_new - p)) < tol:
                warm["p"] = p_new.copy()
                return p_new
            p = p_new
        raise RuntimeError(f"damped stock iteration did not converge at λ0={λ0}")

    def _g(λ0: float) -> float:
        p = _stocks_at(λ0)
        return (equilibrium_row(p, λ0, τ10h_fixed=custom_τ10h)["TR"]
                - required_spending(p, g_per_capita))

    # the revenue curve can cross G twice at high σ: scan λ0 descending from 1
    # and bracket the first crossing — the sweep-selection convention (first
    # grid point from the top with TR ≥ G)
    λ_grid = np.linspace(1.0, 0.75, 26)
    g_prev = _g(λ_grid[0])
    for λ_lo, g_here in ((λ, _g(λ)) for λ in λ_grid[1:]):
        if np.sign(g_here) != np.sign(g_prev):
            λ0 = optimize.brentq(_g, λ_lo, λ_lo + (λ_grid[0] - λ_grid[1]),
                                 xtol=1e-12)
            return np.concatenate([_stocks_at(λ0), [λ0]])
        g_prev = g_here
    raise RuntimeError("no budget-clearing crossing of TR = G in [0.75, 1.0]")


def sweep_diff_dest(λ_vals, custom_τ10h: Optional[float] = None, x0=None,
                    θ: float = θ_default) -> pd.DataFrame:
    """Laffer sweep over λ0 (no budget clearing), warm-started along the sweep."""
    curr = (np.asarray(x0, dtype=float)[:4] if x0 is not None
            else np.asarray(p_foreign, dtype=float).copy())
    rows = []
    for λ in λ_vals:
        sol = optimize.minimize(
            lambda v: err_diff_dest(v, out=1, custom_λ0=λ, custom_τ10h=custom_τ10h, θ=θ),
            curr, method='SLSQP', bounds=bounds_lambda, tol=1e-12,
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        if not sol.success or sol.fun > 1e-8:
            continue
        p0 = sol.x[:4]
        rows.append(equilibrium_row(p0, λ, τ10h_fixed=custom_τ10h, fun=sol.fun))
        curr = p0.copy()
    return pd.DataFrame(rows)


# --- Progressivity adjustment: λ0 fixed at λ*, τl clears the budget ---
def err_diff_dest_taul(params, out: int = 1, tau_l_param: Optional[float] = None,
                       custom_τ10h: Optional[float] = None, θ: float = θ_default):
    """Stock residuals under the progressivity schedule 1 − λ*·y^(−tau_l_param)."""
    from .equilibrium import _taul_param
    p0 = np.asarray(params, dtype=float)[:4]

    w00l = wage(0, 0, 'l', p0); w10l = wage(1, 0, 'l', p0)
    w00h = wage(0, 0, 'h', p0); w10h = wage(1, 0, 'h', p0)
    τ0 = np.array([
        _taul_param(w00l, λ_star, tau_l_param),
        _taul_param(w10l, λ_star, tau_l_param),
        _taul_param(w00h, λ_star, tau_l_param),
        (float(custom_τ10h) if custom_τ10h is not None
         else _taul_param(w10h, λ_star, tau_l_param)),
    ])

    res = [
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * origin0_low,
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * origin0_high,
        *_immigrant_residuals(p0, τ0, θ),
    ]
    return np.sum(np.square(res)) if out == 1 else res


def solve_diff_dest_taul(custom_τ10h: Optional[float] = None,
                         tau_bracket=(0.05, 0.40), θ: float = θ_default):
    """
    Budget clearing through the progressivity parameter (λ0 = λ*).
    Returns (p0, tau_l_star, τ0) — stocks, the clearing τl, and the tax rates.
    """
    from .equilibrium import _taul_param
    warm = {"x": np.asarray(p_foreign, dtype=float).copy()}

    def _solve_at(tau_l: float) -> np.ndarray:
        sol = optimize.minimize(
            lambda v: err_diff_dest_taul(v, out=1, tau_l_param=tau_l,
                                         custom_τ10h=custom_τ10h, θ=θ),
            warm["x"], method='SLSQP', bounds=bounds_lambda, tol=1e-12,
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        assert sol.success and sol.fun < 1e-8, \
            f"taul solve failed at tau_l={tau_l}: fun={sol.fun:.2e}"
        warm["x"] = sol.x.copy()
        return sol.x[:4]

    def _taxes(p0: np.ndarray, tau_l: float) -> np.ndarray:
        w = [wage(0, 0, 'l', p0), wage(1, 0, 'l', p0),
             wage(0, 0, 'h', p0), wage(1, 0, 'h', p0)]
        τ0 = np.array([_taul_param(wi, λ_star, tau_l) for wi in w])
        if custom_τ10h is not None:
            τ0[3] = float(custom_τ10h)
        return τ0

    def _g(tau_l: float) -> float:
        p0 = _solve_at(tau_l)
        row = equilibrium_row(p0, λ_star, fun=np.nan, τ_override=_taxes(p0, tau_l))
        return row["TR"] - G0

    tau_l_star = optimize.brentq(_g, *tau_bracket, xtol=1e-12)
    p0 = _solve_at(tau_l_star)
    return p0, tau_l_star, _taxes(p0, tau_l_star)


# --- Cross-check only: binary denominators, foreign net wages frozen at baseline ---
# No land income terms (the notebook computation predates land; valid at γ = 0).
_netw1 = {
    ('nat', 'l'): foreign_w[(0, 'l')] * (1.0 - τ1[0]),
    ('imm', 'l'): foreign_w[(1, 'l')] * (1.0 - τ1[1]),
    ('nat', 'h'): foreign_w[(0, 'h')] * (1.0 - τ1[2]),
    ('imm', 'h'): foreign_w[(1, 'h')] * (1.0 - τ1[3]),
}


def err_fixed_foreign_netw(params, out: int = 1, custom_λ0: Optional[float] = None,
                           custom_τ10h: Optional[float] = None):
    p0 = np.asarray(params, dtype=float)[:4]
    λ0 = params[4] if custom_λ0 is None else custom_λ0

    res = []
    for a, mass, κ_a, idx in (('l', origin0_low, κl, (0, 1)),
                              ('h', origin0_high, κh, (2, 3))):
        w00 = wage(0, 0, a, p0)
        u_stay0 = w00 * (1.0 - taul(w00, λ0))
        u_abroad = _netw1[('imm', a)]
        P00 = u_stay0**σ / (u_stay0**σ + (u_abroad / κ_a)**σ)

        w10 = wage(1, 0, a, p0)
        τ10 = (float(custom_τ10h) if (a == 'h' and custom_τ10h is not None)
               else taul(w10, λ0))
        u_move0 = w10 * (1.0 - τ10)
        u_stay1 = _netw1[('nat', a)]
        P10 = (u_move0 / κ_a)**σ / ((u_move0 / κ_a)**σ + u_stay1**σ)

        res += [p0[idx[0]] - P00 * mass, p0[idx[1]] - P10 * mass]

    if custom_λ0 is None:
        TR = equilibrium_row(p0, λ0, τ10h_fixed=custom_τ10h)["TR"]
        res.insert(0, TR - G0)

    return np.sum(np.square(res)) if out == 1 else res


def solve_fixed_foreign_netw(custom_τ10h: Optional[float] = None, x0_init=None) -> np.ndarray:
    """Cross-check only; not part of the reported equilibrium."""
    x0i = (np.asarray(x0_init, dtype=float) if x0_init is not None
           else np.concatenate([p_foreign, [λ_star]]))
    sol = optimize.minimize(
        lambda v: err_fixed_foreign_netw(v, out=1, custom_τ10h=custom_τ10h),
        x0i, method='SLSQP', bounds=bnds, tol=1e-12,
        options={"ftol": 1e-12, "maxiter": 1000},
    )
    assert sol.success, f"fixed-foreign-netw solve failed: {sol.message}"
    assert sol.fun < 1e-8, f"fixed-foreign-netw residual too large: {sol.fun:.2e}"
    return sol.x
