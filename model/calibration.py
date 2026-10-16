# model/calibration.py
"""
Calibration of the symmetric baseline (sec:calib, app:calibration).

The targets fix the baseline allocation directly, and every parameter then follows
in closed form, in sequence:
  1. allocation: immigrant stocks π_{10}^a and native stocks π_{00}^a = Π^a − π_{10}^a;
  2. S from the college wage premium, SP = S^{(ρ−1)/ρ} (N^l/N^h)^{1/ρ};
  3. G = G/Y target × output; λ from the government budget constraint,
     λ = (Σ π w + R_k − G) / Σ π w^{1−τ}, with R_k = τ_k r Σ_a ℓ^a Π^a the
     residence-based tax on land income (at the symmetric baseline every resident's
     origin rent is r; zero at γ = 0);
  4. κ^a from the staying probability. At a symmetric allocation all foreign
     payoffs are equal, so D_θ = V_move^σ for every θ and
     P^a = V_stay^σ / (V_stay^σ + V_move^σ), V_stay = c_{jj} + C^a,
     V_move = (c_{j,-j} + C^a)/κ^a (C^a = r ℓ^a (1 − τ_k), origin land income; κ
     divides consumption, eq:ut), which inverts to
     κ^a = (c_{j,-j} + C^a) / (V_stay [(1 − P^a)/P^a]^{1/σ}); at γ = 0,
     κ^a = (c_{j,-j}/c_{jj}) [P^a/(1 − P^a)]^{1/σ}.
Only step 4 involves σ; θ enters nowhere. Every specification (σ, ε, γ) is calibrated
to the same data targets.

Self-contained (no config import): config.py calls it at import. The wage formulas
replicate model/equations.py with A = L = 1 (tests/test_calibration.py checks that
the two agree).
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Targets:
    foreign_share: float      # foreign-born share of the 25–64 population
    tertiary_native: float    # tertiary share among natives
    tertiary_foreign: float   # tertiary share among the foreign-born
    skill_premium: float      # college / non-college average earnings
    g_share: float            # G / Y

    def allocation(self):
        """Symmetric allocation [π00l, π10l, π00h, π10h] and origin masses (Π^l, Π^h)."""
        f = self.foreign_share
        p10h = f * self.tertiary_foreign
        p10l = f * (1.0 - self.tertiary_foreign)
        p00h = (1.0 - f) * self.tertiary_native
        p00l = (1.0 - f) * (1.0 - self.tertiary_native)
        # symmetric equilibrium: the born population by skill equals the resident one
        return np.array([p00l, p10l, p00h, p10h]), (p00l + p10l, p00h + p10h)


@dataclass(frozen=True)
class Calibration:
    S: float
    G: float
    lam: float
    kappa_l: float
    kappa_h: float
    alloc: tuple          # [π00l, π10l, π00h, π10h]
    origin_low: float
    origin_high: float


def _Nia(p_nat, p_imm, eps, phi):
    return (p_nat**((eps-1)/eps) + phi * p_imm**((eps-1)/eps))**(eps/(eps-1))


def symmetric_wages(p, S, rho, eps, phi, gamma=0.0):
    """
    Wages [w00l, w10l, w00h, w10h], output F, and the land return r at allocation p
    (A = L = 1; F = N^{1−γ}); as model/equations.py.
    """
    p00l, p10l, p00h, p10h = p
    Nl, Nh = _Nia(p00l, p10l, eps, phi), _Nia(p00h, p10h, eps, phi)
    N = (Nl**((rho-1)/rho) + (S*Nh)**((rho-1)/rho))**(rho/(rho-1))
    F = N**(1.0 - gamma)
    w00l = ((1-gamma) * F * Nl**((rho-1)/rho - (eps-1)/eps)
            * N**((1-rho)/rho) * p00l**(-1/eps))
    w00h = ((1-gamma) * F * S**((rho-1)/rho) * Nh**((rho-1)/rho - (eps-1)/eps)
            * N**((1-rho)/rho) * p00h**(-1/eps))
    w10l = phi * (p10l/p00l)**(-1/eps) * w00l
    w10h = phi * (p10h/p00h)**(-1/eps) * w00h
    r = gamma * N**(1.0 - gamma)
    return np.array([w00l, w10l, w00h, w10h]), F, r


def immigrant_wage_ratio(p, eps, phi):
    """Mean over skills of the immigrant/native wage ratio, w_{10}^a/w_{00}^a =
    φ (π_{10}^a/π_{00}^a)^{−1/ε} (it does not depend on S, ρ or γ)."""
    return 0.5 * phi * ((p[1] / p[0])**(-1.0/eps) + (p[3] / p[2])**(-1.0/eps))


def phi_for_wage_ratio(p, eps, ratio):
    """φ at which the mean immigrant/native wage ratio equals `ratio`."""
    return ratio / immigrant_wage_ratio(p, eps, 1.0)


def skill_premium_S(p, premium, rho, eps, phi):
    """S such that S^{(ρ−1)/ρ}(N^l/N^h)^{1/ρ} equals the target premium."""
    Nl, Nh = _Nia(p[0], p[1], eps, phi), _Nia(p[2], p[3], eps, phi)
    return (premium * (Nh/Nl)**(1/rho))**(rho/(rho-1))


def budget_lambda(p, w, G, tau_prog, rev_capital=0.0):
    """λ clearing Σ π (w − λ w^{1−τ}) + R_k = G."""
    return (np.dot(p, w) + rev_capital - G) / np.dot(p, w**(1.0 - tau_prog))


def kappa_from_stay(p, w, lam, tau_prog, origin_low, origin_high, sigma,
                    land_income=(0.0, 0.0)):
    """
    κ^a = (c_{j,-j} + C^a) / (V_stay [(1−P^a)/P^a]^{1/σ}), V_stay = c_{jj} + C^a,
    c = λ w^{1−τ}, C^a origin land income; κ^a = (c_{j,-j}/c_{jj})[P^a/(1−P^a)]^{1/σ}
    when C^a = 0.
    """
    c = lam * w**(1.0 - tau_prog)
    out = []
    for i_nat, i_imm, mass, C in ((0, 1, origin_low, land_income[0]),
                                  (2, 3, origin_high, land_income[1])):
        P = p[i_nat] / mass
        V_move = (c[i_nat] + C) * ((1.0 - P) / P)**(1.0/sigma)
        out.append((c[i_imm] + C) / V_move)
    return tuple(out)


def calibrate(targets: Targets, sigma, rho, eps, phi, tau_prog,
              gamma=0.0, tau_k=0.0, ell=(0.0, 0.0)) -> Calibration:
    """
    Full calibration (S, G, λ, κ) of the specification (σ, ρ, ε, φ, τ, γ) to the data
    targets. ell = (ℓ^l, ℓ^h): land endowments per origin worker by skill.
    """
    p, (Pi_l, Pi_h) = targets.allocation()
    # the premium S^{(ρ−1)/ρ}(N^l/N^h)^{1/ρ} is a ratio of marginal products of
    # N^h and N^l: the factor (1−γ)F/N cancels, so γ does not enter step 2
    S = skill_premium_S(p, targets.skill_premium, rho, eps, phi)
    w, F, r = symmetric_wages(p, S, rho, eps, phi, gamma)
    G = targets.g_share * F
    rev_k = tau_k * r * (ell[0] * Pi_l + ell[1] * Pi_h)
    lam = budget_lambda(p, w, G, tau_prog, rev_k)
    C = (r * ell[0] * (1.0 - tau_k), r * ell[1] * (1.0 - tau_k))
    k_l, k_h = kappa_from_stay(p, w, lam, tau_prog, Pi_l, Pi_h, sigma, C)
    return Calibration(S=float(S), G=float(G), lam=float(lam),
                       kappa_l=float(k_l), kappa_h=float(k_h),
                       alloc=tuple(float(x) for x in p),
                       origin_low=float(Pi_l), origin_high=float(Pi_h))
