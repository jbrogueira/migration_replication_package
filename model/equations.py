# model/equations.py

from config import σ, ϵ, φ, S, ρ, γ, A0, A1, L0, L1, tau_k, ell_l, ell_h
from .parameters import κl, κh


# ------------------------- Core functions -------------------------
def Nia(piia, pija):
    return (piia**((ϵ-1)/ϵ) + φ * pija**((ϵ-1)/ϵ))**(ϵ/(ϵ-1))

def Ni(pi):
    piil, pijl, piih, pijh = pi
    return (Nia(piil,pijl)**((ρ-1)/ρ) + (S * Nia(piih,pijh))**((ρ-1)/ρ))**(ρ/(ρ-1))

def Fi(i, pi):
    A = A0 if i==0 else A1
    L = L0 if i==0 else L1
    return A * L**γ * Ni(pi)**(1-γ)

def wiil(i, pi):
    piil, pjil, piih, pjih = pi
    return (1-γ) * Fi(i,pi) * Nia(piil,pjil)**((ρ-1)/ρ-(ϵ-1)/ϵ) * Ni(pi)**((1-ρ)/ρ) * piil**(-1/ϵ)

def wjil(i, pi):
    piil, pjil, *_ = pi
    return φ * (pjil/piil)**(-1/ϵ) * wiil(i,pi)

def wiih(i, pi):
    piil,pjil,piih,pjih = pi
    return (1-γ) * Fi(i,pi) * S**((ρ-1)/ρ) * Nia(piih,pjih)**((ρ-1)/ρ-(ϵ-1)/ϵ) * Ni(pi)**((1-ρ)/ρ) * piih**(-1/ϵ)

def wjih(i, pi):
    *_,piih,pjih = pi
    return φ * (pjih/piih)**(-1/ϵ) * wiih(i,pi)

def wage(j, i, a, pi):
    native = wiil if a=='l' and j==i else wiih if a=='h' and j==i else None
    foreign = wjil if a=='l' and j!=i else wjih if a=='h' and j!=i else None
    return native(i,pi) if j==i else foreign(i,pi)

def ri(i, pi):
    Ai = A0 if i==0 else A1
    Li = L0 if i==0 else L1
    return Ai * γ * Li**(γ-1) * Ni(pi)**(1-γ)

def taul(y, λ):
    from config import τl
    return 1 - λ * y**(-τl)

# ---------------- Migration probability
def cap_tax_revenue(p_res, r_own: float, r_other: float) -> float:
    """
    Residence-based tax on land income (app:land): country i taxes the land income of
    its residents — natives on their own country's land (rent r_own), immigrants on
    their origin country's land (rent r_other). p_res = [natives l, immigrants l,
    natives h, immigrants h] resident in i.
    """
    return tau_k * ((p_res[0]*ell_l + p_res[2]*ell_h) * r_own
                    + (p_res[1]*ell_l + p_res[3]*ell_h) * r_other)


def land_income_origin(j: int, a: str, pj) -> float:
    ell = ell_l if a == 'l' else ell_h
    return ri(j, pj) * ell * (1.0 - tau_k)

def Pi(j, i, a, pj, p_j, τj, τ_j):
    """
    Migration probability with origin-based land income inside powered terms.

    Args:
        j: origin country (0/1)
        i: destination country (0/1)
        a: 'l' or 'h' (skill type)
        pj: composition in origin j
        p_j: composition in destination i
        τj: taxes in j (native), ordering [τ00l, τ10l, τ00h, τ10h]
        τ_j: taxes in i (immigrant), same ordering
    Returns:
        P(j→i | a)
    """
    if j != i:
        κ_a   = κl if a == 'l' else κh
        tau_i = τ_j[1] if a == 'l' else τ_j[3]   # immigrant tax in i
        tau_j0= τj [0] if a == 'l' else τj [2]   # native tax in j

        # Net wage components
        w_move_net = wage(j, i, a, p_j) * (1.0 - tau_i)
        w_stay_net = wage(j, j, a, pj ) * (1.0 - tau_j0)

        # Origin land income (same if move or stay): C_{j,a} = r_j(pj) * ell_a * (1 - tau_k)
        ell = ell_l if a == 'l' else ell_h
        C   = ri(j, pj) * ell * (1.0 - tau_k)

        # Power-on-levels discrete choice; κ divides consumption, land income included
        eps    = 1e-12
        V_move = max((w_move_net + C) / κ_a, eps)
        V_stay = max(w_stay_net + C, eps)

        num = V_move ** σ
        den = num + V_stay ** σ
        return num / den
    else:
        # staying probability = 1 - prob(move to 1-j)
        return 1.0 - Pi(j, 1 - j, a, pj, p_j, τj, τ_j)
