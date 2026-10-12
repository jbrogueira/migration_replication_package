import numpy as np
import pandas as pd
from scipy import optimize
from .equations import wage, taul, Pi, ri, cap_tax_revenue       
from config import τl as TAU_L_BASE, G0, λ0_vals, bounds_lambda, G1, tau_k, ell_l, ell_h, x0, bnds, origin0_low, origin0_high, bounds_lambda, bounds_lambda_8
import inspect

# --- Symmetric GE equilibrium solver (baseline) ---

def err_ge_eq(params, out=1, custom_λ0=None):
    p0 = params[0:4]
    λ0 = params[4] if custom_λ0 is None else custom_λ0

    # home wages & taxes (country 0 perspective)
    w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
    w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)
    τ00l = taul(w00l, λ0);   τ10l = taul(w10l, λ0)
    τ00h = taul(w00h, λ0);   τ10h = taul(w10h, λ0)
    τ0   = np.array([τ00l, τ10l, τ00h, τ10h])

    # Government revenue in country 0:
    #  (i) labor-tax on residents in 0 (natives + immigrants)
    rev_labor0 = (
        p0[0]*τ00l*w00l + p0[1]*τ10l*w10l +
        p0[2]*τ00h*w00h + p0[3]*τ10h*w10h
    )
    #  (ii) residence-based tax on land income; symmetric world: the immigrants' origin
    #       rent equals the home rent
    r0 = ri(0, p0)
    rev_cap0 = cap_tax_revenue(p0, r0, r0)

    err = [
        rev_labor0 + rev_cap0 - G0,                               # govt budget (0)
        p0[0] - Pi(0,0,'l', p0, p0, τ0, τ0)*origin0_low,          # natives low in 0
        p0[2] - Pi(0,0,'h', p0, p0, τ0, τ0)*origin0_high,         # natives high in 0
        p0[0] + p0[1] - origin0_low,                              # residents low in 0
        p0[2] + p0[3] - origin0_high                              # residents high in 0
    ]
    return np.sum(np.square(err)) if out == 1 else err


# run this once at import to get pi_ge_eq
_sym = optimize.minimize(err_ge_eq, x0, method='SLSQP', bounds=bnds, tol=1e-12)
assert _sym.success, f"Baseline equilibrium solver failed: {_sym.message}"
assert _sym.fun < 1e-8, f"Baseline equilibrium residual too large: {_sym.fun:.2e}"
pi_ge_eq = _sym.x
# polish on the residual vector (square system: 5 conditions, 5 unknowns); SLSQP on
# the sum of squares stops at residuals ~1e-6
_pol = optimize.root(lambda v: err_ge_eq(v, out=0), pi_ge_eq, method='hybr', tol=1e-14)
if np.max(np.abs(err_ge_eq(_pol.x, out=0))) < np.max(np.abs(err_ge_eq(pi_ge_eq, out=0))):
    pi_ge_eq = _pol.x

# --- Laffer‐type residual for symmetric sweep (baseline) ---
def err_ge_eq_laffer(params, out=1, custom_λ0=0.95):
    p0 = params[0:4]
    λ0 = custom_λ0

    # home wages & taxes
    w = [
        wage(0,0,'l',p0),
        wage(1,0,'l',p0),
        wage(0,0,'h',p0),
        wage(1,0,'h',p0)
    ]
    τ0 = np.array([taul(wi, λ0) for wi in w])

    # only the population-share residuals, allow bc to change
    res = [
        p0[0] - Pi(0,0,'l', p0, p0, τ0, τ0)*origin0_low,
        p0[2] - Pi(0,0,'h', p0, p0, τ0, τ0)*origin0_high,
        p0[0] + p0[1] - origin0_low,
        p0[2] + p0[3] - origin0_high
    ]
    return np.sum(np.square(res)) if out==1 else res

# --- Precompute foreign equilibrium values ---
p_foreign    = pi_ge_eq[:4]
λ_star       = pi_ge_eq[4]
tau_star     = 1 - λ_star
w10h_star    = wage(1, 0, 'h', pi_ge_eq[:4])  
tau10h_star  = taul(w10h_star, λ_star)

# foreign wages and taxes
foreign_w = { (i,a): wage(i,0,a,p_foreign)
    for i,a in [(0,'l'),(1,'l'),(0,'h'),(1,'h')] }
τ1 = np.array([taul(foreign_w[(i,a)], λ_star)
    for i,a in [(0,'l'),(1,'l'),(0,'h'),(1,'h')]])

# --- SOE residual under Laffer tax ---
def err_ge_soe_laffer(params, out=1, custom_λ0=0.95):
    p0 = params[:4]
    λ0 = custom_λ0

    home_w = {
        (0, 'l'): wage(0, 0, 'l', p0),
        (1, 'l'): wage(1, 0, 'l', p0),
        (0, 'h'): wage(0, 0, 'h', p0),
        (1, 'h'): wage(1, 0, 'h', p0),
    }

    τ0 = np.array([
        taul(home_w[(0, 'l')], λ0),
        taul(home_w[(1, 'l')], λ0),
        taul(home_w[(0, 'h')], λ0),
        taul(home_w[(1, 'h')], λ0),
    ])

    res = [
        # native stay-share conditions
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * origin0_low,
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * origin0_high,

        # immigrant no-arbitrage on wages
        (1 - τ0[1]) * home_w[(1, 'l')] - (1 - τ1[1]) * foreign_w[(1, 'l')],
        (1 - τ0[3]) * home_w[(1, 'h')] - (1 - τ1[3]) * foreign_w[(1, 'h')],
    ]

    return np.sum(np.square(res)) if out == 1 else res


# --- Prefixed high-skill tax residual ---
def err_ge_pref_tax(params, out=1, custom_λ0=None, custom_τ10h=None):
    p0 = params

    home_w = {
        (0, 'l'): wage(0, 0, 'l', p0),
        (1, 'l'): wage(1, 0, 'l', p0),
        (0, 'h'): wage(0, 0, 'h', p0),
        (1, 'h'): wage(1, 0, 'h', p0),
    }

    τ0 = np.array([
        taul(home_w[(0, 'l')], custom_λ0),
        taul(home_w[(1, 'l')], custom_λ0),
        taul(home_w[(0, 'h')], custom_λ0),
        custom_τ10h
    ])

    res = [
        # native stay-share conditions
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * origin0_low,
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * origin0_high,

        # immigrant no-arbitrage on wages only 
        (1 - τ0[1]) * home_w[(1, 'l')] - (1 - τ1[1]) * foreign_w[(1, 'l')],
        (1 - τ0[3]) * home_w[(1, 'h')] - (1 - τ1[3]) * foreign_w[(1, 'h')],
    ]
    return np.sum(np.square(res)) if out == 1 else res



# --- Sweep helper ---
def sweep_equilibrium(err_func, λ_vals, bounds, x0=None, symmetric=False):
    """
    Sweep equilibrium over λ_vals for the given residual function.
    Uses x0[:4] if provided, else globals()['x0'][:4].
    τ10h can be pinned via err_func.tau10h_fixed (float).
    Land income is always origin-based: r0 = ri(0,p0), r1 = ri(1,p_foreign).
    """
    # init (simple)
    base_x0 = np.asarray(x0 if x0 is not None else globals()['x0'])
    assert base_x0.size >= 4, "x0 must have at least 4 elements."
    curr_x = base_x0[:4].astype(float)

    tau10h_fixed = getattr(err_func, "tau10h_fixed", None)

    cols = [
        "λ0",
        "p00l","p10l","p00h","p10h",
        "w00l","w10l","w00h","w10h",
        "tau00l","tau10l","tau00h","tau10h",
        "inc00l","inc10l","inc00h","inc10h",
        "TR","fun"
    ]
    out = {c: [] for c in cols}

    for λ in λ_vals:
        sol = optimize.minimize(
            lambda v: err_func(v, out=1, custom_λ0=λ),
            x0=curr_x, method="SLSQP", bounds=bounds,
            tol=1e-12, options={"ftol": 1e-12, "maxiter": 1000}
        )
        if not sol.success or sol.fun > 1e-8:
            continue

        p0 = sol.x[:4]

        w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
        w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)

        # taxes
        tau00l = taul(w00l, λ)
        tau10l = taul(w10l, λ)
        tau00h = taul(w00h, λ)
        tau10h = (float(tau10h_fixed) if tau10h_fixed is not None else taul(w10h, λ))

        # origin-based land returns (immigrants' origin: the symmetric world in the
        # multilateral sweep, the frozen foreign allocation otherwise)
        r0 = ri(0, p0)
        r1 = r0 if symmetric else ri(1, p_foreign)

        land00l = r0 * ell_l * (1.0 - tau_k)
        land10l = r1 * ell_l * (1.0 - tau_k)
        land00h = r0 * ell_h * (1.0 - tau_k)
        land10h = r1 * ell_h * (1.0 - tau_k)

        # disposable income
        inc00l = w00l * (1 - tau00l) + land00l
        inc10l = w10l * (1 - tau10l) + land10l
        inc00h = w00h * (1 - tau00h) + land00h
        inc10h = w10h * (1 - tau10h) + land10h

        # gov revenue in country 0
        rev_labor0 = (
            p0[0]*tau00l*w00l + p0[1]*tau10l*w10l +
            p0[2]*tau00h*w00h + p0[3]*tau10h*w10h
        )
        rev_cap0 = cap_tax_revenue(p0, r0, r1)
        TR = rev_labor0 + rev_cap0

        row = [
            λ,
            p0[0], p0[1], p0[2], p0[3],
            w00l, w10l, w00h, w10h,
            tau00l, tau10l, tau00h, tau10h,
            inc00l, inc10l, inc00h, inc10h,
            TR, sol.fun
        ]
        for c, v in zip(cols, row):
            out[c].append(v)

        curr_x = p0.copy()

    return pd.DataFrame(out)




###############################################################################
################### LARGE COUNTRY EQUILIBRIA ##################################
###############################################################################
# --- Two-country residual (PE/GE) ---
def err_2country(
    params,
    out=1,
    λ0=None,
    λ1=None,
    tau_high_imm_override=None,   # GLOBAL: sets τ10h (C0) and τ01h (C1)
    tau10h_override=None,         # UNILATERAL: sets τ10h (C0) only if GLOBAL is None
    mode='GE',                    # 'GE' → λ1 from params[-1] and adds C1 GBC; 'PE' → λ1 from arg
    g_per_capita=False            # C1 spending G1 times its resident population (1 at baseline)
):
    """Two-country residual. Country 0's budget is not in this system: callers clear it
    outside (λ0 sweep + root-find), against G0 or, with g_per_capita, against G0 times
    country 0's resident population. g_per_capita here applies to country 1 only."""
    if λ0 is None:
        raise ValueError("Pass λ0 via argument.")

    p0 = params[0:4]   # [p00l, p10l, p00h, p10h]
    p1 = params[4:8]   # [p11l, p01l, p11h, p01h]

    λ0_val = λ0
    λ1_val = (params[-1] if mode == 'GE' else λ1)

    # wages
    w00l = wage(0,0,'l',p0);  w10l = wage(1,0,'l',p0)
    w00h = wage(0,0,'h',p0);  w10h = wage(1,0,'h',p0)
    w11l = wage(1,1,'l',p1);  w01l = wage(0,1,'l',p1)
    w11h = wage(1,1,'h',p1);  w01h = wage(0,1,'h',p1)

    # taxes C0
    τ00l = taul(w00l, λ0_val)
    τ10l = taul(w10l, λ0_val)
    τ00h = taul(w00h, λ0_val)
    τ10h = (float(tau_high_imm_override) if tau_high_imm_override is not None
            else float(tau10h_override)   if tau10h_override      is not None
            else taul(w10h, λ0_val))
    τ0 = np.array([τ00l, τ10l, τ00h, τ10h])

    # taxes C1
    τ11l = taul(w11l, λ1_val)
    τ01l = taul(w01l, λ1_val)
    τ11h = taul(w11h, λ1_val)
    τ01h = (float(tau_high_imm_override) if tau_high_imm_override is not None
            else taul(w01h, λ1_val))
    τ1 = np.array([τ11l, τ01l, τ11h, τ01h])

    # residuals: stay shares + accounting identities
    err = [
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1) * origin0_low,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1) * origin0_high,
        p0[0] + p1[1] - origin0_low,
        p0[2] + p1[3] - origin0_high,

        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0) * origin0_low,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0) * origin0_high,
        p1[0] + p0[1] - origin0_low,
        p1[2] + p0[3] - origin0_high,
    ]

    if mode == 'GE':
        TR1_lab = (p1[0]*τ11l*w11l + p1[1]*τ01l*w01l + p1[2]*τ11h*w11h + p1[3]*τ01h*w01h)
        TR1_cap = cap_tax_revenue(p1, ri(1, p1), ri(0, p0))   # residence-based
        G1_req = G1 * float(np.sum(p1)) if g_per_capita else G1
        err.insert(0, TR1_lab + TR1_cap - G1_req)

    return np.sum(np.square(err)) if out == 1 else err


def _polish(res_vec, x):
    """Newton polish of an SLSQP solution on the residual vector of a square system;
    kept only if it lowers the largest residual."""
    pol = optimize.root(res_vec, x, method='hybr', tol=1e-14)
    size = lambda v: np.max(np.abs(res_vec(v)))
    return pol.x if np.all(np.isfinite(pol.x)) and size(pol.x) < size(x) else x


# --- Two-country sweep (PE/GE) ---
def sweep_equilibrium_2country(
    err_func,
    λ0_vals,
    λ1_fixed,
    bounds,
    x0=None,
    free_lambda1=False,
    **report_kwargs  # optional, only used for *reporting* τ’s
):
    # allow λ1 to float
    if free_lambda1:
        bounds = bounds + ((1e-5, 1.5),)

    # default x0
    if x0 is None:
        p0_star = pi_ge_eq[:4]
        base = np.concatenate([p0_star, p0_star])
        x0 = base if not free_lambda1 else np.concatenate([base, [pi_ge_eq[4]]])

    cols = [
        "λ0","λ1",
        "p00l","p10l","p00h","p10h","p11l","p01l","p11h","p01h",
        "w00l","w10l","w00h","w10h","w11l","w01l","w11h","w01h",
        "tau00l","tau10l","tau00h","tau10h","tau11l","tau01l","tau11h","tau01h",
        "inc00l","inc10l","inc00h","inc10h","inc11l","inc01l","inc11h","inc01h",
        "TR0","TR1","fun"
    ]
    results = {c: [] for c in cols}

    # reporting-only overrides (don’t pass into err_func)
    tau_hi_imm_global = report_kwargs.get("tau_high_imm_override", None)
    tau10h_unilateral = report_kwargs.get("tau10h_override", report_kwargs.get("custom_τ10h", None))
    pick_tau = lambda base, g=None, u=None: (float(g) if g is not None else (float(u) if u is not None else base))

    for λ0 in λ0_vals:
        # call residual with *positional* args only: (params, out, λ0, λ1)
        sol = optimize.minimize(
            lambda params: err_func(params, 1, λ0, (None if free_lambda1 else λ1_fixed)),
            x0, method="SLSQP", bounds=bounds, tol=1e-10
        )
        assert sol.success, f"Minimization failed at λ0={λ0}: {sol.message}"

        p = _polish(lambda v: np.asarray(err_func(v, 0, λ0, (None if free_lambda1 else λ1_fixed))),
                    sol.x)
        sol.x = p
        p0, p1 = p[:4], p[4:8]
        λ1_val = sol.x[8] if free_lambda1 else λ1_fixed

        # wages
        w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
        w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])

        # taxes for reporting (apply reporting override if provided)
        τ0 = np.array([
            taul(w0[0], λ0),
            taul(w0[1], λ0),
            taul(w0[2], λ0),
            pick_tau(taul(w0[3], λ0), g=tau_hi_imm_global, u=tau10h_unilateral)
        ])
        τ1 = np.array([
            taul(w1[0], λ1_val),
            taul(w1[1], λ1_val),
            taul(w1[2], λ1_val),
            pick_tau(taul(w1[3], λ1_val), g=tau_hi_imm_global, u=None)
        ])

        # origin-based land returns
        r0, r1 = ri(0, p0), ri(1, p1)
        land00l = r0 * ell_l * (1 - tau_k); land10l = r1 * ell_l * (1 - tau_k)
        land00h = r0 * ell_h * (1 - tau_k); land10h = r1 * ell_h * (1 - tau_k)
        land11l = r1 * ell_l * (1 - tau_k); land01l = r0 * ell_l * (1 - tau_k)
        land11h = r1 * ell_h * (1 - tau_k); land01h = r0 * ell_h * (1 - tau_k)

        inc0 = np.array([
            w0[0]*(1-τ0[0]) + land00l,
            w0[1]*(1-τ0[1]) + land10l,
            w0[2]*(1-τ0[2]) + land00h,
            w0[3]*(1-τ0[3]) + land10h
        ])
        inc1 = np.array([
            w1[0]*(1-τ1[0]) + land11l,
            w1[1]*(1-τ1[1]) + land01l,
            w1[2]*(1-τ1[2]) + land11h,
            w1[3]*(1-τ1[3]) + land01h
        ])

        TR0 = float(np.dot(p0, τ0 * w0)) + cap_tax_revenue(p0, r0, r1)
        TR1 = float(np.dot(p1, τ1 * w1)) + cap_tax_revenue(p1, r1, r0)

        row = [λ0, λ1_val, *p0, *p1, *w0, *w1, *τ0, *τ1, *inc0, *inc1, TR0, TR1, sol.fun]
        for c, v in zip(cols, row):
            results[c].append(v)

        x0 = sol.x  # warm start

    return pd.DataFrame(results)



# --- (PE) Residual with preferential τ10h in C0 ---
def err_pe_large_tau10h(params, out=1, custom_λ0=None, custom_λ1=None, custom_τ10h=None):
    p0 = params[0:4]
    p1 = params[4:8]

    w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
    w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])

    τ0 = np.array([taul(w0[0], custom_λ0), taul(w0[1], custom_λ0), taul(w0[2], custom_λ0), custom_τ10h])
    τ1 = np.array([taul(w1[0], custom_λ1), taul(w1[1], custom_λ1), taul(w1[2], custom_λ1), taul(w1[3], custom_λ1)])

    err = [
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1) * origin0_low,
        p0[0] + p1[1] - origin0_low,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1) * origin0_high,
        p0[2] + p1[3] - origin0_high,

        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0) * origin0_low,
        p1[0] + p0[1] - origin0_low,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0) * origin0_high,
        p1[2] + p0[3] - origin0_high,
    ]
    return np.sum(np.square(err)) if out == 1 else err


# --- GE global preferential tax glue ---
def ge_pref_residual_unilateral(params, out, custom_λ0, custom_λ1=None, custom_τ10h=None,
                                g_per_capita=False):
    """Unilateral preferential tax in country 0; country 1 has no preferential tax and
    λ1 (params[-1]) clears its budget."""
    return err_2country(
        params=params,
        out=out,
        λ0=custom_λ0,
        λ1=None,                          # GE: λ1 floats in params[-1]
        tau10h_override=custom_τ10h,      # preferential rate in country 0 only
        mode='GE',
        g_per_capita=g_per_capita
    )


def ge_pref_residual_global(params, out, custom_λ0, custom_λ1=None, custom_τ10h=None,
                            g_per_capita=False):
    return err_2country(
        params=params,
        out=out,
        λ0=custom_λ0,
        λ1=None,                          # GE: λ1 floats in params[-1]
        tau_high_imm_override=custom_τ10h,
        mode='GE',
        g_per_capita=g_per_capita
    )


###############################################################################
############### CHANGING PROGRESSIVITY PARAMETER ##############################
###############################################################################

# --- τ_l helpers (progressivity sweeps with λ fixed) -------------------------
def _taul_param(y, lam, tau_l_param):
    return 1.0 - lam * y**(-tau_l_param)


def err_ge_soe_pref_tau(params, out=1, lam=None, tau_l_param=None, τ10h_pref=None):
    """
    SOE residual with preferential τ10h (high-skill immigrants in C0).
    No GBC here; land only matters inside Pi via incentives.
    """
    p0 = params[:4]

    # wages (home, from country-0 perspective)
    w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
    w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)

    # taxes (progressive schedule for all except τ10h, which is pinned)
    τ00l = _taul_param(w00l, lam, tau_l_param)
    τ10l = _taul_param(w10l, lam, tau_l_param)
    τ00h = _taul_param(w00h, lam, tau_l_param)
    τ10h = τ10h_pref
    τ0 = np.array([τ00l, τ10l, τ00h, τ10h])

    # residuals (stay-shares + immigrant no-arbitrage)
    res = [
        p0[0] - Pi(0,0,'l', p0, p_foreign, τ0, τ1) * origin0_low,
        p0[2] - Pi(0,0,'h', p0, p_foreign, τ0, τ1) * origin0_high,
        (1 - τ10l) * w10l - (1 - τ1[1]) * foreign_w[(1,'l')],
        (1 - τ10h) * w10h - (1 - τ1[3]) * foreign_w[(1,'h')],
    ]
    return np.sum(np.square(res)) if out == 1 else res


def sweep_taul_soe_pref(lambda_star, tau_grid, τ10h_pref):
    # initial guess and bounds
    curr_x = np.asarray(x0)[:4].astype(float)
    bounds = tuple(bounds_lambda[:4])

    # foreign land return is fixed in SOE
    r1_foreign = ri(1, p_foreign)

    rows = []
    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda x: err_ge_soe_pref_tau(x, out=1, lam=lambda_star,
                                          tau_l_param=tau_l_param, τ10h_pref=τ10h_pref),
            x0=curr_x, method='SLSQP', bounds=bounds,
            tol=1e-10, options={"ftol": 1e-12, "maxiter": 1000}
        )
        if not sol.success:
            continue  # keep your current behavior; switch to `assert sol.success, ...` if you prefer hard-fail

        p0 = sol.x[:4]
        curr_x = p0.copy()  # warm start

        # wages
        w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
        w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)

        # taxes at this (lambda_star, tau_l_param)
        tau00l = _taul_param(w00l, lambda_star, tau_l_param)
        tau10l = _taul_param(w10l, lambda_star, tau_l_param)
        tau00h = _taul_param(w00h, lambda_star, tau_l_param)
        tau10h = τ10h_pref

        # government revenue in country 0: labor + origin-based land
        TR_lab = (
            p0[0]*tau00l*w00l + p0[1]*tau10l*w10l +
            p0[2]*tau00h*w00h + p0[3]*tau10h*w10h
        )
        r0 = ri(0, p0)
        TR_cap = cap_tax_revenue(p0, r0, r1_foreign)
        TR = TR_lab + TR_cap

        # disposable income (labor net + origin-based net land income)
        land00l = r0 * ell_l * (1 - tau_k)
        land00h = r0 * ell_h * (1 - tau_k)
        land10l = r1_foreign * ell_l * (1 - tau_k)
        land10h = r1_foreign * ell_h * (1 - tau_k)

        inc00l = w00l*(1 - tau00l) + land00l
        inc10l = w10l*(1 - tau10l) + land10l
        inc00h = w00h*(1 - tau00h) + land00h
        inc10h = w10h*(1 - tau10h) + land10h

        rows.append(dict(
            tau_l=tau_l_param, TR=TR,
            p00l=p0[0], p10l=p0[1], p00h=p0[2], p10h=p0[3],
            w00l=w00l, w10l=w10l, w00h=w00h, w10h=w10h,
            tau00l=tau00l, tau10l=tau10l, tau00h=tau00h, tau10h=tau10h,
            inc00l=inc00l, inc10l=inc10l, inc00h=inc00h, inc10h=inc10h,
            fun=sol.fun
        ))

    return pd.DataFrame(rows).sort_values('tau_l').reset_index(drop=True)





def _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lam_star, fun_val,
                           tau_l_param1=None):
    # origin-based land returns
    r0 = ri(0, p0)
    r1 = ri(1, p1)

    # net land per capita (origin-based)
    land00l = r0 * ell_l * (1 - tau_k); land00h = r0 * ell_h * (1 - tau_k)
    land10l = r1 * ell_l * (1 - tau_k); land10h = r1 * ell_h * (1 - tau_k)
    land11l = r1 * ell_l * (1 - tau_k); land11h = r1 * ell_h * (1 - tau_k)
    land01l = r0 * ell_l * (1 - tau_k); land01h = r0 * ell_h * (1 - tau_k)

    # after-tax incomes = net labor + net land (origin-based)
    inc0 = np.array([
        w0[0]*(1-τ0[0]) + land00l,  # inc00l
        w0[1]*(1-τ0[1]) + land10l,  # inc10l
        w0[2]*(1-τ0[2]) + land00h,  # inc00h
        w0[3]*(1-τ0[3]) + land10h,  # inc10h
    ])
    inc1 = np.array([
        w1[0]*(1-τ1[0]) + land11l,  # inc11l
        w1[1]*(1-τ1[1]) + land01l,  # inc01l
        w1[2]*(1-τ1[2]) + land11h,  # inc11h
        w1[3]*(1-τ1[3]) + land01h,  # inc01h
    ])

    # tax revenue = labor + origin-based land
    TR0 = float(np.dot(p0, τ0 * w0)) + cap_tax_revenue(p0, r0, r1)
    TR1 = float(np.dot(p1, τ1 * w1)) + cap_tax_revenue(p1, r1, r0)

    return {
        'tauL0': tau_l_param,
        'tauL1': tau_l_param if tau_l_param1 is None else tau_l_param1,
        # populations
        'p00l': p0[0], 'p10l': p0[1], 'p00h': p0[2], 'p10h': p0[3],
        'p11l': p1[0], 'p01l': p1[1], 'p11h': p1[2], 'p01h': p1[3],
        # wages
        'w00l': w0[0], 'w10l': w0[1], 'w00h': w0[2], 'w10h': w0[3],
        'w11l': w1[0], 'w01l': w1[1], 'w11h': w1[2], 'w01h': w1[3],
        # taxes
        'tau00l': τ0[0], 'tau10l': τ0[1], 'tau00h': τ0[2], 'tau10h': τ0[3],
        'tau11l': τ1[0], 'tau01l': τ1[1], 'tau11h': τ1[2], 'tau01h': τ1[3],
        # incomes
        'inc00l': inc0[0], 'inc10l': inc0[1], 'inc00h': inc0[2], 'inc10h': inc0[3],
        'inc11l': inc1[0], 'inc01l': inc1[1], 'inc11h': inc1[2], 'inc01h': inc1[3],
        # revenues (mirror into G0/G1 for convenience)
        'TR0': TR0, 'TR1': TR1, 'G0': TR0, 'G1': TR1,
        'λ0': lam_star, 'λ1': lam_star, 'fun': fun_val
    }



def _err_two_country_case2(params, out, lam_star, tau_l_param, τ10h_pref):
    p0 = params[:4]; p1 = params[4:8]
    w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
    w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])

    τ0 = np.array([
        _taul_param(w0[0], lam_star, tau_l_param),
        _taul_param(w0[1], lam_star, tau_l_param),
        _taul_param(w0[2], lam_star, tau_l_param),
        τ10h_pref
    ])
    # country 1 has no preferential tax; its progressivity τl1 = params[8] clears its
    # budget (λ* in both countries)
    tau_l1 = params[8]
    τ1 = np.array([
        _taul_param(w1[0], lam_star, tau_l1),
        _taul_param(w1[1], lam_star, tau_l1),
        _taul_param(w1[2], lam_star, tau_l1),
        _taul_param(w1[3], lam_star, tau_l1)
    ])
    TR1 = (float(np.dot(p1, τ1 * w1))
           + cap_tax_revenue(p1, ri(1, p1), ri(0, p0)))

    err = [
        TR1 - G1,
        # C0
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1)*origin0_low,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1)*origin0_high,
        p0[0] + p1[1] - origin0_low,
        p0[2] + p1[3] - origin0_high,
        # C1
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0)*origin0_low,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0)*origin0_high,
        p1[0] + p0[1] - origin0_low,
        p1[2] + p0[3] - origin0_high,
    ]
    return np.sum(np.square(err)) if out == 1 else err



def _err_two_country_case3(params, out, lam_star, tau_l_param, τ10h_pref, τ01h_pref):
    p0 = params[:4]; p1 = params[4:8]
    w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
    w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])

    τ0 = np.array([
        _taul_param(w0[0], lam_star, tau_l_param),
        _taul_param(w0[1], lam_star, tau_l_param),
        _taul_param(w0[2], lam_star, tau_l_param),
        τ10h_pref
    ])
    τ1 = np.array([
        _taul_param(w1[0], lam_star, tau_l_param),
        _taul_param(w1[1], lam_star, tau_l_param),
        _taul_param(w1[2], lam_star, tau_l_param),
        τ01h_pref
    ])

    err = [
        # C0
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1)*origin0_low,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1)*origin0_high,
        p0[0] + p1[1] - origin0_low,
        p0[2] + p1[3] - origin0_high,
        # C1
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0)*origin0_low,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0)*origin0_high,
        p1[0] + p0[1] - origin0_low,
        p1[2] + p0[3] - origin0_high,
    ]
    return np.sum(np.square(err)) if out == 1 else err



def sweep_taul_large_case2(lambda_star, tau_grid, τ10h_pref, x0=None):
    """Country 0: preferential tax, progressivity τl0 swept; country 1: no preferential
    tax, progressivity τl1 (x[8]) clears its budget."""
    if x0 is None:
        pstar = pi_ge_eq[:4]; x0 = np.concatenate([pstar, pstar, [TAU_L_BASE]])
    bounds = tuple(bounds_lambda_8) + ((0.0, 1.0),)

    rows = []
    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda v: _err_two_country_case2(v, out=1, lam_star=lambda_star,
                                             tau_l_param=tau_l_param, τ10h_pref=τ10h_pref),
            x0, method='SLSQP', bounds=bounds, tol=1e-10
        )
        if not sol.success:
            print(f"⚠️ two-country case2 failed at τ_l={tau_l_param:.3f}: {sol.message}")
            continue
        sol.x = _polish(lambda v: np.asarray(_err_two_country_case2(
            v, out=0, lam_star=lambda_star, tau_l_param=tau_l_param, τ10h_pref=τ10h_pref)), sol.x)

        p0, p1 = sol.x[:4], sol.x[4:8]
        w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
        w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])
        τ0 = np.array([
            _taul_param(w0[0], lambda_star, tau_l_param),
            _taul_param(w0[1], lambda_star, tau_l_param),
            _taul_param(w0[2], lambda_star, tau_l_param),
            τ10h_pref
        ])
        tau_l1 = sol.x[8]
        τ1 = np.array([
            _taul_param(w1[0], lambda_star, tau_l1),
            _taul_param(w1[1], lambda_star, tau_l1),
            _taul_param(w1[2], lambda_star, tau_l1),
            _taul_param(w1[3], lambda_star, tau_l1)
        ])
        row = _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lambda_star, sol.fun,
                                     tau_l_param1=tau_l1)
        row['tau10h'] = τ0[3]; row['tau01h'] = τ1[3]
        rows.append(row)

        x0 = sol.x  # warm start

    return pd.DataFrame(rows).sort_values('tauL0').reset_index(drop=True)


def sweep_taul_large_case3(lambda_star, tau_grid, τ_pref, x0=None):
    if x0 is None:
        pstar = pi_ge_eq[:4]; x0 = np.concatenate([pstar, pstar])
    bounds = tuple(bounds_lambda_8)

    rows = []
    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda v: _err_two_country_case3(v, out=1, lam_star=lambda_star,
                                             tau_l_param=tau_l_param, τ10h_pref=τ_pref, τ01h_pref=τ_pref),
            x0, method='SLSQP', bounds=bounds, tol=1e-10
        )
        if not sol.success:
            print(f"⚠️ two-country case3 failed at τ_l={tau_l_param:.3f}: {sol.message}")
            continue
        sol.x = _polish(lambda v: np.asarray(_err_two_country_case3(
            v, out=0, lam_star=lambda_star, tau_l_param=tau_l_param, τ10h_pref=τ_pref,
            τ01h_pref=τ_pref)), sol.x)

        p0, p1 = sol.x[:4], sol.x[4:8]
        w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
        w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])
        τ0 = np.array([
            _taul_param(w0[0], lambda_star, tau_l_param),
            _taul_param(w0[1], lambda_star, tau_l_param),
            _taul_param(w0[2], lambda_star, tau_l_param),
            τ_pref
        ])
        τ1 = np.array([
            _taul_param(w1[0], lambda_star, tau_l_param),
            _taul_param(w1[1], lambda_star, tau_l_param),
            _taul_param(w1[2], lambda_star, tau_l_param),
            τ_pref
        ])
        row = _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lambda_star, sol.fun)
        row['tau10h'] = τ0[3]; row['tau01h'] = τ1[3]
        rows.append(row)

        x0 = sol.x  # warm start

    return pd.DataFrame(rows).sort_values('tauL0').reset_index(drop=True)


