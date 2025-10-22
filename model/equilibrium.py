import numpy as np
import pandas as pd
from scipy import optimize
from .equations import wage, taul, Pi, ri        # + ri
from config import G0, λ0_vals, bounds_lambda, G1, tau_k, ell_l, ell_h
import inspect

# --- Symmetric GE equilibrium solver (baseline) ---
_x0  = np.array([0.4, 0.2, 0.24, 0.1, 0.1])
_bnds = ((1e-5, .75), (1e-5, .75), (1e-5, .4), (1e-5, .4), (1e-5, 10))

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
    #  (ii) capital-tax on origin-0 land income (origin-based, regardless of residence)
    r0 = ri(0, p0)
    origin0_low, origin0_high = 0.64, 0.36   # masses by skill from your constraints
    rev_cap0 = tau_k * r0 * (ell_l * origin0_low + ell_h * origin0_high)

    err = [
        rev_labor0 + rev_cap0 - G0,                               # govt budget (0)
        p0[0] - Pi(0,0,'l', p0, p0, τ0, τ0)*origin0_low,          # natives low in 0
        p0[2] - Pi(0,0,'h', p0, p0, τ0, τ0)*origin0_high,         # natives high in 0
        p0[0] + p0[1] - origin0_low,                              # residents low in 0
        p0[2] + p0[3] - origin0_high                              # residents high in 0
    ]
    return np.sum(np.square(err)) if out == 1 else err


# run this once at import to get pi_ge_eq
_sym = optimize.minimize(err_ge_eq, _x0, method='SLSQP', bounds=_bnds, tol=1e-12)
if not _sym.success:
    raise RuntimeError(f"Symmetric GE solver failed: {_sym.message}")
pi_ge_eq = _sym.x

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

    # only the population-share residuals
    res = [
        p0[0] - Pi(0,0,'l', p0, p0, τ0, τ0)*0.64,
        p0[2] - Pi(0,0,'h', p0, p0, τ0, τ0)*0.36,
        p0[0] + p0[1] - 0.64,
        p0[2] + p0[3] - 0.36
    ]
    return np.sum(np.square(res)) if out==1 else res

# --- Precompute foreign equilibrium values ---
p_foreign = pi_ge_eq[:4]
λ_star = pi_ge_eq[4]
tau_star = 1 - λ_star
w10h_star = wage(1, 0, 'h', pi_ge_eq[:4])   # (1,0) = immigrant, high skill
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

    home_w = {(i, a): wage(i, 0, a, p0) for i, a in [(0, 'l'), (1, 'l'), (0, 'h'), (1, 'h')]}
    τ0 = np.array([taul(home_w[(i, a)], λ0) for i, a in [(0, 'l'), (1, 'l'), (0, 'h'), (1, 'h')]])

    res = [
        # native stay-share conditions (Pi now includes land for staying)
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * 0.64,
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * 0.36,

        # immigrant no-arbitrage on wages
        (1 - τ0[1]) * home_w[(1, 'l')] - (1 - τ1[1]) * foreign_w[(1, 'l')],
        (1 - τ0[3]) * home_w[(1, 'h')] - (1 - τ1[3]) * foreign_w[(1, 'h')],
    ]
    return np.sum(np.square(res)) if out == 1 else res


# --- Prefixed high-skill tax residual ---
def err_ge_pref_tax(params, out=1, custom_λ0=None, custom_τ10h=None):
    p0 = params

    home_w = {(i, a): wage(i, 0, a, p0) for i, a in [(0, 'l'), (1, 'l'), (0, 'h'), (1, 'h')]}
    τ0 = np.array([
        taul(home_w[(0, 'l')], custom_λ0),
        taul(home_w[(1, 'l')], custom_λ0),
        taul(home_w[(0, 'h')], custom_λ0),
        custom_τ10h
    ])

    res = [
        # native stay-share conditions (Pi includes land for staying)
        p0[0] - Pi(0, 0, 'l', p0, p_foreign, τ0, τ1) * 0.64,
        p0[2] - Pi(0, 0, 'h', p0, p_foreign, τ0, τ1) * 0.36,

        # immigrant no-arbitrage on wages only (simple for now)
        (1 - τ0[1]) * home_w[(1, 'l')] - (1 - τ1[1]) * foreign_w[(1, 'l')],
        (1 - τ0[3]) * home_w[(1, 'h')] - (1 - τ1[3]) * foreign_w[(1, 'h')],
    ]
    return np.sum(np.square(res)) if out == 1 else res


# --- Sweep helper ---
def sweep_equilibrium(err_func, λ_vals, bounds, x0=None):
    """
    Symmetric (two identical countries) λ-sweep.
    For each λ in λ_vals:
      - solve err_func(..., custom_λ0=λ)
      - compute wages, taxes τ(y)=taul(y, λ)
      - inc codes (00l,10l,00h,10h) = labor net + origin-based land net of tau_k
      - TR = labor tax + origin-based land tax

    Returns: DataFrame with columns
      λ0, p00l,p10l,p00h,p10h, w00l,w10l,w00h,w10h,
      tau00l,tau10l,tau00h,tau10h,
      inc00l,inc10l,inc00h,inc10h, TR, fun
    """
    if x0 is None:
        x0 = np.array([0.52, 0.10, 0.30, 0.10])

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

    origin0_low, origin0_high = 0.64, 0.36  # fixed origin masses by skill

    for λ in λ_vals:
        sol = optimize.minimize(
            lambda v: err_func(v, out=1, custom_λ0=λ),
            x0=x0, method="SLSQP", bounds=bounds,
            tol=1e-12, options={"ftol": 1e-12, "maxiter": 1000}
        )
        if not sol.success:
            # skip infeasible/failed points to keep the sweep clean
            continue

        p0 = sol.x[:4]

        # wages (country 0 perspective)
        w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
        w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)

        # labor taxes at this λ (override τ10h if err_func pins it)
        tau00l = taul(w00l, λ)
        tau10l = taul(w10l, λ)
        tau00h = taul(w00h, λ)
        tau10h = taul(w10h, λ) if (tau10h_fixed is None) else float(tau10h_fixed)

        # origin-based land returns (symmetric baseline ⇒ r1 == r0)
        r0 = ri(0, p0)
        r1 = r0

        # net land income per code (origin 0 for 00·, origin 1 for 10·)
        land00l = r0 * ell_l * (1.0 - tau_k)
        land10l = r1 * ell_l * (1.0 - tau_k)
        land00h = r0 * ell_h * (1.0 - tau_k)
        land10h = r1 * ell_h * (1.0 - tau_k)

        # consumption (labor net + origin-based land net)
        inc00l = w00l * (1 - tau00l) + land00l
        inc10l = w10l * (1 - tau10l) + land10l
        inc00h = w00h * (1 - tau00h) + land00h
        inc10h = w10h * (1 - tau10h) + land10h

        # government revenue in country 0: labor + land (origin-based)
        rev_labor0 = (
            p0[0]*tau00l*w00l + p0[1]*tau10l*w10l +
            p0[2]*tau00h*w00h + p0[3]*tau10h*w10h
        )
        rev_cap0 = tau_k * r0 * (ell_l*origin0_low + ell_h*origin0_high)
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

        # warm start next iteration
        x0 = p0.copy()

    return pd.DataFrame(out)




###############################################################################
################### LARGE COUNTRY EQUILIBRIA ##################################
###############################################################################
# --- Two-country residual function ---
def err_2country(
    params,
    out=1,
    λ0=None,
    λ1=None,
    tau_high_imm_override=None,  # GLOBAL: τ10h (C0) and τ01h (C1)
    tau10h_override=None,        # UNILATERAL: τ10h (C0) only; ignored if global provided
    mode='GE'                    # 'GE' adds C1 GBC and takes λ1 from params[-1]; 'PE' takes λ1 from arg
):
    if λ0 is None:
        raise ValueError("Pass λ0 via argument.")

    # unpack
    p0 = params[0:4]  # [p00l, p10l, p00h, p10h]
    p1 = params[4:8]  # [p11l, p01l, p11h, p01h]

    λ0_val = λ0
    λ1_val = (params[-1] if mode == 'GE' else λ1)

    # wages
    w00l = wage(0,0,'l',p0);  w10l = wage(1,0,'l',p0)
    w00h = wage(0,0,'h',p0);  w10h = wage(1,0,'h',p0)
    w11l = wage(1,1,'l',p1);  w01l = wage(0,1,'l',p1)
    w11h = wage(1,1,'h',p1);  w01h = wage(0,1,'h',p1)

    # taxes: C0 [τ00l, τ10l, τ00h, τ10h]
    τ00l = taul(w00l, λ0_val)
    τ10l = taul(w10l, λ0_val)
    τ00h = taul(w00h, λ0_val)
    if tau_high_imm_override is not None:
        τ10h = float(tau_high_imm_override)   # GLOBAL
    elif tau10h_override is not None:
        τ10h = float(tau10h_override)         # UNILATERAL
    else:
        τ10h = taul(w10h, λ0_val)
    τ0 = np.array([τ00l, τ10l, τ00h, τ10h])

    # taxes: C1 [τ11l, τ01l, τ11h, τ01h]
    τ11l = taul(w11l, λ1_val)
    τ01l = taul(w01l, λ1_val)
    τ11h = taul(w11h, λ1_val)
    if tau_high_imm_override is not None:
        τ01h = float(tau_high_imm_override)   # GLOBAL
    else:
        τ01h = taul(w01h, λ1_val)
    τ1 = np.array([τ11l, τ01l, τ11h, τ01h])

    # residuals
    err = [
        # C0
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1) * 0.64,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1) * 0.36,
        p0[0] + p1[1] - 0.64,
        p0[2] + p1[3] - 0.36,
        # C1
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0) * 0.64,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0) * 0.36,
        p1[0] + p0[1] - 0.64,
        p1[2] + p0[3] - 0.36,
    ]

    # C1 GBC (only GE)
    if mode == 'GE':
        TR1_lab = (
            p1[0]*τ11l*w11l + p1[1]*τ01l*w01l +
            p1[2]*τ11h*w11h + p1[3]*τ01h*w01h
        )
        origin1_low, origin1_high = 0.64, 0.36
        r1 = ri(1, p1)
        TR1_cap = tau_k * r1 * (ell_l * origin1_low + ell_h * origin1_high)
        gbc1 = TR1_lab + TR1_cap - G1
        err.insert(0, gbc1)

    return np.sum(np.square(err)) if out == 1 else err

def sweep_equilibrium_2country(
    err_func,
    λ0_vals,
    λ1_fixed,
    bounds,
    x0=None,
    free_lambda1=False,
    **err_kwargs
):
    if free_lambda1:
        bounds = bounds + ((1e-5, 1.5),)

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

    origin_low, origin_high = 0.64, 0.36

    # resolve overrides from caller
    tau_hi_imm_global = err_kwargs.get("tau_high_imm_override", None)                         # explicit GLOBAL
    tau10h_unilateral = err_kwargs.get("tau10h_override", err_kwargs.get("custom_τ10h", None))# UNILATERAL

    # if user called the GE-global residual with custom_τ10h, treat it as GLOBAL for reporting
    if (tau_hi_imm_global is None) and (err_func.__name__ == "ge_pref_residual_global"):
        tau_hi_imm_global = err_kwargs.get("custom_τ10h", None)

    # filter kwargs to err_func signature (prevents PE funcs from choking)
    possible_policy_kwargs = {
        "tau_high_imm_override": tau_hi_imm_global,
        "tau10h_override": tau10h_unilateral,
        "custom_τ10h": err_kwargs.get("custom_τ10h", None),
    }
    sig = set(inspect.signature(err_func).parameters.keys())
    filtered_policy_kwargs = {k: v for k, v in possible_policy_kwargs.items() if (k in sig and v is not None)}

    for λ0 in λ0_vals:
        sol = optimize.minimize(
            lambda params: err_func(
                params,
                1,
                λ0,
                (None if free_lambda1 else λ1_fixed),
                **filtered_policy_kwargs
            ),
            x0, method="SLSQP", bounds=bounds, tol=1e-10
        )

        if sol.success:
            p = sol.x
            p0, p1 = p[:4], p[4:8]
            λ1_val = sol.x[8] if free_lambda1 else λ1_fixed

            # wages
            w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
            w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])

            # taxes (reporter mirrors policy semantics)
            # C0: [τ00l, τ10l, τ00h, τ10h]
            τ0 = np.array([
                taul(w0[0], λ0),
                taul(w0[1], λ0),
                taul(w0[2], λ0),
                ( float(tau_hi_imm_global) if tau_hi_imm_global is not None
                  else ( float(tau10h_unilateral) if tau10h_unilateral is not None
                         else taul(w0[3], λ0) ))
            ])

            # C1: [τ11l, τ01l, τ11h, τ01h]
            τ1_list = [taul(w1[i], λ1_val) for i in range(4)]
            if tau_hi_imm_global is not None:
                τ1_list[3] = float(tau_hi_imm_global)  # GLOBAL forces τ01h
            τ1 = np.array(τ1_list)

            # asserts
            if tau_hi_imm_global is not None:
                assert np.isclose(τ0[3], tau_hi_imm_global), "τ10h (C0) not enforced (GLOBAL)"
                assert np.isclose(τ1[3], tau_hi_imm_global), "τ01h (C1) not enforced (GLOBAL)"
            if tau10h_unilateral is not None and tau_hi_imm_global is None:
                assert np.isclose(τ0[3], tau10h_unilateral), "τ10h (C0) not enforced (UNILATERAL)"

            # revenues (labor + origin-based capital; if r=0 → just labor)
            r0 = ri(0, p0); r1 = ri(1, p1)
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

            TR0_lab = np.dot(p0, τ0 * w0)
            TR1_lab = np.dot(p1, τ1 * w1)
            TR0_cap = tau_k * r0 * (ell_l * origin_low + ell_h * origin_high)
            TR1_cap = tau_k * r1 * (ell_l * origin_low + ell_h * origin_high)
            TR0 = TR0_lab + TR0_cap
            TR1 = TR1_lab + TR1_cap

            row = [λ0, λ1_val, *p0, *p1, *w0, *w1, *τ0, *τ1, *inc0, *inc1, TR0, TR1, sol.fun]
        else:
            λ1_val = sol.x[8] if (free_lambda1 and hasattr(sol, "x") and len(sol.x) >= 9) else (np.nan if free_lambda1 else λ1_fixed)
            print(f"⚠️ failed at λ0={λ0}: {sol.message}")
            row = [λ0, λ1_val] + [np.nan] * (len(cols) - 2)

        for c, v in zip(cols, row):
            results[c].append(v)

    return pd.DataFrame(results)



# --- (PE) Residual with preferential τ10h in Country 0 -----------------------
def err_pe_large_tau10h(params, out=1, custom_λ0=None, custom_λ1=None, custom_τ10h=None):
    """
    Two-country PE residual:
      • Country 0 sets τ10h = custom_τ10h (preferential on high-skill immigrants).
      • Country 1 uses its regular schedule.
      • No government budget constraints here (PE), so no land-income terms.
    """
    # unpack populations
    p0 = params[0:4]   # [p00l, p10l, p00h, p10h]
    p1 = params[4:8]   # [p11l, p01l, p11h, p01h]

    λ0 = custom_λ0
    λ1 = custom_λ1
    τ10h_pref = custom_τ10h

    # wages
    w0 = np.array([
        wage(0,0,'l',p0), wage(1,0,'l',p0),
        wage(0,0,'h',p0), wage(1,0,'h',p0)
    ])
    w1 = np.array([
        wage(1,1,'l',p1), wage(0,1,'l',p1),
        wage(1,1,'h',p1), wage(0,1,'h',p1)
    ])

    # taxes (override τ10h in country 0)
    τ0 = np.array([
        taul(w0[0], λ0),
        taul(w0[1], λ0),
        taul(w0[2], λ0),
        τ10h_pref
    ])
    τ1 = np.array([
        taul(w1[0], λ1),
        taul(w1[1], λ1),
        taul(w1[2], λ1),
        taul(w1[3], λ1)
    ])

    # residuals: native stay shares and population sums (both countries)
    err = [
        # Country 0
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1) * 0.64,
        p0[0] + p1[1] - 0.64,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1) * 0.36,
        p0[2] + p1[3] - 0.36,

        # Country 1
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0) * 0.64,
        p1[0] + p0[1] - 0.64,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0) * 0.36,
        p1[2] + p0[3] - 0.36,
    ]

    return np.sum(np.square(err)) if out == 1 else err


# --- equilibrium glue ---


def ge_pref_residual_global(params, out, custom_λ0, custom_λ1=None, custom_τ10h=None):
    """
    GE global preferential tax: force τ10h (Country 0) and τ01h (Country 1) = custom_τ10h.
    λ1 floats inside params[-1] under mode='GE'.
    """
    return err_2country(
        params=params,
        out=out,
        λ0=custom_λ0,
        λ1=None,                           # GE: λ1 floats in params[-1]
        tau_high_imm_override=custom_τ10h, # GLOBAL override → τ10h & τ01h
        mode='GE'
    )







###############################################################################
############### CHANGING PROGRESSIVITY PARAMETER ##############################
###############################################################################

# --- τ_l helpers (progressivity sweeps with λ fixed) -------------------------
def _taul_param(y, lam, tau_l_param):
    return 1.0 - lam * y**(-tau_l_param)


def err_ge_soe_pref_tau(params, out=1, lam=None, tau_l_param=None, τ10h_pref=None):
    """
    SOE residual like above, but override τ10h with a preferential tax.
    Still no GBC; land affects only the stay-vs-move via Pi.
    """
    p0 = params[:4]
    w00l = wage(0,0,'l',p0); w10l = wage(1,0,'l',p0)
    w00h = wage(0,0,'h',p0); w10h = wage(1,0,'h',p0)

    τ00l = _taul_param(w00l, lam, tau_l_param)
    τ10l = _taul_param(w10l, lam, tau_l_param)
    τ00h = _taul_param(w00h, lam, tau_l_param)
    τ10h = τ10h_pref
    τ0 = np.array([τ00l, τ10l, τ00h, τ10h])

    res = [
        p0[0] - Pi(0,0,'l', p0, p_foreign, τ0, τ1)*0.64,
        p0[2] - Pi(0,0,'h', p0, p_foreign, τ0, τ1)*0.36,
        (1-τ10l)*w10l - (1-τ1[1]) * foreign_w[(1,'l')],
        (1-τ10h)*w10h - (1-τ1[3]) * foreign_w[(1,'h')],
    ]
    return np.sum(np.square(res)) if out==1 else res


def sweep_taul_soe_pref(lambda_star, tau_grid, τ10h_pref, x0=None, bnds=None):
    if x0 is None:
        x0 = np.array([0.52, 0.12, 0.30, 0.06])
    if bnds is None:
        # (p00l, p10l, p00h, p10h)
        bnds = ((1e-5, .64), (1e-5, .64), (1e-5, .36), (1e-5, .36))

    rows = []
    origin_low, origin_high = 0.64, 0.36  # origin masses by skill

    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda x: err_ge_soe_pref_tau(
                x, out=1, lam=lambda_star, tau_l_param=tau_l_param, τ10h_pref=τ10h_pref
            ),
            x0, method='SLSQP', bounds=bnds, tol=1e-10
        )
        if not sol.success:
            continue

        p0 = sol.x[:4]
        w  = [wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)]
        τ  = [
            _taul_param(w[0], lambda_star, tau_l_param),
            _taul_param(w[1], lambda_star, tau_l_param),
            _taul_param(w[2], lambda_star, tau_l_param),
            τ10h_pref
        ]

        # --- Government revenue: labor + origin-based land (already correct)
        TR_lab = p0[0]*τ[0]*w[0] + p0[1]*τ[1]*w[1] + p0[2]*τ[2]*w[2] + p0[3]*τ[3]*w[3]
        r0 = ri(0, p0)                      # home land return (origin 0)
        TR_cap = tau_k * r0 * (ell_l * origin_low + ell_h * origin_high)
        TR = TR_lab + TR_cap

        # --- Add origin-based net land income to consumption (SOE: foreign side fixed)
        r1_foreign = ri(1, p_foreign)       # foreign land return (origin 1, fixed at baseline)

        land00l = r0 * ell_l * (1 - tau_k)  # origin 0, low
        land00h = r0 * ell_h * (1 - tau_k)  # origin 0, high
        land10l = r1_foreign * ell_l * (1 - tau_k)  # origin 1, low
        land10h = r1_foreign * ell_h * (1 - tau_k)  # origin 1, high

        inc00l = w[0]*(1-τ[0]) + land00l
        inc10l = w[1]*(1-τ[1]) + land10l
        inc00h = w[2]*(1-τ[2]) + land00h
        inc10h = w[3]*(1-τ[3]) + land10h

        rows.append(dict(
            tau_l=tau_l_param, TR=TR,
            p00l=p0[0], p10l=p0[1], p00h=p0[2], p10h=p0[3],
            w00l=w[0], w10l=w[1], w00h=w[2], w10h=w[3],
            tau00l=τ[0], tau10l=τ[1], tau00h=τ[2], tau10h=τ[3],
            inc00l=inc00l, inc10l=inc10l, inc00h=inc00h, inc10h=inc10h,
            fun=sol.fun
        ))
    return pd.DataFrame(rows).sort_values('tau_l').reset_index(drop=True)




def _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lam_star, fun_val):
    """
    Build a row for the two-country sweeps.

    - Incomes (inc··) = net labor + origin-based net land income
      * Country 0 origin land return: r0 = ri(0, p0)
      * Country 1 origin land return: r1 = ri(1, p1)
      * Net land per capita for skill a ∈ {l,h} is r_origin * ell_a * (1 - tau_k),
        regardless of residence.
    - Tax revenue (TR0/TR1) = labor tax + origin-based land tax (already correct).
    """
    # land returns by origin
    r0 = ri(0, p0)
    r1 = ri(1, p1)

    # origin masses by skill (used only for revenue)
    origin_low, origin_high = 0.64, 0.36

    # --- net land components (origin-based), per capita
    # Country 0 display codes (origin 0 for 00·, origin 1 for 10·):
    land00l = r0 * ell_l * (1 - tau_k)
    land00h = r0 * ell_h * (1 - tau_k)
    land10l = r1 * ell_l * (1 - tau_k)
    land10h = r1 * ell_h * (1 - tau_k)

    # Country 1 display codes (origin 1 for 11·, origin 0 for 01·):
    land11l = r1 * ell_l * (1 - tau_k)
    land11h = r1 * ell_h * (1 - tau_k)
    land01l = r0 * ell_l * (1 - tau_k)
    land01h = r0 * ell_h * (1 - tau_k)

    # --- net incomes = net labor + net land (origin-based)
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

    # --- tax revenues: labor + origin-based land (unchanged logic)
    TR0_lab = float(np.dot(p0, τ0 * w0))
    TR1_lab = float(np.dot(p1, τ1 * w1))
    TR0_cap = tau_k * r0 * (ell_l * origin_low + ell_h * origin_high)
    TR1_cap = tau_k * r1 * (ell_l * origin_low + ell_h * origin_high)
    TR0 = TR0_lab + TR0_cap
    TR1 = TR1_lab + TR1_cap

    return {
        'tauL0': tau_l_param, 'tauL1': tau_l_param,
        # populations
        'p00l': p0[0], 'p10l': p0[1], 'p00h': p0[2], 'p10h': p0[3],
        'p11l': p1[0], 'p01l': p1[1], 'p11h': p1[2], 'p01h': p1[3],
        # wages
        'w00l': w0[0], 'w10l': w0[1], 'w00h': w0[2], 'w10h': w0[3],
        'w11l': w1[0], 'w01l': w1[1], 'w11h': w1[2], 'w01h': w1[3],
        # taxes
        'tau00l': τ0[0], 'tau10l': τ0[1], 'tau00h': τ0[2], 'tau10h': τ0[3],
        'tau11l': τ1[0], 'tau01l': τ1[1], 'tau11h': τ1[2], 'tau01h': τ1[3],
        # after-tax incomes (now include land)
        'inc00l': inc0[0], 'inc10l': inc0[1], 'inc00h': inc0[2], 'inc10h': inc0[3],
        'inc11l': inc1[0], 'inc01l': inc1[1], 'inc11h': inc1[2], 'inc01h': inc1[3],
        # revenues (also duplicated as G0/G1)
        'TR0': TR0, 'TR1': TR1, 'G0': TR0, 'G1': TR1,
        'λ0': lam_star, 'λ1': lam_star, 'fun': fun_val
    }


def _err_two_country_case2(params, out, lam_star, tau_l_param, τ10h_pref):
    """Case 2: pref only on τ10h (Country 0). τ_l same in both countries."""
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
        _taul_param(w1[3], lam_star, tau_l_param)
    ])
    err = [
        # Country 0 stay/sum
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1)*0.64,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1)*0.36,
        p0[0] + p1[1] - 0.64,
        p0[2] + p1[3] - 0.36,
        # Country 1 stay/sum
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0)*0.64,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0)*0.36,
        p1[0] + p0[1] - 0.64,
        p1[2] + p0[3] - 0.36,
    ]
    return np.sum(np.square(err)) if out==1 else err


def _err_two_country_case3(params, out, lam_star, tau_l_param, τ10h_pref, τ01h_pref):
    """Case 3: pref on τ10h (Country 0) and τ01h (Country 1). τ_l same in both."""
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
        # Country 0 stay/sum
        p0[0] - Pi(0,0,'l', p0, p1, τ0, τ1)*0.64,
        p0[2] - Pi(0,0,'h', p0, p1, τ0, τ1)*0.36,
        p0[0] + p1[1] - 0.64,
        p0[2] + p1[3] - 0.36,
        # Country 1 stay/sum
        p1[0] - Pi(1,1,'l', p1, p0, τ1, τ0)*0.64,
        p1[2] - Pi(1,1,'h', p1, p0, τ1, τ0)*0.36,
        p1[0] + p0[1] - 0.64,
        p1[2] + p0[3] - 0.36,
    ]
    return np.sum(np.square(err)) if out==1 else err


def sweep_taul_large_case2(lambda_star, tau_grid, τ10h_pref, x0=None, bnds=None):
    """Country 0 adopts preferential τ10h; Country 1 does not. τ_l varies (same in both)."""
    if x0 is None:
        pstar = pi_ge_eq[:4]
        x0 = np.concatenate([pstar, pstar])
    if bnds is None:
        b4 = ((1e-5, .64), (1e-5, .36), (1e-5, .64), (1e-5, .36))
        bnds = b4 + b4

    rows = []
    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda v: _err_two_country_case2(v, out=1, lam_star=lambda_star, tau_l_param=tau_l_param, τ10h_pref=τ10h_pref),
            x0, method='SLSQP', bounds=bnds, tol=1e-10
        )
        if not sol.success:
            print(f"⚠️ two-country case2 failed at τ_l={tau_l_param:.3f}: {sol.message}")
            continue

        p0, p1 = sol.x[:4], sol.x[4:8]
        w0 = np.array([wage(0,0,'l',p0), wage(1,0,'l',p0), wage(0,0,'h',p0), wage(1,0,'h',p0)])
        w1 = np.array([wage(1,1,'l',p1), wage(0,1,'l',p1), wage(1,1,'h',p1), wage(0,1,'h',p1)])
        τ0 = np.array([
            _taul_param(w0[0], lambda_star, tau_l_param),
            _taul_param(w0[1], lambda_star, tau_l_param),
            _taul_param(w0[2], lambda_star, tau_l_param),
            τ10h_pref
        ])
        τ1 = np.array([
            _taul_param(w1[0], lambda_star, tau_l_param),
            _taul_param(w1[1], lambda_star, tau_l_param),
            _taul_param(w1[2], lambda_star, tau_l_param),
            _taul_param(w1[3], lambda_star, tau_l_param)
        ])
        # NOTE: TR0/TR1 include land-tax revenue inside _build_row_two_country
        row = _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lambda_star, sol.fun)
        row['tau10h'] = τ0[3]; row['tau01h'] = τ1[3]
        rows.append(row)

    return pd.DataFrame(rows).sort_values('tauL0').reset_index(drop=True)


def sweep_taul_large_case3(lambda_star, tau_grid, τ_pref, x0=None, bnds=None):
    """Both countries adopt preferential tax on immigrant high-skill (τ10h = τ01h = τ_pref)."""
    if x0 is None:
        pstar = pi_ge_eq[:4]
        x0 = np.concatenate([pstar, pstar])
    if bnds is None:
        b4 = ((1e-5, .64), (1e-5, .36), (1e-5, .64), (1e-5, .36))
        bnds = b4 + b4

    rows = []
    for tau_l_param in tau_grid:
        sol = optimize.minimize(
            lambda v: _err_two_country_case3(v, out=1, lam_star=lambda_star, tau_l_param=tau_l_param,
                                             τ10h_pref=τ_pref, τ01h_pref=τ_pref),
            x0, method='SLSQP', bounds=bnds, tol=1e-10
        )
        if not sol.success:
            print(f"⚠️ two-country case3 failed at τ_l={tau_l_param:.3f}: {sol.message}")
            continue

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
        # NOTE: TR0/TR1 include land-tax revenue inside _build_row_two_country
        row = _build_row_two_country(p0, p1, w0, w1, τ0, τ1, tau_l_param, lambda_star, sol.fun)
        row['tau10h'] = τ0[3]; row['tau01h'] = τ1[3]
        rows.append(row)

    return pd.DataFrame(rows).sort_values('tauL0').reset_index(drop=True)
