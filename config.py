### config.py
import dataclasses
import os
import numpy as np

from model.calibration import Targets, calibrate, immigrant_wage_ratio, phi_for_wage_ratio


def _env_param(name, default):
    # variant runs (σ=10, ε=20, γ=0.05) set these in the environment; defaults
    # reproduce the baseline calibration bit-for-bit
    v = os.environ.get(name)
    if v is None:
        return default
    f = float(v)
    return int(f) if f.is_integer() else f


def _fmt(x):
    # filename-safe parameter value: 1 → "1", 0.5 → "0p5", inf → "inf"
    if np.isinf(x):
        return "inf"
    return str(int(x)) if float(x).is_integer() else str(x).replace(".", "p")


# --- Calibration targets (app:data, app:calibration)
# EU 2024, ages 25–64: population by country of birth (Eurostat migr_pop3ctb) and
# tertiary shares by country of birth (edat_lfs_9915); G = 20% of output.
# College premium: Structure of Earnings Survey 2022 (Eurostat earn_ses22_30, mean
# annual gross earnings, and earn_ses22_04, employees), EU27, enterprises with 10
# or more employees, NACE B–S except O, by ISCED 2011 level; the non-college value
# is the employee-weighted mean of levels 0–2 and 3–4 (data/eurostat/README.md).
POP_25_64 = 238_238_721
FOREIGN_25_64 = 42_027_781
SES22_EARNINGS = {"ED0-2": 27_165, "ED3_4": 34_282, "ED5-8": 50_948}        # EUR per year
SES22_EMPLOYEES = {"ED0-2": 21_567_667, "ED3_4": 60_849_551, "ED5-8": 43_914_560}
_noncollege_earnings = (
    (SES22_EMPLOYEES["ED0-2"] * SES22_EARNINGS["ED0-2"]
     + SES22_EMPLOYEES["ED3_4"] * SES22_EARNINGS["ED3_4"])
    / (SES22_EMPLOYEES["ED0-2"] + SES22_EMPLOYEES["ED3_4"])
)
TARGETS_BASE = Targets(
    foreign_share=FOREIGN_25_64 / POP_25_64,
    tertiary_native=0.368,
    tertiary_foreign=0.329,
    skill_premium=SES22_EARNINGS["ED5-8"] / _noncollege_earnings,
    g_share=0.20,
)
# EU-born variant (MIG_EU_BORN=1): immigrants are the residents born in another EU
# country. Immigrant share = EU-born residents 25–64 over all residents 25–64 (migr_pop3ctb
# EU27_2020_FOR, EU27_2020, 1 Jan 2024; the non-EU-born count among the stayers);
# tertiary share of the EU-born 0.348 (edat_lfs_9915, EU27_2020_FOR, 2024, release of
# 2026-09-10); natives' tertiary share unchanged (data/eurostat/eu_born/README.md).
EU_BORN = os.environ.get("MIG_EU_BORN", "0") == "1"
EU_BORN_25_64 = 10_919_080
if EU_BORN:
    TARGETS = dataclasses.replace(TARGETS_BASE, foreign_share=EU_BORN_25_64 / POP_25_64,
                                  tertiary_foreign=0.348)
else:
    TARGETS = TARGETS_BASE

# --- Preferences & technology
γ = _env_param("MIG_GAMMA", 0)
ρ = _env_param("MIG_RHO", 3)
ϵ = _env_param("MIG_EPS", 7)
# φ: relative productivity of immigrants. Baseline 0.72 at ε = 7: immigrants earn
# 0.89 (low skill) and 0.91 (high skill) of the native wage, 0.90 on average. In a
# specification with another ε, φ is recalibrated so that the mean immigrant/native
# wage ratio at the targeted allocation is the same as at the baseline (closed form,
# model/calibration.py). MIG_PHI overrides it.
PHI_BASE, EPS_BASE = 0.72, 7
# (the target is the wage ratio at the baseline targets, also in the EU-born variant)
WAGE_RATIO_TARGET = immigrant_wage_ratio(TARGETS_BASE.allocation()[0], EPS_BASE, PHI_BASE)
if os.environ.get("MIG_PHI") is not None:
    φ = _env_param("MIG_PHI", PHI_BASE)
elif ϵ == EPS_BASE and not EU_BORN:
    φ = PHI_BASE
else:
    φ = round(phi_for_wage_ratio(TARGETS.allocation()[0], ϵ, WAGE_RATIO_TARGET), 3)

# σ: stay/leave margin (home vs. best offer abroad); θ: choice among foreign
# destinations (θ ≥ σ; θ = inf is the perfect-substitutes limit). θ enters only the
# small-country (differentiated-destinations) solver. Both are calibrated so that the
# model reproduces two estimates (app:shape-parameters; analysis/calibrate_shape.py):
# θ, the arc elasticity of high-skill immigrants over the policy cut equals 2.392
# (Kleven et al. 2014, Table II, A4, col. 3, controlling for imperfect take-up);
# σ, the pooled native logit coefficient with wages adjusting equals 1.323 (Kleven,
# Landais, Saez 2013, Table 2, panel A, col. 1). Calibrated values are stored by
# specification (ε, γ, ρ, φ) in model/shape_params.json; a specification without an
# entry uses the baseline entry, except the EU-born variant. MIG_SIGMA / MIG_THETA override them.
_SHAPE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "model", "shape_params.json")
SHAPE_KEY = (f"eps={float(ϵ):g},gamma={float(γ):g},rho={float(ρ):g},phi={float(φ):g}"
             + (",targets=eu_born" if EU_BORN else ""))
_SHAPE_BASE_KEY = "eps=7,gamma=0,rho=3,phi=0.72"


def _shape_default(name):
    import json
    try:
        with open(_SHAPE_FILE, encoding="utf-8") as f:
            params = json.load(f)
    except FileNotFoundError:
        return None
    # the EU-born variant has no fallback: its σ, θ must be calibrated to its own targets
    entry = (params.get(SHAPE_KEY) if EU_BORN
             else params.get(SHAPE_KEY, params.get(_SHAPE_BASE_KEY)))
    return None if entry is None else entry[name]


σ = _env_param("MIG_SIGMA", _shape_default("sigma"))
θ = (float("inf") if os.environ.get("MIG_THETA", "").lower() == "inf"
     else _env_param("MIG_THETA", _shape_default("theta")))
assert σ is not None and θ is not None, \
    f"no calibrated σ, θ in {_SHAPE_FILE}: run analysis.calibrate_shape or set MIG_SIGMA / MIG_THETA"
assert θ >= σ, f"θ = {θ} < σ = {σ}: the nested taste-shock distribution requires θ ≥ σ"

# output-filename tags: θ-independent outputs carry σ only
SIG_TAG = f"sig_{_fmt(σ)}"
DD_TAG = f"sig_{_fmt(σ)}_theta_{_fmt(θ)}"


# --- Factor endowments & TFP
A0, A1 = 1.0, 1.0
L0, L1 = 1.0, 1.0


# --- Tax Progressivity parameter
# Average progressivity of the estimated tax functions of married couples in 17 OECD
# countries, Holter, Krueger and Stepanchuk ("Until the IRS Do Us Part", 2023, p. 13);
# the net income they fit deducts income tax and employee social security contributions
# and adds benefits (their Appendix A.7.1). The calibrated σ, θ in model/shape_params.json
# are for this value.
τl = 0.218

# --- Capital tax
tau_k = 0.18

# --- Preferential tax on high skilled foreigners
pref_tax = 0.114   # half the baseline average rate of high-skill immigrants (0.228)
# Share of government spending proportional to the resident population (rival share):
# individual consumption P31_S13 / government final consumption P3_S13, EU27, Eurostat
# nama_10_gdp, current prices (0.621 in 2024; 0.620-0.630 over 2015-2025)
g_rival_share = 0.62

# --- Lambda grid parameters
λ0_vals = np.linspace(1, 0.75, 250)
# solver bounds on [natives l, immigrants l, natives h, immigrants h]: natives cannot
# exceed their origin mass; immigrant stocks can exceed it (they come from a
# measure-one world), and do in the perfect-substitutes limit at high ϵ
bounds_lambda = ((1e-5, 1.0), (1e-5, 5.0), (1e-5, 1.0), (1e-5, 5.0))
bounds_lambda_8 = bounds_lambda + bounds_lambda 

# --- Tau grid
tau_grid = np.linspace(0.05, 0.40, 40)

# --- Preferential tau grid (laffer curve for changing tau10h)
tau10h_vals = np.linspace(0, 0.5, 25)

# --- Initial conditions
x0  = np.array([0.4, 0.2, 0.24, 0.1, 0.1])
bnds = bounds_lambda + ((1e-5, 10),)

# --- Land endowments by skill (app:land): ℓ^h/ℓ^l equals the skill premium target
if os.environ.get("MIG_LEGACY_CALIB", "0") == "1":
    ell_l = 1.0 / (1.0 + 1.46)
    ell_h = 1.46 / (1.0 + 1.46)
else:
    # level: natives own the whole land stock, Π^l ℓ^l + Π^h ℓ^h = L (= 1), so land
    # rents are fully distributed and the land market clears (app:land)
    _, (_Pi_l, _Pi_h) = TARGETS.allocation()
    ell_l = L0 / (_Pi_l + TARGETS.skill_premium * _Pi_h)
    ell_h = TARGETS.skill_premium * ell_l

# MIG_LEGACY_CALIB=1 restores the previous rounded calibration (σ = 2 table κ,
# S = 1.32, G = 0.35, origin masses 0.64/0.36) for regression checks only.
LEGACY_CALIB = os.environ.get("MIG_LEGACY_CALIB", "0") == "1"

if LEGACY_CALIB:
    _LEGACY_KAPPA = {(2, 7): (2.3, 2.6), (2, 20): (2.05, 2.25), (10, 7): (1.26, 1.30)}
    S = 1.32
    G0 = G1 = 0.35
    origin0_low, origin0_high = 0.64, 0.36
    κl_cal, κh_cal = _LEGACY_KAPPA.get((σ, ϵ), _LEGACY_KAPPA[(σ, 7)])
    CALIBRATION = None
else:
    # every specification (σ, ϵ, γ) is recalibrated to the same data targets
    CALIBRATION = calibrate(TARGETS, sigma=σ, rho=ρ, eps=ϵ, phi=φ, tau_prog=τl,
                            gamma=γ, tau_k=tau_k, ell=(ell_l, ell_h))
    S = CALIBRATION.S
    G0 = G1 = CALIBRATION.G
    origin0_low, origin0_high = CALIBRATION.origin_low, CALIBRATION.origin_high
    κl_cal, κh_cal = CALIBRATION.kappa_l, CALIBRATION.kappa_h









