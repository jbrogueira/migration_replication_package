### config.py
import numpy as np

# --- Preferences & technology
σ = 2
γ = 0.05
ρ = 3
ϵ = 7
#ϵ = 20
#S = 1.85
S = 1.32
φ = 0.9


# --- Factor endowments & TFP
A0, A1 = 1.0, 1.0
L0, L1 = 1.0, 1.0

# --- Land income shares by skill
ell_l = 1.0 / (1.0 + 1.46) #1.46 is the skill premium
ell_h = 1.46    / (1.0 + 1.46)

# --- Public expenditure
G0, G1 = 0.35, 0.35

# --- Tax Progressivity parameter
τl = 0.216


tau_k = 0.18


pref_tax = 0.13

# --- Sweep parameters
#λ0_vals = np.linspace(0.975, 0.80, 50)
λ0_vals = np.linspace(0.975, 0.6, 50)  # wider & denser
#λ0_vals = np.linspace(1.2, 0.5, 50)
#λ0_vals = np.linspace(1.05, 0.75, 50)

    
#λ0_vals = np.linspace(1.05, 0.85, 10) 
bounds_lambda = ((1e-5, .64), (1e-5, 10), (1e-5, .4), (1e-5, 10))
bounds_lambda_8 = bounds_lambda + bounds_lambda 










