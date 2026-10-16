# analysis/calibrate_shape.py
"""
Calibration of the shape parameters σ and θ (app:shape-parameters).

For a given specification (ε, γ, ρ, φ, set through MIG_EPS / MIG_GAMMA / MIG_RHO /
MIG_PHI; without MIG_PHI, φ follows the rule in config.py), find σ such that the
model's pooled native coefficient equals 1.323, with θ solved for each σ so that the arc elasticity of high-skill immigrants over the policy
cut equals 2.392 (analysis/shape_moments.py). σ is set at import, so each σ is
evaluated in its own process; the outer loop is a secant iteration on σ.

The results are written to model/shape_params.json, keyed by specification, rounded
to three decimals (the values used in the paper and in output filenames), with the
unrounded solution and the matched moments at the rounded values.

Usage, from the package root:
    python -m analysis.calibrate_shape                       # baseline (ε = 7, γ = 0)
    MIG_EPS=20 python -m analysis.calibrate_shape            # ε = 20 variant
    MIG_GAMMA=0.05 python -m analysis.calibrate_shape        # land variant
"""
import json
import os
import subprocess
import sys

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARAMS = os.path.join(PKG, "model", "shape_params.json")
ALPHA_TARGET = 1.323


def _moments(sigma: float, theta: float = None) -> dict:
    """Run analysis.shape_moments in a fresh process at σ (θ solved there unless given)."""
    env = os.environ.copy()
    env["MIG_SIGMA"] = repr(float(sigma))
    env["MIG_THETA"] = "15"          # θ ≥ σ at import; θ is passed explicitly
    code = ("import json, analysis.shape_moments as m; "
            + ("t = m.theta_star(); " if theta is None else f"t = {theta!r}; ")
            + "print(json.dumps({'sigma': m.σ, 'theta': t, 'arc_eF': m.arc_eF(t), "
              "'marginal_eF': m.marginal_eF(t), 'alpha': m.pooled_alpha(t)}))")
    r = subprocess.run([sys.executable, "-c", code], cwd=PKG, env=env,
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def calibrate(s0: float = 1.25, s1: float = 1.30, tol: float = 1e-5, maxit: int = 20) -> dict:
    m0, m1 = _moments(s0), _moments(s1)
    f0, f1 = m0["alpha"] - ALPHA_TARGET, m1["alpha"] - ALPHA_TARGET
    for _ in range(maxit):
        s2 = s1 - f1 * (s1 - s0) / (f1 - f0)
        m2 = _moments(s2)
        f2 = m2["alpha"] - ALPHA_TARGET
        print(f"  sigma = {s2:.6f}  theta = {m2['theta']:.6f}  alpha = {m2['alpha']:.6f}", flush=True)
        if abs(s2 - s1) < tol and abs(f2) < 1e-6:
            return m2
        s0, f0, s1, f1, m1 = s1, f1, s2, f2, m2
    raise RuntimeError("shape calibration did not converge")


def main():
    # the specification key as config builds it (φ is recalibrated when ε differs
    # from the baseline), read in a process where σ, θ are set explicitly
    env = os.environ.copy()
    env.update({"MIG_SIGMA": "1.3", "MIG_THETA": "15"})
    key = subprocess.run([sys.executable, "-c", "import config; print(config.SHAPE_KEY)"],
                         cwd=PKG, env=env, capture_output=True, text=True,
                         check=True).stdout.strip().splitlines()[-1]
    print(f"== shape calibration: {key} ==", flush=True)
    sol = calibrate()
    sigma_r, theta_r = round(sol["sigma"], 3), round(sol["theta"], 3)
    at_r = _moments(sigma_r, theta_r)
    entry = {"sigma": sigma_r, "theta": theta_r,
             "sigma_exact": sol["sigma"], "theta_exact": sol["theta"],
             "alpha_target": ALPHA_TARGET, "e_F_target": 2.392,
             "alpha_at_rounded": at_r["alpha"], "arc_eF_at_rounded": at_r["arc_eF"],
             "marginal_eF_at_rounded": at_r["marginal_eF"]}
    params = {}
    if os.path.exists(PARAMS):
        with open(PARAMS, encoding="utf-8") as f:
            params = json.load(f)
    params[key] = entry
    with open(PARAMS, "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2, sort_keys=True)
    print(json.dumps({key: entry}, indent=2))


if __name__ == "__main__":
    main()
