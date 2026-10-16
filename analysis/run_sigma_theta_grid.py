# analysis/run_sigma_theta_grid.py
"""
Sensitivity of the main results to σ (stay/leave margin) and θ (choice among foreign
destinations), at ρ = 3 and φ = 0.72 (appendix tables).

Table 1: native CEV in the small country with the preferential tax, for each (σ, θ)
with θ ≥ σ, and the high-skill immigrant share.
Table 2: results by σ: κ, the θ at which the native CEV of each skill turns positive,
the self-financing θ, and the CEVs in the perfect-substitutes limit, the large
country and multilateral adoption.

Every σ is fully recalibrated to the same data targets (model/calibration.py). Each σ
runs run_all.py in a scratch copy of the package (output filenames carry σ only, so
σ = 1 would overwrite the main run); the results are collected under
output/variants/sigma_theta_grid/.

Usage, from the package root:  python -m analysis.run_sigma_theta_grid
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PKG, "output", "variants", "sigma_theta_grid")

# baseline σ, θ: the calibrated values in model/shape_params.json (config.py defaults)
with open(os.path.join(PKG, "model", "shape_params.json"), encoding="utf-8") as _f:
    _BASE = json.load(_f)["eps=7,gamma=0,rho=3,phi=0.72"]
SIGMA_BASE, THETA_BASE = _BASE["sigma"], _BASE["theta"]   # SIGMA_BASE's column carries the immigrant share
SIGMAS = (1.04, SIGMA_BASE, 1.8, 2.5)
# θ grid: the baseline, the bounds of its 95 percent confidence interval (3.07, 4.35;
# app:shape-parameters), and θ = 2 to 4.5, the range implied by the Dutch estimates
THETAS = (2, 2.5, 3, 3.07, THETA_BASE, 4.35, 4.5)
RHO, PHI = 3, 0.72


def _cev_row(tex_path: str) -> list:
    # CEV row of a baseline-vs-policy table written by build_equilibrium_tables
    with open(tex_path, encoding="utf-8") as f:
        line = next(l for l in f if l.startswith("{CEV}"))
    return [float(v) for v in re.findall(r"([+-]\d+\.\d)\\%", line)]


def cell():
    """Run inside a copy after run_all.py: write the σ-specific results."""
    import config
    from config import σ, pref_tax, τl, SIG_TAG
    from model.equilibrium import pi_ge_eq
    from model.equilibrium_diff_dest import equilibrium_row
    from analysis.run_theta_sweep import solve_theta_path

    row_b = equilibrium_row(pi_ge_eq[:4], pi_ge_eq[4])
    df, _ = solve_theta_path([t for t in THETAS if t >= σ], row_b)
    df["sigma"] = σ
    df["p10h_base"] = row_b["p10h"]
    df.to_csv(f"output/tables/csv/sigma_theta_grid_{SIG_TAG}.csv", index=False)

    tag = f"pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}"
    thr = pd.read_csv(f"output/tables/csv/theta_thresholds_{tag}.csv").iloc[0]
    tex = "output/tables/tex/baseline_vs_{}_{}_taul_{}_combined.tex"
    soe = _cev_row(tex.format(f"soe_pref_tax_{int(pref_tax*100)}", SIG_TAG, τl))
    two = _cev_row(tex.format(f"two_country_pref_{int(pref_tax*100)}", SIG_TAG, τl))
    glo = _cev_row(tex.format(f"global_pref_{int(pref_tax*100)}", SIG_TAG, τl))
    # CEV columns: low-skill natives, low-skill immigrants, high-skill natives, high-skill immigrants
    pd.DataFrame([{
        "sigma": σ, "rho": config.ρ, "phi": config.φ,
        "kappa_l": config.κl_cal, "kappa_h": config.κh_cal,
        "lam_star": config.CALIBRATION.lam,
        "theta_cev0_l": thr["theta_cev0_l"], "theta_cev0_h": thr["theta_cev0_h"],
        "theta_star": thr["theta_star"],
        "ps_cev00l": soe[0], "ps_cev00h": soe[2],
        "lc_cev00l": two[0], "lc_cev00h": two[2], "lc_cev10h": two[3],
        "ml_cevl": glo[0], "ml_cevh": glo[2],
    }]).to_csv(f"output/tables/csv/sigma_theta_summary_{SIG_TAG}.csv", index=False)


def _run_sigma(σ) -> tuple:
    root = tempfile.mkdtemp(prefix="migration_sigma_theta_")
    pkg = os.path.join(root, "pkg")
    shutil.copytree(PKG, pkg, ignore=shutil.ignore_patterns(
        "output", "__pycache__", ".pytest_cache", "tests", ".DS_Store"))
    env = os.environ.copy()
    # θ only has to satisfy θ ≥ σ here: run_all's θ-dependent outputs are not used
    env.update({"MIG_SIGMA": str(σ), "MIG_THETA": str(max(3, σ)),
                "MIG_RHO": str(RHO), "MIG_PHI": str(PHI)})
    subprocess.run([sys.executable, "run_all.py"], cwd=pkg, env=env, check=True)
    subprocess.run([sys.executable, "-m", "analysis.run_sigma_theta_grid", "--cell"],
                   cwd=pkg, env=env, check=True)
    csv = os.path.join(pkg, "output", "tables", "csv")
    grid = [f for f in os.listdir(csv) if f.startswith("sigma_theta_grid_")][0]
    summ = [f for f in os.listdir(csv) if f.startswith("sigma_theta_summary_")][0]
    for f in (grid, summ):
        shutil.copy(os.path.join(csv, f), os.path.join(OUT, f))
    shutil.rmtree(root)
    return pd.read_csv(os.path.join(OUT, grid)), pd.read_csv(os.path.join(OUT, summ))


def _sig(x) -> str:
    return f"{x:g}"


def _tex_grid(grid: pd.DataFrame) -> str:
    from config import pref_tax
    sig = sorted(grid["sigma"].unique())
    p_b = grid["p10h_base"].iloc[0]
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\begin{tabular}{ c | " + " ".join("c" for _ in sig) + r" | c }",
        r"\toprule",
        r" & \multicolumn{" + str(len(sig)) + r"}{c|}{CEV, natives: low skill / high skill (\%)} & \\",
        r"$\theta$ & " + " & ".join(fr"$\sigma = {_sig(s)}$" for s in sig) + r" & $\pi_{F0}^h$ \\",
        r"\midrule",
    ]
    for θ in THETAS:
        cells = []
        for s in sig:
            r = grid[(grid["sigma"] == s) & (grid["theta"] == θ)]
            if r.empty:
                cells.append("")
            elif not bool(r["solved"].iloc[0]):
                cells.append("n.a.")
            else:
                r = r.iloc[0]
                cells.append(f"{100*r['cev00l']:+.2f} / {100*r['cev00h']:+.2f}")
        p = grid[(grid["sigma"] == SIGMA_BASE) & (grid["theta"] == θ)]["p10h"]
        lines.append(fr"{_sig(θ)} & " + " & ".join(cells) + f" & {p.iloc[0]:.3f} \\\\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{Unilateral preferential tax ($\tau^P = {pref_tax}$) in the small country, for "
        fr"values of $\sigma$ and $\theta \geq \sigma$ ($\rho = {RHO}$, $\varphi = {PHI}$). Each "
        fr"$\sigma$ is recalibrated to the same targets. $\pi_{{F0}}^h$: high-skill immigrants "
        fr"with the preferential tax at the baseline $\sigma = {_sig(SIGMA_BASE)}$ ({p_b:.3f} without it); it is "
        fr"the same for every $\sigma$ to three decimals.}}",
        r"\label{tab:sigma_theta_grid}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def _tex_summary(summ: pd.DataFrame) -> str:
    summ = summ.sort_values("sigma")
    def row(label, f):
        return label + " & " + " & ".join(f(r) for _, r in summ.iterrows()) + r" \\"
    pair = lambda a, b, fmt="{:+.1f}": (lambda r: fmt.format(r[a]) + " / " + fmt.format(r[b]))
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\begin{tabular}{ l | " + " ".join("c" for _ in range(len(summ))) + r" }",
        r"\toprule",
        " & " + " & ".join(fr"$\sigma = {_sig(s)}$" for s in summ["sigma"]) + r" \\",
        r"\midrule",
        row(r"$\kappa^l$ / $\kappa^h$", pair("kappa_l", "kappa_h", "{:.2f}")),
        row(r"$\theta$ at which low-skill natives gain", lambda r: f"{r['theta_cev0_l']:.2f}"),
        row(r"$\theta$ at which high-skill natives gain", lambda r: f"{r['theta_cev0_h']:.1f}"),
        row(r"Self-financing $\theta$", lambda r: f"{r['theta_star']:.1f}"),
        r"\midrule",
        r"\multicolumn{" + str(len(summ) + 1) + r"}{l}{CEV (\%), low skill / high skill} \\",
        row(r"Perfect substitutes, natives", pair("ps_cev00l", "ps_cev00h")),
        row(r"Large country, natives", pair("lc_cev00l", "lc_cev00h")),
        row(r"Large country, high-skill immigrants", lambda r: f"{r['lc_cev10h']:+.1f}"),
        row(r"Multilateral adoption", pair("ml_cevl", "ml_cevh")),
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{Results by $\sigma$ ($\rho = {RHO}$, $\varphi = {PHI}$). Each $\sigma$ is "
        fr"recalibrated to the same targets ($1-\lambda = {1 - summ['lam_star'].iloc[0]:.3f}$ for "
        fr"every $\sigma$). Large country: the adopting country, with the other country "
        fr"balancing its budget. Self-financing $\theta$: the preferential tax clears the "
        fr"budget at the baseline tax level.}}",
        r"\label{tab:sigma_summary}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def main():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    grids, summs = [], []
    for σ in SIGMAS:
        print(f"== sigma = {σ} (rho = {RHO}, phi = {PHI}) ==")
        g, s = _run_sigma(σ)
        grids.append(g)
        summs.append(s)
    grid, summ = pd.concat(grids), pd.concat(summs)
    grid.to_csv(os.path.join(OUT, "sigma_theta_grid.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "sigma_theta_summary.csv"), index=False)
    with open(os.path.join(OUT, "sigma_theta_grid.tex"), "w", encoding="utf-8") as f:
        f.write(_tex_grid(grid))
    with open(os.path.join(OUT, "sigma_summary.tex"), "w", encoding="utf-8") as f:
        f.write(_tex_summary(summ))
    print(f"done: {OUT}")


if __name__ == "__main__":
    cell() if "--cell" in sys.argv else main()
