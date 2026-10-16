# run_paper.py
"""
Regenerate every numerical output in the paper.

- Main run (σ, θ calibrated in model/shape_params.json: 1.275, 3.675; ρ=3, φ=0.72, ε=7, γ=0):
  run_all.py, in place; includes the θ sweep
  (app:theta, with θ* and the θ at which native CEV turns positive) and the
  perfect-substitutes (θ = inf) and large-country outputs.
- Temporary tax (app:temporary-tax): analysis.run_temporary_tax, in place.
- σ=2, θ=2 (appendix: independent shocks across destinations), ε=20
  (app:elasticity) and γ=0.05 (app:land): full run_all in a scratch copy of the
  package, collected under output/variants/{sig2,eps20,gamma005}/.
- EU-born immigrants (MIG_EU_BORN=1: immigrant share and tertiary share of the EU-born;
  φ, σ, θ recalibrated): full run_all in a scratch copy, collected under
  output/variants/eu_born/.
- σ × θ sensitivity at ρ = 3, φ = 0.72 (appendix tables tab:sigma_theta_grid,
  tab:sigma_summary): analysis.run_sigma_theta_grid, one scratch copy per σ,
  collected under output/variants/sigma_theta_grid/.
- Paper tables and figures: analysis.build_paper_tables writes the tabular block of every
  table of the draft (one file per table or panel) and analysis.build_paper_figures copies
  every figure under the draft's file names, both from the outputs above (and from the
  dynamic exercise's runs when present); destination text/tables and text/figs when the
  package sits in the paper repository, else output/paper_tables and output/paper_figures
  ($MIG_PAPER_TABLES, $MIG_PAPER_FIGS override).

Parameter variants are set through the MIG_SIGMA / MIG_THETA / MIG_EPS / MIG_GAMMA /
MIG_RHO / MIG_PHI environment variables read by config.py; defaults are the baseline calibration.
"""
import os
import shutil
import subprocess
import sys
import tempfile

PKG = os.path.dirname(os.path.abspath(__file__))


def _run(env_extra: dict, cwd: str, module: str = None):
    env = os.environ.copy()
    env.update(env_extra)
    cmd = ([sys.executable, "-m", module] if module
           else [sys.executable, "run_all.py"])
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def _scratch_copy() -> str:
    root = tempfile.mkdtemp(prefix="migration_variant_")
    pkg = os.path.join(root, "pkg")
    shutil.copytree(PKG, pkg, ignore=shutil.ignore_patterns(
        "output", "__pycache__", ".pytest_cache", "tests", ".DS_Store"))
    return pkg


def _harvest(src_pkg: str, name: str):
    dst = os.path.join(PKG, "output", "variants", name)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(os.path.join(src_pkg, "output"), dst)
    shutil.rmtree(os.path.dirname(src_pkg))


def main():
    print("== main run (calibrated sigma, theta; eps=7, gamma=0) ==")
    _run({}, PKG)

    print("== temporary tax (app:temporary-tax) ==")
    _run({}, PKG, module="analysis.run_temporary_tax")

    print("== sigma=2, theta=2 variant (appendix) ==")
    pkg = _scratch_copy()
    _run({"MIG_SIGMA": "2", "MIG_THETA": "2"}, pkg)
    _harvest(pkg, "sig2")

    print("== eps=20 variant (app:elasticity) ==")
    pkg = _scratch_copy()
    _run({"MIG_EPS": "20"}, pkg)
    _harvest(pkg, "eps20")

    print("== gamma=0.05 variant (app:land) ==")
    pkg = _scratch_copy()
    _run({"MIG_GAMMA": "0.05"}, pkg)
    _harvest(pkg, "gamma005")

    print("== EU-born immigrants variant (calibration targets) ==")
    pkg = _scratch_copy()
    _run({"MIG_EU_BORN": "1"}, pkg)
    _harvest(pkg, "eu_born")

    print("== sigma x theta sensitivity (rho=3, phi=0.72) ==")
    _run({}, PKG, module="analysis.run_sigma_theta_grid")

    print("== paper tables and figures (text/tables, text/figs or output/paper_*) ==")
    _run({}, PKG, module="analysis.build_paper_tables")
    _run({}, PKG, module="analysis.build_paper_figures")

    print("done: output/, output/variants/{sig2,eps20,gamma005,eu_born,sigma_theta_grid}/")


if __name__ == "__main__":
    main()
