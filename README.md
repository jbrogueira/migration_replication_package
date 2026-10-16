# “Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants”: Replication Package

This repository contains the code and materials to reproduce the results, figures, and tables for the project **“Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants.”** It is the **replication package** for the working paper by *[João Sousa](https://jbsousa.com)* and *[Pedro Teles](https://sites.google.com/view/pedroteles/)* (CEPR Discussion Paper **DP18414**).

## What’s inside

```
migration_replication_package/
├── analysis/                  # Scripts that solve each experiment and write tables
│   ├── run_baseline.py        # Multilateral (symmetric) Laffer curve
│   ├── run_diff_dest.py       # Small country, differentiated destinations (theta finite)
│   ├── run_soe.py, run_pref_tax.py   # Small country, perfect-substitutes limit
│   ├── run_two_country.py, run_pref_large.py, run_pref_large_ge.py  # Large country, multilateral
│   ├── run_taul_*.py          # Progressivity parameter clears the budget (appendix)
│   ├── run_theta_sweep.py     # Small country for several values of theta (appendix)
│   ├── run_temporary_tax.py   # Temporary preferential tax (appendix)
│   ├── run_g_per_capita.py    # Government spending proportional to the population, in full or in its rival share (appendix)
│   ├── run_sigma_theta_grid.py  # Sensitivity to sigma and theta
│   ├── calibrate_shape.py, shape_moments.py  # Calibration of sigma and theta (appendix)
│   ├── build_equilibrium_tables.py, exact_refine.py, run_variant.py
├── model/                     # Calibration, equations, equilibrium solvers
│   └── shape_params.json      # Calibrated sigma and theta, by specification
├── viz/                       # Plotting
├── tests/                     # pytest checks of the calibration and the solvers
├── data/eurostat/             # Eurostat extracts used for the calibration targets
│   └── eu_born/               # Extracts for the calibration with EU-born immigrants (appendix)
├── config.py                  # Parameters and calibration targets
├── run_all.py                 # Main experiments at the baseline calibration
├── run_paper.py               # Every table and figure of the static model
├── cops/                     # Dynamic multi-country exercise (appendix), see cops/README.md
├── environment.yml, requirements.txt
```

Baseline calibration (`config.py`): rho = 3, epsilon = 7, phi = 0.72, tax progressivity 0.218, preferential rate 0.114, sigma = 1.275 and theta = 3.675. The parameters sigma and theta are calibrated so that the model reproduces the pooled native logit coefficient of Kleven, Landais and Saez (2013), 1.323, and the elasticity of high-skill immigrants of Kleven et al. (2014), 2.392 (`analysis/calibrate_shape.py`); the calibrated values are stored by specification in `model/shape_params.json`, which `config.py` reads. With epsilon different from 7, phi is recalibrated to the same immigrant/native wage ratio. The remaining parameters (S, kappa, lambda, G) are calibrated in `model/calibration.py` to the data targets. The parameters sigma, theta, epsilon, gamma, rho and phi can be overridden with the environment variables `MIG_SIGMA`, `MIG_THETA`, `MIG_EPS`, `MIG_GAMMA`, `MIG_RHO` and `MIG_PHI`; `MIG_EU_BORN=1` sets the calibration targets with immigrants defined as the residents born in another EU country.

> All tables (CSV and LaTeX) and figures are written to a new `output/` folder.

## Reference working paper

Please cite the paper when using this package:

> **Brogueira de Sousa, J. and Teles, P.** “Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants.” CEPR Discussion Paper **DP18414**.

## Getting started

### 1) Environment (conda / mamba)

This project ships an **`environment.yml`**. It pins **Python 3.11** and the core dependencies (numpy, scipy, matplotlib, pandas). Create the environment with **conda**:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate migration_replication_package
```

### 2) Configuration (optional)

Open `config.py` to adjust paths or calibration options. Defaults should work out of the box.

## Reproducing everything (figures + tables)

From the repo root (with the environment activated):

```bash
python run_paper.py
```

This runs:

- `run_all.py` at the baseline calibration: the multilateral and unilateral Laffer curves, the small country with differentiated destinations and in the perfect-substitutes limit, the large country, multilateral adoption, the progressivity variant, government spending proportional to the population, in full and in its rival share 0.62 (`analysis.run_g_per_capita`, `--share`), and the theta sweep (`output/tables/`, `output/figures/`);
- the temporary preferential tax (`analysis.run_temporary_tax`);
- the variants sigma = theta = 2, epsilon = 20, land (gamma = 0.05) and EU-born immigrants (`MIG_EU_BORN=1`), each in a temporary copy of the package, collected in `output/variants/{sig2,eps20,gamma005,eu_born}/`; in the epsilon = 20 variant the theta sweep solves the perfect-substitutes limit on the decreasing side of the Laffer curve, where the budget clears (`solve_high_tax_branch_ps` in `analysis/run_theta_sweep.py`);
- the sensitivity to sigma and theta (`analysis.run_sigma_theta_grid`), in `output/variants/sigma_theta_grid/`.

The dynamic exercise in the appendix, based on the model of Caliendo, Opromolla, Parro and Sforza (2021), is in `cops/` (Octave or MATLAB); `cops/README.md` lists the files to take from their replication package and the steps to run it. The full run of `run_paper.py` takes about 6.5 minutes on an Apple M2 laptop. Output file names carry the values of sigma and theta (for example `sig_1p275_theta_3p675`). Every reported equilibrium satisfies the government budget constraint to numerical precision.

To run only the main experiments:

```bash
python run_all.py
```

To run the tests (a few seconds; one test reads an output of `run_all.py` and is skipped without it):

```bash
python -m pytest -q
```

The calibrated sigma and theta are stored in `model/shape_params.json`, so the steps above do not recompute them. To recompute them for a specification:

```bash
python -m analysis.calibrate_shape                  # baseline
MIG_EPS=20 python -m analysis.calibrate_shape       # epsilon = 20
MIG_GAMMA=0.05 python -m analysis.calibrate_shape   # land
MIG_EU_BORN=1 python -m analysis.calibrate_shape    # EU-born immigrants
```

## How to cite

If you use this code or reproduce results, please cite both the paper and this repository. Example:

```
Brogueira de Sousa, J. and Teles, P. (2023). "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants". CEPR DP18414.
Replication package: "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants: Replication Package".
```

## Contact

- Maintainer: **[Joao Brogueira de Sousa](https://jbsousa.com)** (Nova SBE) – joao.sousa@novasbe.pt  
- Issues and questions: please open a GitHub issue or send an email.
