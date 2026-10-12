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
│   ├── run_g_per_capita.py    # Government spending proportional to the population (appendix)
│   ├── run_sigma_theta_grid.py  # Sensitivity to sigma and theta
│   ├── build_equilibrium_tables.py, exact_refine.py, run_variant.py
├── model/                     # Calibration, equations, equilibrium solvers
├── viz/                       # Plotting
├── tests/                     # pytest checks of the calibration and the solvers
├── data/eurostat/             # Eurostat extracts used for the calibration targets
├── config.py                  # Parameters and calibration targets
├── run_all.py                 # Main experiments at the baseline calibration
├── run_paper.py               # Every table and figure of the static model
├── cops/                     # Dynamic multi-country exercise (appendix), see cops/README.md
├── environment.yml, requirements.txt
```

Baseline calibration (`config.py`): sigma = 1.323, theta = 3.202, rho = 3, epsilon = 7, phi = 0.72, tax progressivity 0.216, preferential rate 0.13. The remaining parameters (S, kappa, lambda, G) are calibrated in `model/calibration.py` to the data targets. The parameters sigma, theta, epsilon, gamma, rho and phi can be overridden with the environment variables `MIG_SIGMA`, `MIG_THETA`, `MIG_EPS`, `MIG_GAMMA`, `MIG_RHO` and `MIG_PHI`.

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

- `run_all.py` at the baseline calibration: the multilateral and unilateral Laffer curves, the small country with differentiated destinations and in the perfect-substitutes limit, the large country, multilateral adoption, the progressivity variant, government spending proportional to the population (`analysis.run_g_per_capita`), and the theta sweep (`output/tables/`, `output/figures/`);
- the temporary preferential tax (`analysis.run_temporary_tax`);
- the variants sigma = theta = 2, epsilon = 20 and land (gamma = 0.05), each in a temporary copy of the package, collected in `output/variants/{sig2,eps20,gamma005}/`;
- the sensitivity to sigma and theta (`analysis.run_sigma_theta_grid`), in `output/variants/sigma_theta_grid/`.

The dynamic exercise in the appendix, based on the model of Caliendo, Opromolla, Parro and Sforza (2021), is in `cops/` (Octave or MATLAB); `cops/README.md` lists the files to take from their replication package and the steps to run it. The full run takes a few minutes. Output file names carry the values of sigma and theta (for example `sig_1p323_theta_3p202`). Every reported equilibrium satisfies the government budget constraint to numerical precision.

To run only the main experiments:

```bash
python run_all.py
```

To run the tests:

```bash
python -m pytest -q
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
