# “Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants”: Replication Package

This repository contains the code and materials to reproduce the results, figures, and tables for the project **“Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants.”** It is the **replication package** for the working paper by *[João Sousa](https://jbsousa.com)* and *[Pedro Teles](https://sites.google.com/view/pedroteles/)* (CEPR Discussion Paper **DP18414**).

## What’s inside

```
migration_replication_package/
├── analysis/                  # Run scripts that generate figures/tables
│   ├── build_equilibrium_tables.py
│   ├── run_baseline.py
│   ├── run_pref_large.py
│   ├── run_pref_large_ge.py
│   ├── run_pref_tax.py
│   ├── run_soe.py
│   ├── run_taul_soe.py
│   ├── run_taul_sweeps.py
│   └── run_two_country.py
├── model/                     # Core model code (equations, equilibrium, parameters)
│   ├── equations.py
│   ├── equilibrium.py
│   ├── parameters.py
│   └── __init__.py
├── viz/                       # Plotting helpers
├── config.py                  # Global configuration (paths, toggles, calibration hooks)
├── run_all.py                 # Entry point: runs the full replication pipeline
├── environment.yml            # Conda environment spec (Python 3.11, numpy, scipy, matplotlib, pandas)
└── read_me.txt                # Notes / quick reminders
```

> By design, **all charts and tables** are created in a new `output/` folder (auto-created if missing).

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
python run_all.py
```

This will:

- run the baseline and policy counterfactuals (small vs. large economy; unilateral vs. multilateral; alternative elasticities if configured),
- build tables summarizing migration, wages, taxes, and welfare (CEV),
- save **all outputs** under `./output/` (plots and CSV/Excel tables, depending on your settings).

## Running components individually

Execute specific parts as needed:

```bash
python -m analysis.run_baseline
python -m analysis.run_soe                 # Small open economy experiments
python -m analysis.run_pref_tax            # Preferential tax for high-skill immigrants
python -m analysis.run_pref_large          # Large-country unilateral adoption
python -m analysis.run_pref_large_ge       # Large-country general equilibrium variants
python -m analysis.run_two_country         # Two-country model
python -m analysis.build_equilibrium_tables
```

(Each script writes its own artifacts into `output/`.)

## Expected outputs

- **Figures**: Laffer curves; consumption/welfare comparisons (CEV) by type and scenario; wage/tax incidence visuals.  
- **Tables**: Baseline calibration, equilibrium allocations, migration shares, wages, tax rates, and welfare for each scenario.

## How to cite

If you use this code or reproduce results, please cite both the paper and this repository. Example:

```
Brogueira de Sousa, J. and Teles, P. (2023). "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants". CEPR DP18414.
Replication package: "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants: Replication Package".
```

## Contact

- Maintainer: **[Joao Brogueira de Sousa](https://jbsousa.com)** (Nova SBE) – joao.sousa@novasbe.pt  
- Issues and questions: please open a GitHub issue or send an email.
