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
│   ├── build_paper_tables.py, build_paper_figures.py  # The paper's table files and figures
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
- the sensitivity to sigma and theta (`analysis.run_sigma_theta_grid`), in `output/variants/sigma_theta_grid/`;
- the paper's tables and figures (`analysis.build_paper_tables`, `analysis.build_paper_figures`), from the outputs above, in `output/paper_tables/` and `output/paper_figures/` (see "Tables and figures of the paper" below).

The dynamic exercise in the appendix, based on the model of Caliendo, Opromolla, Parro and Sforza (2021), is in `cops/` (Octave or MATLAB); `cops/README.md` lists the files to take from their replication package and the steps to run it. The full run of `run_paper.py` takes about 6.5 minutes on an Apple M2 laptop. Output file names carry the values of sigma and theta (for example `sig_1p275_theta_3p675`). Every reported equilibrium satisfies the government budget constraint to numerical precision.

To run only the main experiments:

```bash
python run_all.py
```

To run the tests (a few seconds; one test reads an output of `run_all.py`, and the tests of the paper's tables and figures read the outputs of `run_paper.py`; they are skipped without them):

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

## Tables and figures of the paper

`python -m analysis.build_paper_tables` writes, for every table of the paper, a file with its `tabular` block(s), exactly as the paper inputs it (`\input{tables/<name>}`); captions and notes are in the paper. A table with two panels has one file per panel. `python -m analysis.build_paper_figures` copies every figure of the paper under the file name the paper uses. Both read only existing outputs, so they run in seconds after `run_paper.py` (which calls them at the end). The destinations are `output/paper_tables/` and `output/paper_figures/` (environment variables `MIG_PAPER_TABLES` and `MIG_PAPER_FIGS` override them). The two tables and four figures of the dynamic exercise are written only when the runs of `cops/` are in `cops/runs/` and their figures in `cops/figures/` (`cops/README.md`, steps 3 and 4); otherwise the scripts say so and skip them.

| File | Paper table | Source |
|---|---|---|
| `tab_model_param.tex` | `tab:model-param` | `config.py`, `model/calibration.py` |
| `tab_diff_dest.tex` | `tab:diff_dest_pref_tax_13_sig_2_taul_0.216_vs_baseline` (small country) | `baseline_vs_diff_dest_pref_tax_11_*_combined.tex` |
| `tab_two_country_a.tex`, `tab_two_country_b.tex` | `tab:two_country_pref_16_sig_2_large_vs_baseline` (large country, panels A and B) | `baseline_vs_two_country_pref_11_*_combined{,_cty1}.tex` |
| `tab_multilateral.tex` | `tab:pref_16_sig_2_large_gbc_vs_baseline` (multilateral) | `baseline_vs_global_pref_11_*_combined.tex` |
| `tab_welfare_summary.tex` | `tab:baseline_welfare_comparison` | CEV rows of the small country, theta sweep (theta = 5, 10), perfect-substitutes limit, large country, multilateral |
| `sigma_theta_grid.tex`, `sigma_summary.tex` | `tab:sigma_theta_grid`, `tab:sigma_summary` | `output/variants/sigma_theta_grid/` (sigma labels at two decimals, the baseline theta at one) |
| `tab_ps_robust.tex` | `tab:ps-robust` | `baseline_vs_soe_pref_tax_11_*` and `baseline_vs_grival_0p62_soe_pref_tax_11_*` |
| `tab_nested_elasticities.tex` | `tab:nested-elasticities` | formulas only |
| `tab_education.tex` | `tab:education` | `data/eurostat/edat_lfs_9915_2024_paper_release.csv` (see the note below) |
| `tab_population.tex` | `tab:population` | `data/eurostat/eu_born/migr_pop3ctb_EU27_2024.json`, `config.py` |
| `tab_earnings.tex` | `tab:earnings` | `config.py` (Structure of Earnings Survey 2022) |
| `tab_theta_sweep.tex` | `tab:theta_sweep_sig_1p323` | `theta_sweep_pref_tax_11_*.tex` |
| `tab_case1_soe.tex` | `tab:case1_soe_pref_tax_13_sig_2_taul_0.216_vs_baseline` | `baseline_vs_case1_soe_pref_tax_11_*` |
| `tab_diff_dest_taul.tex` | `tab:diff_dest_taul_theta_2` | `baseline_vs_diff_dest_taul_pref_tax_11_*` |
| `tab_two_country_taul.tex` | `tab:two_country_taul_case2_pref_13_sig_2_taul_0.216_vs_baseline` | `baseline_vs_two_country_taul_case2_pref_11_*` |
| `tab_multilateral_taul.tex` | `tab:global_pref_13_sig_2_taul_0.216_vs_baseline_taul` | `baseline_vs_two_country_taul_case3_pref_11_*` |
| `tab_gpc_diff_dest.tex`, `tab_gpc_ps.tex` | `tab:g-per-capita`, `tab:g-per-capita-ps` | `baseline_vs_diff_dest_gpc_*`, `baseline_vs_gpc_soe_*` |
| `tab_gpc_two_country_a.tex`, `_b.tex` | `tab:g-per-capita-lc` (panels A and B) | `baseline_vs_gpc_two_country_pref_11_*` |
| `tab_grival_diff_dest.tex` | `tab:g-rival` | `baseline_vs_diff_dest_grival_0p62_*` |
| `tab_grival_two_country_a.tex`, `_b.tex` | `tab:g-rival-lc` (panels A and B) | `baseline_vs_grival_0p62_two_country_pref_11_*` |
| `tab_model_param_eps20.tex` | `tab:model-param-eps-20` | `config.py` with `MIG_EPS=20` |
| `tab_baseline_vars_eps20.tex` | `tab:baseline_sig_2_taul_0.216_baseline_model_vars_i0` | `variants/eps20/.../baseline_model_vars_i0.tex` |
| `tab_two_country_eps20.tex`, `tab_multilateral_eps20.tex`, `tab_diff_dest_eps20.tex`, `tab_theta_sweep_eps20.tex` | `tab:two_country_pref_13_sig_2_taul_0.216_vs_baseline`, `tab:global_pref_13_sig_2_taul_0.216_vs_baseline_eps20`, `tab:diff_dest_eps20`, `tab:theta_sweep_eps20` | `output/variants/eps20/` |
| `tab_diff_dest_land.tex`, `tab_ps_land.tex`, `tab_two_country_land.tex`, `tab_multilateral_land.tex`, `tab_welfare_summary_land.tex` | `tab:diff_dest_land`, `tab:soe_pref_tax_13_sig_2_taul_0.216_vs_baseline_gamma`, `tab:two_country_pref_13_sig_2_taul_0.216_vs_baseline_gamma`, `tab:global_pref_13_sig_2_taul_0.216_vs_baseline_gamma`, `tab:land_welfare_comparison` | `output/variants/gamma005/` |
| `tab_temporary_tax.tex` | `tab:temporary-tax` | `analysis.run_temporary_tax` outputs (`pref_tax_{11,20,18}`) |
| `tab_cops_transition.tex`, `tab_cops_cev.tex` | `tab:cops-transition`, `tab:cops-cev` | `cops/summarize_cops_run.py` applied to `cops/runs/` (CEV in exact percent) |
| `tab_eu_born.tex` | `tab:eu-born` | equilibrium CSVs and CEV rows of the main run and of `output/variants/eu_born/` |

Figures: `base_*`, `comparison_*`, `pref_t_*`, `two_country_pref_*`, `global_pref_*`, `theta_laffer_*`, `diff_dest_laffer_*`, `g_per_capita_laffer_*`, `g_rival_0p62_laffer_*` from `output/figures/`; the land figures `*_sig_1p38_taul_0.218_gamma.png` from `output/variants/gamma005/figures/`; `cops_stocks.png`, `cops_taxes.png`, `cops_welfare.png`, `cops_homotopy.png` from `cops/figures/` (`fig2_response_composition.png`, `fig4_taxes_budget.png`, `fig6_welfare_cev.png`, `fig7_homotopy.png`).

The cells of `tab:education` other than the two tertiary shares used in the calibration are not in any extract of the package: they are from an earlier release of Eurostat `edat_lfs_9915` that the API no longer serves, and `data/eurostat/edat_lfs_9915_2024_paper_release.csv` records them with their source (`data/eurostat/eu_born/README.md` compares them with the current release).

## How to cite

If you use this code or reproduce results, please cite both the paper and this repository. Example:

```
Brogueira de Sousa, J. and Teles, P. (2023). "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants". CEPR DP18414.
Replication package: "Taxes and Migration Flows: Preferential Tax Schemes for High-Skill Immigrants: Replication Package".
```

## Contact

- Maintainer: **[Joao Brogueira de Sousa](https://jbsousa.com)** (Nova SBE) – joao.sousa@novasbe.pt  
- Issues and questions: please open a GitHub issue or send an email.
