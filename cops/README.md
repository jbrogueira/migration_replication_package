# Preferential tax in the dynamic multi-country model (appendix)

This folder computes the transition induced by the preferential tax in the dynamic model of
Caliendo, Opromolla, Parro and Sforza (2021, *Journal of Political Economy* 129(12),
"Goods and Factor Market Integration: A Quantitative Assessment of the EU Enlargement";
COPS below), for Portugal, Spain and Germany. It produces the appendix tables
`tab:cops-transition` and `tab:cops-cev` and the figures `fig:cops-stocks`, `fig:cops-taxes`,
`fig:cops-welfare` and `fig:cops-homotopy`.

The folder contains only our code. The COPS functions and data it uses are not
redistributed here; they come from the COPS replication package published with the article.
Requirements: GNU Octave (tested with 9.4) or MATLAB, and Python 3 with numpy, scipy and matplotlib.

## Our files

| File | Role |
| --- | --- |
| `run_cops.m` | Driver: zero-policy check, open-budget solve at tau_P = 0.2127, closed-budget homotopy 0.2127 → 0.19 → 0.16 → 0.13 → 0.114, stock paths at each step and welfare paths at 0.114 |
| `anderson_solve.m`, `sweep_phi.m` | Solver of the transition (Anderson acceleration on the value hats) |
| `sweep_engine.m` | One pass of the dynamic hat algebra with the progressive tax, the preferential schedule, temporary eligibility and period-by-period budget clearing |
| `dump_pt_stocks.m`, `compute_cev_paths.m` | Stocks of the adopting country by skill and tax status; consumption-equivalent variation from the value functions |
| `summarize_cops_run.py` | Numbers of the two appendix tables |
| `plot_transition.py` | Figures |

The baseline average tax rates (0.2436 for high-skill and 0.1668 for low-skill workers) are
the native rates of the static model at its baseline calibration (`baseline_*_equilibrium.csv`
of `run_paper.py`); `run_cops.m` sets them. The preferential rate of the paper, 0.114, is half
the baseline average rate of high-skill immigrants in the static model (0.2284); the
homotopy reaches it in steps from 0.2127.

## Steps

1. **COPS files.** From the folder `Benchmark_results_EU_enlargement` of the COPS replication
   package, copy into this folder the functions `Dinprime_eu.m`, `GMCeu_skill.m`,
   `GMCeu_unskill.m`, `P_h_om_eu_skill.m`, `expenditure_eu_skill.m`,
   `expenditure_eu_unskilled.m`, `tempeq_skill.m` and the data files `Base_year2002_skill.mat`,
   `Baseline_europe_data_skills.mat`, `data_xbilat_tariffs_tau_skills.mat`. We use them unchanged.

2. **Baseline without public goods.** In the COPS folder `Public_goods/step3`, open
   `Step_3_Baseline_forward_europe_skills.m` and replace the vector assigned to `alpha`
   (the shares of public goods, just before the lines that compute `realwages_h` and
   `realwages_l`) by `alpha=0;`. Run the script, copy the resulting
   `Baseline_europe_forward_skills.mat` into `Public_goods/step4`, and run
   `Step_4_europe_skills.m` there. Copy the resulting `Baseline_europe_Factual_skills.mat`
   into this folder. We checked that this reproduces our baseline file exactly (up to the
   timestamp in the file header).

3. **Transition.** From this folder:

   ```bash
   mkdir runs && cd runs
   octave --no-gui --eval "addpath('..'); run_cops('pt')"
   octave --no-gui --eval "addpath('..'); run_cops('es')"
   octave --no-gui --eval "addpath('..'); run_cops('de')"
   ```

   Each country took two to three hours on our machine, with the three running in parallel.
   Each run stops with an error if the zero-policy check does not reproduce the baseline or
   if a step of the homotopy does not converge.

4. **Tables and figures.** From this folder:

   ```bash
   python summarize_cops_run.py --dir runs --prefix v6 --taup 0.114
   python plot_transition.py --dir runs --prefix v6_pt --taup 0.114 --tau-h 0.2436 --tau-l 0.1668 --out figures
   ```

   The paper uses `figures/fig2_response_composition.png` (stocks),
   `fig4_taxes_budget.png` (taxes), `fig6_welfare_cev.png` (welfare) and `fig7_homotopy.png`.
