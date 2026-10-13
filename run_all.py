# run_all.py
import os
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))

from config import σ, θ, G0, pref_tax, τl, γ, SIG_TAG, DD_TAG, g_rival_share
from model.equilibrium import pi_ge_eq, tau10h_star, tau_star
from analysis.run_baseline import main as run_base
from analysis.run_soe import main as run_soe
from analysis.run_pref_tax import main as run_pref_tax
from analysis.run_diff_dest import main as run_diff_dest
from analysis.run_g_per_capita import main as run_g_per_capita
from analysis.run_two_country import main as run_two_country_pe
from analysis.run_pref_large import main as run_pref_large
from analysis.run_pref_large_ge import main as run_pref_large_ge
from analysis.run_taul_sweeps import main as run_taul_sweeps
from analysis.run_taul_soe import main as run_taul_soe
from analysis.build_equilibrium_tables import (
    build_equilibrium_tables,
    build_equilibrium_tables1,
    build_equilibrium_tables_diff_dest,
    build_baseline_model_vars_table,
)
from analysis.run_taul_diff_dest import main as run_taul_diff_dest
from analysis.run_theta_sweep import main as run_theta_sweep
from viz.plotting import (
    plot_tax_revenue,
    plot_comparison,
    plot_diff_dest_laffer,
    plot_unilateral_laffer,
    plot_tax_vs_tau10h,
    plot_two_country_comparison,
    plot_two_country_unilateral,
    plot_tax_rates_vs_lambda,
    plot_net_incomes_vs_lambda,
    plot_taul_mlt_vs_soe,
    plot_taul_soe_pref,
    plot_two_country_taul_case2,
    plot_two_country_taul_case3,
)


# Define output paths
OUTPUT_DIR   = os.path.join(os.getcwd(), "output")
FIGURES_DIR  = os.path.join(OUTPUT_DIR, "figures")
TABLES_DIR   = os.path.join(OUTPUT_DIR, "tables")


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)
    
    # Run all analysis scripts, which handle saving CSV & TeX themselves
    df_base = run_base()
    df_soe = run_soe()
    df_pref, df_tau10h = run_pref_tax()
    _, _, df_dd_base, df_dd_pref, df_dd_cross = run_diff_dest()
    run_g_per_capita()   # spending proportional to the population: table + figure
    run_g_per_capita(share=g_rival_share)   # rival share of spending proportional to the population
    df_two_country_pe = run_two_country_pe()
    df_pref_large = run_pref_large()
    df_pref_large_gbc = run_pref_large_ge()
    df_tau_pe, df_tau_ge = run_taul_sweeps()
    df_tau_soe_pref = run_taul_soe()

    # Generate figures
    # - Baseline
    plot_tax_revenue(
        df_base, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"base_{SIG_TAG}_taul_{τl}.png")
    )
    # - Multilateral vs the two unilateral cases (perfect substitutes,
    #   small country at the active θ)
    plot_unilateral_laffer(
        df_base, df_soe, df_dd_base, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"comparison_{DD_TAG}_taul_{τl}.png"),
        dd_label=fr"Unilateral, $\theta = {θ:g}$",
    )
    # - SOE vs Preferential tax
    plot_comparison(
        df_soe, df_pref, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"pref_t_{SIG_TAG}_taul_{τl}.png")
    )
    
    # - Differentiated destinations: revenue with and without the preferential tax,
    #   SOE preferential curve overlaid (perfect-substitutes limit)
    plot_diff_dest_laffer(
        df_dd_base, df_dd_pref, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"diff_dest_laffer_{DD_TAG}_taul_{τl}.png"),
        df_soe_pref=df_pref,
    )

    # - SOE vs Preferential tax (Marked equilibrium)
    # plot_comparison(df_soe, df_pref, tau_star, G0,
    #             filename=os.path.join(FIGURES_DIR, f"pref_t_{SIG_TAG}_taul_{τl}_marked.png"), mark_equilibrium=True,eq_kwargs={'xytext': (-0.05, 0.36)})
    
    
    # - Tax revenue vs high-skill tax
    # plot_tax_vs_tau10h(
    #     df_tau10h, G0,
    #     tau10h_star=tau10h_star,
    #     filename=os.path.join(FIGURES_DIR, f"tau10h_sweep_{SIG_TAG}_taul_{τl}.png")
    # )

    # Large economy (unilateral preferential tax)
    plot_two_country_unilateral(
        df_pref_large,
        pi_ge_eq,
        G0,
        filename=os.path.join(FIGURES_DIR, f"two_country_pref_{SIG_TAG}_taul_{τl}.png")
        )
    
    # Large economy (global preferential tax)
    plot_two_country_comparison(
        df_pref_large_gbc,
        df_pref_large,
        pi_ge_eq,
        G0,
        filename=os.path.join(FIGURES_DIR, f"global_pref_{SIG_TAG}_taul_{τl}.png")
    )
    
    # Compute paper tables
    # Baseline vs SOE preferential tax 
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    )
    
    # Baseline vs differentiated-destinations preferential tax (both rows exact)
    build_equilibrium_tables_diff_dest()

    # Baseline vs two-country preferential
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    )
    
    # Baseline vs two-country preferential (country 1)
    build_equilibrium_tables1(
        pref_eq_path=f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    )
    
    # Baseline vs large-GE budget constraint 
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/global_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    )
    
    
    # Sweeping tau instead of sweeping lambda for appendix exercise (no land income)

    if γ == 0:
        # SOE (λ fixed, preferential tax in country 0)
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/case1_soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        )

        # Differentiated destinations (λ fixed at λ*, τl clears the budget)
        _, taul_dd_csv = run_taul_diff_dest()
        build_equilibrium_tables_diff_dest(pref_eq_path=taul_dd_csv)
    
        # Two large economies — Case 2: pref only in Country 0
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/two_country_taul_case2_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        )
    
        # Two large economies — Case 3: bilateral pref (both countries)
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/two_country_taul_case3_pref_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
        )
        
        
   
    # Small country for several θ (appendix app:theta)
    if γ == 0:
        run_theta_sweep()

    # Build Table 2 of the paper
    build_baseline_model_vars_table()
    



if __name__ == "__main__":
    main()
