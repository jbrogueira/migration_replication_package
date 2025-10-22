# run_all.py
import os
import pandas as pd

os.chdir(r"C:\Users\csi\Documents\PhD\RA_Migration\migration\migration_project")

from config import σ, G0, pref_tax, τl, γ
from model.equilibrium import pi_ge_eq, tau10h_star, tau_star
from analysis.run_baseline import main as run_base
from analysis.run_soe import main as run_soe
from analysis.run_pref_tax import main as run_pref_tax
from analysis.run_two_country import main as run_two_country_pe
from analysis.run_pref_large import main as run_pref_large
from analysis.run_pref_large_ge import main as run_pref_large_ge
from analysis.run_taul_sweeps import main as run_taul_sweeps
from analysis.run_taul_soe import main as run_taul_soe
from analysis.build_equilibrium_tables import build_equilibrium_tables, build_equilibrium_tables1

from viz.plotting import plot_tax_revenue, plot_comparison, plot_tax_vs_tau10h, plot_two_country_comparison, plot_two_country_unilateral, plot_tax_rates_vs_lambda, plot_net_incomes_vs_lambda
from viz.plotting import plot_taul_mlt_vs_soe, plot_taul_soe_pref, plot_two_country_taul_case2, plot_two_country_taul_case3


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
    df_two_country_pe = run_two_country_pe()
    df_pref_large = run_pref_large()
    df_pref_large_gbc = run_pref_large_ge()
    if γ == 0:
        df_tau_pe, df_tau_ge = run_taul_sweeps()
        df_tau_soe_pref = run_taul_soe()
        return df_tau_pe, df_tau_ge, df_tau_soe_pref


    
    # Compute paper tables
    # Baseline vs SOE preferential tax 
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    )
    
    # Baseline vs two-country preferential 
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    )
    
    build_equilibrium_tables1(
        pref_eq_path=f"output/tables/csv/two_country_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    )
    
    # Baseline vs large-GE budget constraint 
    build_equilibrium_tables(
        pref_eq_path=f"output/tables/csv/global_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    )
    
    # === τ_l sweep tables (no Country 1 version, no captions/labels) ===

    if γ == 0:
        # SOE (λ fixed, preferential tax in country 0)
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/case1_soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
        )

        # Two large economies — Case 2: pref only in Country 0
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/two_country_taul_case2_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
        )

        # Two large economies — Case 3: bilateral pref (both countries)
        build_equilibrium_tables(
            pref_eq_path=f"output/tables/csv/two_country_taul_case3_pref_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
        )




    # Generate figures
    # - Baseline
    plot_tax_revenue(
        df_base, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"base_sig_{int(σ)}_taul_{τl}.png")
    )
    # - SOE vs Symmetric
    plot_comparison(
        df_base, df_soe, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"comparison_sig_{int(σ)}_taul_{τl}.png")
    )
    # - SOE vs Preferential tax
    plot_comparison(
        df_soe, df_pref, tau_star, G0,
        filename=os.path.join(FIGURES_DIR, f"pref_t_sig_{int(σ)}_taul_{τl}.png")
    )
    
    plot_comparison(df_soe, df_pref, tau_star, G0,
                filename=os.path.join(FIGURES_DIR, f"pref_t_sig_{int(σ)}_taul_{τl}_marked.png"), mark_equilibrium=True,eq_kwargs={'xytext': (-0.05, 0.36)})
    
    # - Tax revenue vs high-skill tax
    plot_tax_vs_tau10h(
        df_tau10h, G0,
        tau10h_star=tau10h_star,
        filename=os.path.join(FIGURES_DIR, f"tau10h_sweep_sig_{int(σ)}_taul_{τl}.png")
    )


    # plot_two_country_comparison(
    #     df_two_country_pe,
    #     df_base,               # symmetric baseline
    #     pi_ge_eq,
    #     G0,
    #     filename=os.path.join(FIGURES_DIR, f"two_country_pe_sig_{int(σ)}_taul_{τl}.png")
    # )
    
    
    plot_two_country_unilateral(
        df_pref_large,
        pi_ge_eq,
        G0,
        filename=os.path.join(FIGURES_DIR, f"two_country_pref_sig_{int(σ)}_taul_{τl}.png")
        )
    
    
    plot_two_country_comparison(
        df_pref_large_gbc,  # PE sweep with TR0, TR1
        df_pref_large,      # symmetric baseline with TR column
        pi_ge_eq,
        G0,
        filename=os.path.join(FIGURES_DIR, f"global_pref_sig_{int(σ)}_taul_{τl}.png")
    )
    



if __name__ == "__main__":
    main()
