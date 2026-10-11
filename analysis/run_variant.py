# analysis/run_variant.py
# One parameter variant of the small-country exercises: multilateral baseline,
# unilateral perfect-substitutes (SOE) and differentiated-destinations runs, the
# three-curve unilateral Laffer figure, and the SOE + diff-dest tables.
# The variant is set by the MIG_SIGMA / MIG_THETA / MIG_EPS / MIG_GAMMA environment
# variables read by config.py; all output filenames embed the active σ, θ (and τl).
import os

import matplotlib
matplotlib.use("Agg")

from config import σ, θ, τl, G0, SIG_TAG, DD_TAG
from model.equilibrium import tau_star
from analysis.run_baseline import main as run_base
from analysis.run_soe import main as run_soe
from analysis.run_pref_tax import main as run_pref_tax
from analysis.run_diff_dest import main as run_diff_dest
from analysis.build_equilibrium_tables import (
    build_equilibrium_tables,
    build_equilibrium_tables_diff_dest,
)
from viz.plotting import (
    plot_comparison,
    plot_diff_dest_laffer,
    plot_unilateral_laffer,
)

FIGURES_DIR = "output/figures"


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    df_base = run_base()
    df_soe = run_soe()
    df_pref, _ = run_pref_tax()
    _, _, df_dd_base, df_dd_pref, _ = run_diff_dest()

    # fig 2-family: multilateral vs the two unilateral cases
    plot_unilateral_laffer(
        df_base, df_soe, df_dd_base, tau_star, G0,
        filename=f"{FIGURES_DIR}/comparison_{DD_TAG}_taul_{τl}.png",
        dd_label=fr"Unilateral, $\theta = {θ:g}$")
    # SOE preferential-tax Laffer (perfect substitutes)
    plot_comparison(
        df_soe, df_pref, tau_star, G0,
        filename=f"{FIGURES_DIR}/pref_t_{SIG_TAG}_taul_{τl}.png")
    # differentiated-destinations Laffer with the perfect-substitutes overlay
    plot_diff_dest_laffer(
        df_dd_base, df_dd_pref, tau_star, G0,
        filename=f"{FIGURES_DIR}/diff_dest_laffer_{DD_TAG}_taul_{τl}.png",
        df_soe_pref=df_pref)

    build_equilibrium_tables()
    build_equilibrium_tables_diff_dest()


if __name__ == "__main__":
    main()
