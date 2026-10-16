# analysis/build_paper_figures.py
"""
Copy the figures of the paper draft (text/main_final.tex, text/appendix_final.tex) from the
outputs of run_paper.py, under the file names the draft's \\includegraphics use.

Sources: output/figures/ (main run), output/variants/gamma005/figures/ (land), and the
figures of the dynamic exercise (Appendix M) written by plot_transition.py; the latter are
skipped with a message when absent.

Destination: $MIG_PAPER_FIGS; default ../../text/figs relative to the package root if that
directory exists, else output/paper_figures/.

Usage (from the package root): python -m analysis.build_paper_figures
"""
import os
import shutil

from analysis.build_paper_tables import MAIN as _M, LAND as _L, TAUL, PREF_TAG, OUT, PKG

# file tags: σ (and θ) of each specification, from model/shape_params.json
MAIN, DD = f"{_M.sig}_taul_{TAUL}", f"{_M.dd}_taul_{TAUL}"
LAND = f"{_L.sig}_taul_{TAUL}"

# paper file name -> source, relative to output/
STATIC = {
    f"base_{MAIN}.png": f"figures/base_{MAIN}.png",                              # fig:multilateral-laffer
    f"comparison_{DD}.png": f"figures/comparison_{DD}.png",                      # fig:unilateral-laffer
    f"pref_t_{MAIN}.png": f"figures/pref_t_{MAIN}.png",                          # fig:laffer-soe-noem
    f"two_country_pref_{MAIN}.png": f"figures/two_country_pref_{MAIN}.png",      # fig:laffer-lc
    f"global_pref_{MAIN}.png": f"figures/global_pref_{MAIN}.png",                # fig:laffer-lc-multi
    f"theta_laffer_{MAIN}.png": f"figures/theta_laffer_{MAIN}.png",              # fig:laffer-theta
    f"diff_dest_laffer_{DD}.png": f"figures/diff_dest_laffer_{DD}.png",          # fig:laffer-diff-dest
    f"g_per_capita_laffer_pref_tax_{PREF_TAG}_{DD}.png": f"figures/g_per_capita_laffer_pref_tax_{PREF_TAG}_{DD}.png",
    f"g_per_capita_laffer_ps_pref_tax_{PREF_TAG}_{MAIN}.png": f"figures/g_per_capita_laffer_ps_pref_tax_{PREF_TAG}_{MAIN}.png",
    f"g_per_capita_laffer_lc_pref_tax_{PREF_TAG}_{MAIN}.png": f"figures/g_per_capita_laffer_lc_pref_tax_{PREF_TAG}_{MAIN}.png",
    f"g_rival_0p62_laffer_pref_tax_{PREF_TAG}_{DD}.png": f"figures/g_rival_0p62_laffer_pref_tax_{PREF_TAG}_{DD}.png",
    f"g_rival_0p62_laffer_ps_pref_tax_{PREF_TAG}_{MAIN}.png": f"figures/g_rival_0p62_laffer_ps_pref_tax_{PREF_TAG}_{MAIN}.png",
    f"g_rival_0p62_laffer_lc_pref_tax_{PREF_TAG}_{MAIN}.png": f"figures/g_rival_0p62_laffer_lc_pref_tax_{PREF_TAG}_{MAIN}.png",
    # land (app:land): the variant's figures, renamed with the suffix _gamma
    f"pref_t_{LAND}_gamma.png": f"variants/gamma005/figures/pref_t_{LAND}.png",                    # fig:laffer-sc-land
    f"two_country_pref_{LAND}_gamma.png": f"variants/gamma005/figures/two_country_pref_{LAND}.png",  # fig:laffer-lc-land
    f"global_pref_{LAND}_gamma.png": f"variants/gamma005/figures/global_pref_{LAND}.png",          # fig:laffer-multi-land
}

# dynamic exercise: plot_transition.py --prefix v8_pt --taup 0.114 --tau-h 0.2440
# --tau-l 0.1664 --cev-exact; candidate folders (in-repo run, replication package)
COPS_DIRS = (
    os.path.join(PKG, "..", "..", "COPS_replication", "Replication_JPE_Matlab_files",
                 "Direct_preferential_tax_PT", "recalib_v8", "figures_v8"),
    os.path.join(PKG, "cops", "figures"),
)
COPS = {
    "cops_stocks.png": "fig2_response_composition.png",     # fig:cops-stocks
    "cops_taxes.png": "fig4_taxes_budget.png",              # fig:cops-taxes
    "cops_welfare.png": "fig6_welfare_cev.png",             # fig:cops-welfare
    "cops_homotopy.png": "fig7_homotopy.png",               # fig:cops-homotopy
}


def destination() -> str:
    d = os.environ.get("MIG_PAPER_FIGS")
    if d:
        return d
    text_figs = os.path.normpath(os.path.join(PKG, "..", "..", "text", "figs"))
    return text_figs if os.path.isdir(text_figs) else os.path.join(OUT, "paper_figures")


def sources() -> dict:
    """Paper file name -> source path (dynamic-exercise figures only if found)."""
    src = {k: os.path.join(OUT, v) for k, v in STATIC.items()}
    for d in COPS_DIRS:
        if all(os.path.exists(os.path.join(d, f)) for f in COPS.values()):
            src.update({k: os.path.join(d, v) for k, v in COPS.items()})
            break
    return src


def build(dest: str = None) -> dict:
    dest = dest or destination()
    os.makedirs(dest, exist_ok=True)
    print(f"build_paper_figures: writing to {dest}")
    src = sources()
    if not any(k in src for k in COPS):
        print("build_paper_figures: dynamic-exercise figures not found; cops_*.png not written")
    for name, path in src.items():
        shutil.copyfile(path, os.path.join(dest, name))
    print(f"build_paper_figures: {len(src)} files")
    return src


if __name__ == "__main__":
    build()
