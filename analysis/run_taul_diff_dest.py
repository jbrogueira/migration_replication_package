# analysis/run_taul_diff_dest.py
# Differentiated-destinations small country with the progressivity adjustment:
# λ0 stays at the baseline λ*, τl clears the budget (app:progressivity exercise).
import numpy as np
import pandas as pd

from model.equilibrium import pi_ge_eq, λ_star
from model.equilibrium_diff_dest import (
    equilibrium_row, err_diff_dest_taul, solve_diff_dest_taul,
)
from config import pref_tax, G0, σ, τl, DD_TAG

CSV_DIR = "output/tables/csv"


def main(pref_tax=pref_tax):
    p0, tau_l_star, τ0 = solve_diff_dest_taul(custom_τ10h=pref_tax)
    fun = err_diff_dest_taul(p0, out=1, tau_l_param=tau_l_star, custom_τ10h=pref_tax)
    row = equilibrium_row(p0, λ_star, fun=fun, τ_override=τ0)
    assert abs(row["TR"] - G0) < 1e-8, "budget clearing failed in the taul exercise"
    row["tau_l"] = tau_l_star

    out_csv = f"{CSV_DIR}/diff_dest_taul_pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}_equilibrium.csv"
    pd.DataFrame([row]).to_csv(out_csv, index=False)
    print(f"diff-dest taul adjustment: tau_l {τl} -> {tau_l_star:.4f} "
          f"(lambda0 fixed at {λ_star:.4f})")
    return row, out_csv


if __name__ == "__main__":
    main()
