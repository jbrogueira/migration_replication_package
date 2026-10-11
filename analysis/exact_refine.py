# analysis/exact_refine.py
"""
Exact budget clearing for the grid-swept experiments.

The sweeps select the first grid row whose revenue reaches G0, so the selected row
over-collects by up to one grid step. `refine_crossing` locates the grid bracket
(the selected row and its predecessor in sweep order, on opposite sides of G0) and
solves revenue = G0 in the swept instrument by brentq, re-solving the allocation at
each trial value with the caller's one-point solver.
"""
import numpy as np
import pandas as pd
from scipy import optimize


def refine_crossing(df: pd.DataFrame, xcol: str, trcol: str, G0: float, solve_at,
                    select: str = "ge", xtol: float = 1e-10) -> pd.DataFrame:
    """
    df       sweep (rows in sweep order), with the instrument in `xcol` and revenue in `trcol`
    solve_at x -> one-row DataFrame at instrument value x (same columns as df)
    select   "ge": the grid selection is the first row with revenue >= G0; "gt": > G0
    Returns the exact one-row DataFrame; raises if the selected row has no bracketing
    predecessor (the revenue curve does not cross G0 on the grid there).
    """
    hit = df[trcol] >= G0 if select == "ge" else df[trcol] > G0
    if not hit.any():
        raise RuntimeError(f"revenue never reaches G0 on the {xcol} grid")
    i = int(np.argmax(hit.to_numpy()))
    if i == 0:
        raise RuntimeError(f"revenue already reaches G0 at the first {xcol} grid point: no bracket")
    x_hi, x_lo = df[xcol].iloc[i], df[xcol].iloc[i - 1]

    f = lambda x: float(solve_at(x)[trcol].iloc[0]) - G0
    x_star = optimize.brentq(f, min(x_lo, x_hi), max(x_lo, x_hi), xtol=xtol)
    row = solve_at(x_star)
    assert abs(float(row[trcol].iloc[0]) - G0) < 1e-7, "exact refinement did not clear the budget"
    return row
