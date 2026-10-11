# analysis/build_equilibrium_tables.py
import os
from pathlib import Path
import pandas as pd

from config import pref_tax, σ, τl, tau_k, ell_l, ell_h, SIG_TAG, DD_TAG
from model.parameters import κl, κh
from model.equations import ri

# ---------- helpers ----------
def _pct_str(x: float) -> str:
    if round(x * 100.0, 1) == 0.0:
        x = 0.0   # no signed zero ("-0.0%") for changes below display precision
    return f"{x:+.1%}".replace("%", r"\%")

def _ell_for(skill: str) -> float:
    return ell_l if skill == "l" else ell_h

def _origin_from_code(code: str) -> int:
    """
    Codes in {00l,10l,00h,10h,01l,11l,01h,11h}; first digit = origin country.
    """
    return 0 if code[0] == "0" else 1

def _pvec_for_country(row: pd.Series, country: int):
    """
    Return [p_il, p_jl, p_ih, p_jh] for `country` suitable for ri(country, pvec).
    If some keys are missing (e.g., no country-1 keys in SOE CSV), mirror country-0.
    """
    if country == 0:
        return [
            float(row.get("p00l", 0.0)), float(row.get("p10l", 0.0)),
            float(row.get("p00h", 0.0)), float(row.get("p10h", 0.0)),
        ]
    else:
        return [
            float(row.get("p11l", row.get("p00l", 0.0))),
            float(row.get("p01l", row.get("p10l", 0.0))),
            float(row.get("p11h", row.get("p00h", 0.0))),
            float(row.get("p01h", row.get("p10h", 0.0))),
        ]

def _inc(row: pd.Series, code: str):
    """
    Net income (your 'INC' = consumption) for `code` ∈ {00l,10l,00h,10h,01l,11l,01h,11h}.
    Preferred: use CSV column 'inc{code}' if present (already labor+land net).
    Fallback:  w_{code}*(1 - tau_{code}) + r_{origin} * ell_{skill} * (1 - tau_k).
    """
    col = f"inc{code}"
    if col in row and pd.notna(row[col]):
        return float(row[col])

    wk, tk = f"w{code}", f"tau{code}"
    if wk not in row or tk not in row or pd.isna(row[wk]) or pd.isna(row[tk]):
        return None

    w = float(row[wk])
    tau = float(row[tk])
    skill = code[-1]              # 'l' or 'h'
    origin = _origin_from_code(code)

    # labor net
    c = w * (1.0 - tau)

    # land net (origin-based)
    pvec = _pvec_for_country(row, origin)
    r_origin = ri(origin, pvec)
    c += r_origin * _ell_for(skill) * (1.0 - tau_k)

    return c

def _land_net(row, code: str) -> float:
    """Net land income contained in inc{code}: inc − w·(1 − tau)."""
    return float(row[f"inc{code}"]) - float(row[f"w{code}"]) * (1.0 - float(row[f"tau{code}"]))


def native_abroad_income(row_b, row_p, a: str) -> float:
    """
    Small country: consumption of a country-0 native abroad in the policy equilibrium —
    the foreign net wage (unchanged: country 0 has measure zero) plus country-0 land
    income, which moves with the policy (land income is origin-based).
    """
    return (float(row_b[f"inc10{a}"]) - _land_net(row_b, f"10{a}")
            + _land_net(row_p, f"00{a}"))


def _kappa(code: str) -> float:
    pair, skill = code[:2], code[-1]
    if pair in ("00", "11"):  # stayers’ utility shifter = 1
        return 1.0
    return κl if skill == "l" else κh

def _needs_policy_mirroring(tag: str) -> bool:
    """
    Old SOE-preferential CSVs may lack (01*,11*) columns. If so, mirror for display.
    Two-country files should NOT mirror. Differentiated-destinations files are handled
    only by build_equilibrium_tables_diff_dest (foreign CEV is analytic there).
    """
    t = tag.lower()
    return ("pref_tax_" in t) and ("large" not in t) and ("two_country" not in t) \
        and ("diff_dest" not in t)

# CES numerator for a given column i (i=0 uses 00/01; i=1 uses 10/11)
def _ces_num_for_i(inc_map, a: str, i: int):
    codes = [f"00{a}", f"01{a}"] if i == 0 else [f"10{a}", f"11{a}"]
    v = []
    for code in codes:
        c = inc_map[code]
        v.append((c / _kappa(code)) ** σ)
    return (sum(v)) ** (1.0 / σ)


# ---------- main ----------
def build_equilibrium_tables(
    baseline_eq_path=f"output/tables/csv/baseline_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    out_dir="output/tables/csv",
    tex_dir="output/tables/tex",
    caption=None,
    label=None,
):
    """
    Compare BASELINE (B) vs POLICY (P). Outputs CSVs and the two-tabular LaTeX.

    Refinement logic (for SOE preferential vs baseline):
      - Compute immigrant %ΔINC (10l and 10h) using current CSVs.
      - If either does NOT round to +0.0% at one decimal, densify λ0 grid by +50 points,
        then re-run BASELINE, SOE, and PREFERENTIAL; reload CSVs; repeat.
      - Stop when both round to +0.0% or MAX_REFINES is hit.
    """
    import numpy as np

    # ensure folders
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tex_dir, exist_ok=True)

    stem = Path(pref_eq_path).stem
    scenario_tag = stem.replace("_equilibrium", "")
    if caption is None:
        caption = "{Comparison of migration measures, wages, taxes, and consumption across baseline (B) and policy (P) scenarios.}"
    if label is None:
        label = f"tab:{scenario_tag}_vs_baseline"

    highlight_tau10h = _needs_policy_mirroring(scenario_tag)

    # ---- Read equilibria
    b = pd.read_csv(baseline_eq_path).iloc[0]
    p = pd.read_csv(pref_eq_path).iloc[0]
    if "budget_crossing" in p and not bool(p["budget_crossing"]):
        print(f"DIAGNOSTIC FLAG: {pref_eq_path} is not a budget-clearing equilibrium; "
              f"table not built")
        return None, None, None

    # ==============================================================
    # AUTO-REFINE (applies to the SOE preferential table; re-runs BASELINE+SOE+PREF)
    # ==============================================================
    MAX_REFINES = 5
    if scenario_tag.startswith("soe_pref_tax_"):
        import importlib
        import config as _cfg
        import analysis.run_baseline as _rbase
        import analysis.run_soe as _rsoe
        import analysis.run_pref_tax as _rpref

        for _iter in range(MAX_REFINES):
            # immigrant percent changes (display precision = 1 decimal)
            inc_b_l = _inc(b, "10l"); inc_p_l = _inc(p, "10l")
            inc_b_h = _inc(b, "10h"); inc_p_h = _inc(p, "10h")
            d10l = (inc_p_l / inc_b_l - 1.0) if inc_b_l else 0.0
            d10h = (inc_p_h / inc_b_h - 1.0) if inc_b_h else 0.0

            if round(d10l * 100.0, 1) == 0.0 and round(d10h * 100.0, 1) == 0.0:
                break

            # densify λ0 grid by +50 points (keep start/stop), then re-run baseline + soe + pref
            try:
                n = len(_cfg.λ0_vals)
                start = float(_cfg.λ0_vals[0]); stop = float(_cfg.λ0_vals[-1])
                _cfg.λ0_vals = np.linspace(start, stop, n + 200)

                # Important: don't reload config (would reset defaults). Reload runners so they re-bind λ0_vals.
                importlib.reload(_rbase)
                importlib.reload(_rsoe)
                importlib.reload(_rpref)

                _rbase.main()
                _rsoe.main()
                _rpref.main()

                # reload updated CSV rows
                b = pd.read_csv(baseline_eq_path).iloc[0]
                p = pd.read_csv(pref_eq_path).iloc[0]
            except Exception as e:
                print(f"⚠️ Auto-refine failed (iter {_iter+1}): {e}")
                break

    # ---------- Table 1 (shares + sum + wages)
    cols_t1 = ["p00l","p10l","p00h","p10h","Σπ","w00l","w10l","w00h","w10h"]
    sum_b = float(b["p00l"] + b["p10l"] + b["p00h"] + b["p10h"])
    sum_p = float(p["p00l"] + p["p10l"] + p["p00h"] + p["p10h"])
    t1 = pd.DataFrame(
        [
            [float(b["p00l"]), float(b["p10l"]), float(b["p00h"]), float(b["p10h"]), sum_b,
             float(b["w00l"]), float(b["w10l"]), float(b["w00h"]), float(b["w10h"])],
            [float(p["p00l"]), float(p["p10l"]), float(p["p00h"]), float(p["p10h"]), sum_p,
             float(p["w00l"]), float(p["w10l"]), float(p["w00h"]), float(p["w10h"])],
        ],
        columns=cols_t1, index=["B","P"],
    )

    # ---------- Incomes (INC) maps for CEV and display
    inc_b = {"00l": _inc(b,"00l"), "10l": _inc(b,"10l"), "00h": _inc(b,"00h"), "10h": _inc(b,"10h")}
    inc_b["01l"] = inc_b["10l"]; inc_b["11l"] = inc_b["00l"]
    inc_b["01h"] = inc_b["10h"]; inc_b["11h"] = inc_b["00h"]

    inc_p = {"00l": _inc(p,"00l"), "10l": _inc(p,"10l"), "00h": _inc(p,"00h"), "10h": _inc(p,"10h")}
    if _needs_policy_mirroring(scenario_tag):
        inc_p["01l"] = native_abroad_income(b, p, "l"); inc_p["11l"] = inc_b["00l"]
        inc_p["01h"] = native_abroad_income(b, p, "h"); inc_p["11h"] = inc_b["00h"]
    else:
        def need_inc(code):
            v = _inc(p, code)
            if v is None:
                raise ValueError(f"{scenario_tag}: missing {code} (needed for CEV) in {pref_eq_path}")
            return v
        inc_p["01l"] = need_inc("01l"); inc_p["11l"] = need_inc("11l")
        inc_p["01h"] = need_inc("01h"); inc_p["11h"] = need_inc("11h")

    # ---------- CEV (on INC)
    num_b_l_i0 = _ces_num_for_i(inc_b, "l", i=0)
    num_b_l_i1 = _ces_num_for_i(inc_b, "l", i=1)
    num_p_l_i0 = _ces_num_for_i(inc_p, "l", i=0)
    num_p_l_i1 = _ces_num_for_i(inc_p, "l", i=1)

    num_b_h_i0 = _ces_num_for_i(inc_b, "h", i=0)
    num_b_h_i1 = _ces_num_for_i(inc_b, "h", i=1)
    num_p_h_i0 = _ces_num_for_i(inc_p, "h", i=0)
    num_p_h_i1 = _ces_num_for_i(inc_p, "h", i=1)

    cev_vals_pct = [
        _pct_str(num_p_l_i0/num_b_l_i0 - 1.0),
        _pct_str(num_p_l_i1/num_b_l_i1 - 1.0),
        _pct_str(num_p_h_i0/num_b_h_i0 - 1.0),
        _pct_str(num_p_h_i1/num_b_h_i1 - 1.0),
    ]

    # ---------- Taxes + "consumption" table (INC used for c-columns)
    disp_B = {"c00l": inc_b["00l"], "c10l": inc_b["10l"], "c00h": inc_b["00h"], "c10h": inc_b["10h"]}
    disp_P_pct = {
        "c00l": _pct_str(inc_p["00l"]/inc_b["00l"] - 1.0),
        "c10l": _pct_str(inc_p["10l"]/inc_b["10l"] - 1.0),
        "c00h": _pct_str(inc_p["00h"]/inc_b["00h"] - 1.0),
        "c10h": _pct_str(inc_p["10h"]/inc_b["10h"] - 1.0),
    }

    tax_cols = ["tau00l","tau10l","tau00h","tau10h"]
    t2 = pd.DataFrame(
        [
            [float(b[c]) for c in tax_cols] + [disp_B["c00l"], disp_B["c10l"], disp_B["c00h"], disp_B["c10h"]],
            [float(p[c]) for c in tax_cols] + [disp_P_pct["c00l"], disp_P_pct["c10l"], disp_P_pct["c00h"], disp_P_pct["c10h"]],
        ],
        columns=tax_cols + ["c00l","c10l","c00h","c10h"], index=["B","P"],
    )

    # ---------- LaTeX (two tabulars + CEV row)
    hdr_t1 = r"& $\pi_{00}^l$ & $\pi_{10}^l$ & $\pi_{00}^h$ & $\pi_{10}^h$ & $\Sigma\pi$ & $w_{00}^l$ & $w_{10}^l$ & $w_{00}^h$ & $w_{10}^h$ \\"
    hdr_t2 = r"& $\tau_{00}^l$ & $\tau_{10}^l$ & $\tau_{00}^h$ & $\tau_{10}^h$ & $c_{00}^l$ & $c_{10}^l$ & $c_{00}^h$ & $c_{10}^h$ \\"
    r2 = lambda x: f"{x:.2f}"

    rowB_t1 = "B & " + " & ".join([
        r2(t1.loc["B","p00l"]), r2(t1.loc["B","p10l"]), r2(t1.loc["B","p00h"]), r2(t1.loc["B","p10h"]),
        r2(t1.loc["B","Σπ"]),
        r2(t1.loc["B","w00l"]), r2(t1.loc["B","w10l"]), r2(t1.loc["B","w00h"]), r2(t1.loc["B","w10h"])
    ]) + r" \\"
    rowP_t1 = "P & " + " & ".join([
        r2(t1.loc["P","p00l"]), r2(t1.loc["P","p10l"]), r2(t1.loc["P","p00h"]), r2(t1.loc["P","p10h"]),
        r2(t1.loc["P","Σπ"]),
        r2(t1.loc["P","w00l"]), r2(t1.loc["P","w10l"]), r2(t1.loc["P","w00h"]), r2(t1.loc["P","w10h"])
    ]) + r" \\"

    tau_cells_B = [r2(float(b[c])) for c in ["tau00l","tau10l","tau00h","tau10h"]]
    tau10hP = r2(float(p["tau10h"]))
    if highlight_tau10h:
        tau10hP = r"\textbf{" + tau10hP + "}"

    rowB_t2 = "B & " + " & ".join([
        *tau_cells_B,
        r2(disp_B["c00l"]), r2(disp_B["c10l"]), r2(disp_B["c00h"]), r2(disp_B["c10h"])
    ]) + r" \\"
    rowP_t2 = "P & " + " & ".join([
        r2(float(p["tau00l"])), r2(float(p["tau10l"])), r2(float(p["tau00h"])), tau10hP,
        disp_P_pct["c00l"], disp_P_pct["c10l"], disp_P_pct["c00h"], disp_P_pct["c10h"]
    ]) + r" \\"

    cev_row = r"{CEV} &  &  &  &  & " + " & ".join(cev_vals_pct) + r" \\"

    tex_lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\vspace{1cm}",
        r"\begin{tabular}{ c | ccccccccc }",
        r"\toprule",
        hdr_t1,
        r"\midrule",
        rowB_t1,
        rowP_t1,
        r"\bottomrule",
        r"\end{tabular}",
        " ",
        r"\vspace{1cm}",
        " ",
        r"\begin{tabular}{ c | cccc c cccc }",
        r"\toprule",
        hdr_t2,
        r"\midrule",
        rowB_t2,
        rowP_t2,
        cev_row,
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{{caption}}}",
        fr"\label{{{label}}}",
        r"\end{table}",
    ]

    # save .tex into tex_dir
    tex_path = os.path.join(tex_dir, f"baseline_vs_{scenario_tag}_combined.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("\n".join([ln for ln in tex_lines if ln != ""]))

    return t1, t2, tex_path



def build_equilibrium_tables1(
    baseline_eq_path=f"output/tables/csv/baseline_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    out_dir="output/tables/csv",
    tex_dir="output/tables/tex",
    caption=None,
    label=None,
):
    """
    Country 1 perspective. Baseline CSV has only Country 0 variables (symmetric),
    so we mirror:
        00 -> 11 (natives in 1),    10 -> 01 (immigrants into 1 from 0)
    """
    import os
    from pathlib import Path
    import pandas as pd

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tex_dir, exist_ok=True)

    stem = Path(pref_eq_path).stem
    scenario_tag = stem.replace("_equilibrium", "")
    if caption is None:
        caption = "{Comparison of migration measures, wages, taxes, and consumption across baseline (B) and policy (P) scenarios.}"
    if label is None:
        label = f"tab:{scenario_tag}_vs_baseline_cty1"


    # ---- Read equilibria
    b = pd.read_csv(baseline_eq_path).iloc[0].copy()
    p = pd.read_csv(pref_eq_path).iloc[0].copy()

    # ---- Mirror BASELINE to create Country-1 keys (since baseline file is Country-0 only)
    mirror_map = {
        # shares
        "p11l": "p00l", "p11h": "p00h",
        "p01l": "p10l", "p01h": "p10h",
        # wages
        "w11l": "w00l", "w11h": "w00h",
        "w01l": "w10l", "w01h": "w10h",
        # taxes
        "tau11l": "tau00l", "tau11h": "tau00h",
        "tau01l": "tau10l", "tau01h": "tau10h",
        # after-tax incomes
        "inc11l": "inc00l", "inc11h": "inc00h",
        "inc01l": "inc10l", "inc01h": "inc10h",
    }
    for new_key, old_key in mirror_map.items():
        b[new_key] = b[old_key]

    # ---------- TABLE 1 (Country 1 view): p11*, p01*, w11*, w01*
    cols_t1 = ["p11l","p01l","p11h","p01h","Σπ","w11l","w01l","w11h","w01h"]
    sum_b = float(b["p11l"] + b["p01l"] + b["p11h"] + b["p01h"])
    sum_p = float(p["p11l"] + p["p01l"] + p["p11h"] + p["p01h"])
    t1 = pd.DataFrame(
        [
            [float(b["p11l"]), float(b["p01l"]), float(b["p11h"]), float(b["p01h"]), sum_b,
             float(b["w11l"]), float(b["w01l"]), float(b["w11h"]), float(b["w01h"])],
            [float(p["p11l"]), float(p["p01l"]), float(p["p11h"]), float(p["p01h"]), sum_p,
             float(p["w11l"]), float(p["w01l"]), float(p["w11h"]), float(p["w01h"])],
        ],
        columns=cols_t1, index=["B","P"],
    )
    #t1.to_csv(os.path.join(out_dir, f"baseline_vs_{scenario_tag}_levels_cty1.csv"), index=True)

    # ---------- Incomes for display (Country 1: 11/01)
    disp_B = {
        "inc11l": float(b["inc11l"]), "inc01l": float(b["inc01l"]),
        "inc11h": float(b["inc11h"]), "inc01h": float(b["inc01h"]),
    }
    disp_P_pct = {
        "inc11l": _pct_str(float(p["inc11l"])/disp_B["inc11l"] - 1.0),
        "inc01l": _pct_str(float(p["inc01l"])/disp_B["inc01l"] - 1.0),
        "inc11h": _pct_str(float(p["inc11h"])/disp_B["inc11h"] - 1.0),
        "inc01h": _pct_str(float(p["inc01h"])/disp_B["inc01h"] - 1.0),
    }

    # ---------- Taxes table (Country 1: tau11*, tau01*)
    tax_cols = ["tau11l","tau01l","tau11h","tau01h"]
    t2 = pd.DataFrame(
        [
            [float(b[c]) for c in tax_cols] + [disp_B["inc11l"], disp_B["inc01l"], disp_B["inc11h"], disp_B["inc01h"]],
            [float(p[c]) for c in tax_cols] + [disp_P_pct["inc11l"], disp_P_pct["inc01l"], disp_P_pct["inc11h"], disp_P_pct["inc01h"]],
        ],
        columns=tax_cols + ["inc11l","inc01l","inc11h","inc01h"], index=["B","P"],
    )
    #t2.to_csv(os.path.join(out_dir, f"baseline_vs_{scenario_tag}_taxes_incomes_cty1.csv"), index=True)

    # ---------- CEV (fill full income maps so your CES helper works for i=0 and i=1)
    inc_b = {
        "00l": float(b["inc00l"]), "10l": float(b["inc10l"]),
        "00h": float(b["inc00h"]), "10h": float(b["inc10h"]),
        "01l": float(b["inc01l"]), "11l": float(b["inc11l"]),
        "01h": float(b["inc01h"]), "11h": float(b["inc11h"]),
    }
    inc_p = {
        "00l": float(p["inc00l"]), "10l": float(p["inc10l"]),
        "00h": float(p["inc00h"]), "10h": float(p["inc10h"]),
        "01l": float(p["inc01l"]), "11l": float(p["inc11l"]),
        "01h": float(p["inc01h"]), "11h": float(p["inc11h"]),
    }

    if _needs_policy_mirroring(scenario_tag):
        # keep your original mirroring rule if you rely on it elsewhere
        inc_p["10l"] = inc_b["10l"]; inc_p["11l"] = inc_b["11l"]
        inc_p["10h"] = inc_b["10h"]; inc_p["11h"] = inc_b["11h"]

    num_b_l_i0 = _ces_num_for_i(inc_b, "l", i=0)
    num_b_l_i1 = _ces_num_for_i(inc_b, "l", i=1)
    num_p_l_i0 = _ces_num_for_i(inc_p, "l", i=0)
    num_p_l_i1 = _ces_num_for_i(inc_p, "l", i=1)

    num_b_h_i0 = _ces_num_for_i(inc_b, "h", i=0)
    num_b_h_i1 = _ces_num_for_i(inc_b, "h", i=1)
    num_p_h_i0 = _ces_num_for_i(inc_p, "h", i=0)
    num_p_h_i1 = _ces_num_for_i(inc_p, "h", i=1)

    cev_vals_pct = [
        _pct_str(num_p_l_i0/num_b_l_i0 - 1.0),
        _pct_str(num_p_l_i1/num_b_l_i1 - 1.0),
        _pct_str(num_p_h_i0/num_b_h_i0 - 1.0),
        _pct_str(num_p_h_i1/num_b_h_i1 - 1.0),
    ]
    
    cev_vals_pct_tex = [cev_vals_pct[1], cev_vals_pct[0], cev_vals_pct[3], cev_vals_pct[2]]

    # ---------- LaTeX (Country 1 perspective: show 1 and 0→1 explicitly)
    hdr_t1 = r"& $\pi_{11}^l$ & $\pi_{01}^l$ & $\pi_{11}^h$ & $\pi_{01}^h$ & $\Sigma\pi$ & $w_{11}^l$ & $w_{01}^l$ & $w_{11}^h$ & $w_{01}^h$ \\"
    hdr_t2 = r"& $\tau_{11}^l$ & $\tau_{01}^l$ & $\tau_{11}^h$ & $\tau_{01}^h$ & $c_{11}^l$ & $c_{01}^l$ & $c_{11}^h$ & $c_{01}^h$ \\"
    r2 = lambda x: f"{x:.2f}"

    rowB_t1 = "B & " + " & ".join([
        r2(t1.loc["B","p11l"]), r2(t1.loc["B","p01l"]), r2(t1.loc["B","p11h"]), r2(t1.loc["B","p01h"]),
        r2(t1.loc["B","Σπ"]),
        r2(t1.loc["B","w11l"]), r2(t1.loc["B","w01l"]), r2(t1.loc["B","w11h"]), r2(t1.loc["B","w01h"])
    ]) + r" \\"
    rowP_t1 = "P & " + " & ".join([
        r2(t1.loc["P","p11l"]), r2(t1.loc["P","p01l"]), r2(t1.loc["P","p11h"]), r2(t1.loc["P","p01h"]),
        r2(t1.loc["P","Σπ"]),
        r2(t1.loc["P","w11l"]), r2(t1.loc["P","w01l"]), r2(t1.loc["P","w11h"]), r2(t1.loc["P","w01h"])
    ]) + r" \\"

    rowB_t2 = "B & " + " & ".join([
        f"{float(b['tau11l']):.2f}", f"{float(b['tau01l']):.2f}",
        f"{float(b['tau11h']):.2f}", f"{float(b['tau01h']):.2f}",
        r2(disp_B["inc11l"]), r2(disp_B["inc01l"]), r2(disp_B["inc11h"]), r2(disp_B["inc01h"])
    ]) + r" \\"
    rowP_t2 = "P & " + " & ".join([
        f"{float(p['tau11l']):.2f}", f"{float(p['tau01l']):.2f}",
        f"{float(p['tau11h']):.2f}", f"{float(p['tau01h']):.2f}",
        disp_P_pct["inc11l"], disp_P_pct["inc01l"], disp_P_pct["inc11h"], disp_P_pct["inc01h"]
    ]) + r" \\"

    cev_row = r"{CEV} &  &  &  &  & " + " & ".join(cev_vals_pct_tex) + r" \\"

    tex_lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\vspace{1cm}",
        r"\begin{tabular}{ c | ccccccccc }",
        r"\toprule",
        hdr_t1,
        r"\midrule",
        rowB_t1,
        rowP_t1,
        r"\bottomrule",
        r"\end{tabular}",
        " ",
        r"\vspace{1cm}",
        " ",
        r"\begin{tabular}{ c | cccc c cccc }",
        r"\toprule",
        hdr_t2,
        r"\midrule",
        rowB_t2,
        rowP_t2,
        cev_row,
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{{caption}}}",
        fr"\label{{{label}}}",
        r"\end{table}",
    ]

    tex_path = os.path.join(tex_dir, f"baseline_vs_{scenario_tag}_combined_cty1.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("\n".join(tex_lines))

    return t1, t2, tex_path


def build_equilibrium_tables_diff_dest(
    baseline_eq_path=f"output/tables/csv/diff_dest_baseline_{DD_TAG}_taul_{τl}_equilibrium.csv",
    pref_eq_path=f"output/tables/csv/diff_dest_pref_tax_{int(pref_tax*100)}_{DD_TAG}_taul_{τl}_equilibrium.csv",
    out_dir="output/tables/csv",
    tex_dir="output/tables/tex",
    caption=None,
    label=None,
):
    """
    BASELINE (B) vs POLICY (P) for the differentiated-destinations small country.

    Both CSV rows are exact budget-clearing solves (run_diff_dest), so there is no
    auto-refine. Native CEV uses the CES numerator with the move option frozen at
    baseline (exact under this structure: the foreign continuum is unaffected).
    Foreign CEV is analytic +0.0%: country 0 has measure zero in foreign choice
    denominators, and welfare conditional on choosing country 0 equals unconditional
    welfare (selection).
    """
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tex_dir, exist_ok=True)

    stem = Path(pref_eq_path).stem
    scenario_tag = stem.replace("_equilibrium", "")
    if caption is None:
        caption = "{Comparison of migration measures, wages, taxes, and consumption across baseline (B) and policy (P) scenarios.}"
    if label is None:
        label = f"tab:{scenario_tag}_vs_baseline"

    # ---- Read equilibria (both exact)
    b = pd.read_csv(baseline_eq_path).iloc[0]
    p = pd.read_csv(pref_eq_path).iloc[0]

    # ---------- Table 1 (shares + sum + wages)
    cols_t1 = ["p00l","p10l","p00h","p10h","Σπ","w00l","w10l","w00h","w10h"]
    sum_b = float(b["p00l"] + b["p10l"] + b["p00h"] + b["p10h"])
    sum_p = float(p["p00l"] + p["p10l"] + p["p00h"] + p["p10h"])
    t1 = pd.DataFrame(
        [
            [float(b["p00l"]), float(b["p10l"]), float(b["p00h"]), float(b["p10h"]), sum_b,
             float(b["w00l"]), float(b["w10l"]), float(b["w00h"]), float(b["w10h"])],
            [float(p["p00l"]), float(p["p10l"]), float(p["p00h"]), float(p["p10h"]), sum_p,
             float(p["w00l"]), float(p["w10l"]), float(p["w00h"]), float(p["w10h"])],
        ],
        columns=cols_t1, index=["B","P"],
    )

    # ---------- Incomes: move option frozen at baseline (mirroring rule)
    inc_b = {"00l": _inc(b,"00l"), "10l": _inc(b,"10l"), "00h": _inc(b,"00h"), "10h": _inc(b,"10h")}
    inc_b["01l"] = inc_b["10l"]; inc_b["01h"] = inc_b["10h"]

    inc_p = {"00l": _inc(p,"00l"), "10l": _inc(p,"10l"), "00h": _inc(p,"00h"), "10h": _inc(p,"10h")}
    inc_p["01l"] = native_abroad_income(b, p, "l"); inc_p["01h"] = native_abroad_income(b, p, "h")

    # ---------- CEV: natives via the CES numerator; foreign analytic +0.0%
    cev_vals_pct = [
        _pct_str(_ces_num_for_i(inc_p, "l", i=0) / _ces_num_for_i(inc_b, "l", i=0) - 1.0),
        _pct_str(0.0),
        _pct_str(_ces_num_for_i(inc_p, "h", i=0) / _ces_num_for_i(inc_b, "h", i=0) - 1.0),
        _pct_str(0.0),
    ]

    # ---------- Taxes + consumption table (INC used for c-columns)
    disp_B = {"c00l": inc_b["00l"], "c10l": inc_b["10l"], "c00h": inc_b["00h"], "c10h": inc_b["10h"]}
    disp_P_pct = {
        "c00l": _pct_str(inc_p["00l"]/inc_b["00l"] - 1.0),
        "c10l": _pct_str(inc_p["10l"]/inc_b["10l"] - 1.0),
        "c00h": _pct_str(inc_p["00h"]/inc_b["00h"] - 1.0),
        "c10h": _pct_str(inc_p["10h"]/inc_b["10h"] - 1.0),
    }

    tax_cols = ["tau00l","tau10l","tau00h","tau10h"]
    t2 = pd.DataFrame(
        [
            [float(b[c]) for c in tax_cols] + [disp_B["c00l"], disp_B["c10l"], disp_B["c00h"], disp_B["c10h"]],
            [float(p[c]) for c in tax_cols] + [disp_P_pct["c00l"], disp_P_pct["c10l"], disp_P_pct["c00h"], disp_P_pct["c10h"]],
        ],
        columns=tax_cols + ["c00l","c10l","c00h","c10h"], index=["B","P"],
    )

    # ---------- LaTeX (two tabulars + CEV row)
    hdr_t1 = r"& $\pi_{00}^l$ & $\pi_{10}^l$ & $\pi_{00}^h$ & $\pi_{10}^h$ & $\Sigma\pi$ & $w_{00}^l$ & $w_{10}^l$ & $w_{00}^h$ & $w_{10}^h$ \\"
    hdr_t2 = r"& $\tau_{00}^l$ & $\tau_{10}^l$ & $\tau_{00}^h$ & $\tau_{10}^h$ & $c_{00}^l$ & $c_{10}^l$ & $c_{00}^h$ & $c_{10}^h$ \\"
    r2 = lambda x: f"{x:.2f}"

    rowB_t1 = "B & " + " & ".join([
        r2(t1.loc["B","p00l"]), r2(t1.loc["B","p10l"]), r2(t1.loc["B","p00h"]), r2(t1.loc["B","p10h"]),
        r2(t1.loc["B","Σπ"]),
        r2(t1.loc["B","w00l"]), r2(t1.loc["B","w10l"]), r2(t1.loc["B","w00h"]), r2(t1.loc["B","w10h"])
    ]) + r" \\"
    rowP_t1 = "P & " + " & ".join([
        r2(t1.loc["P","p00l"]), r2(t1.loc["P","p10l"]), r2(t1.loc["P","p00h"]), r2(t1.loc["P","p10h"]),
        r2(t1.loc["P","Σπ"]),
        r2(t1.loc["P","w00l"]), r2(t1.loc["P","w10l"]), r2(t1.loc["P","w00h"]), r2(t1.loc["P","w10h"])
    ]) + r" \\"

    tau_cells_B = [r2(float(b[c])) for c in tax_cols]
    tau10hP = r"\textbf{" + r2(float(p["tau10h"])) + "}"

    rowB_t2 = "B & " + " & ".join([
        *tau_cells_B,
        r2(disp_B["c00l"]), r2(disp_B["c10l"]), r2(disp_B["c00h"]), r2(disp_B["c10h"])
    ]) + r" \\"
    rowP_t2 = "P & " + " & ".join([
        r2(float(p["tau00l"])), r2(float(p["tau10l"])), r2(float(p["tau00h"])), tau10hP,
        disp_P_pct["c00l"], disp_P_pct["c10l"], disp_P_pct["c00h"], disp_P_pct["c10h"]
    ]) + r" \\"

    cev_row = r"{CEV} &  &  &  &  & " + " & ".join(cev_vals_pct) + r" \\"

    tex_lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\vspace{1cm}",
        r"\begin{tabular}{ c | ccccccccc }",
        r"\toprule",
        hdr_t1,
        r"\midrule",
        rowB_t1,
        rowP_t1,
        r"\bottomrule",
        r"\end{tabular}",
        " ",
        r"\vspace{1cm}",
        " ",
        r"\begin{tabular}{ c | cccc c cccc }",
        r"\toprule",
        hdr_t2,
        r"\midrule",
        rowB_t2,
        rowP_t2,
        cev_row,
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{{caption}}}",
        fr"\label{{{label}}}",
        r"\end{table}",
    ]

    tex_path = os.path.join(tex_dir, f"baseline_vs_{scenario_tag}_combined.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("\n".join([ln for ln in tex_lines if ln != ""]))

    return t1, t2, tex_path


def build_baseline_model_vars_table(
    baseline_eq_path=f"output/tables/csv/baseline_{SIG_TAG}_taul_{τl}_equilibrium.csv",
    tex_dir="output/tables/tex",
    country: int = 0,
    caption=r"{Model variables in the baseline equilibrium.}",
    label=None,
):
    """
    Builds the LaTeX table of baseline model variables for country `i` (0 or 1).
    Uses baseline CSV shares and computes wages, N_i^l, N_i^h, and F_i from equations.py.
    Writes the .tex file and returns its path.
    """
    import os
    import pandas as pd
    from pathlib import Path
    from model.equations import wage, Fi as _Fi, Nia as _Nia

    os.makedirs(tex_dir, exist_ok=True)
    row = pd.read_csv(baseline_eq_path).iloc[0]

    # ----- Build composition vector for country i: pi = [p_iil, p_jil, p_iih, p_jih]
    if country == 0:
        p_iil = float(row["p00l"])
        p_jil = float(row.get("p10l", 0.0))
        p_iih = float(row["p00h"])
        p_jih = float(row.get("p10h", 0.0))
        tau_ll_key, tau_lm_key, tau_hl_key, tau_hm_key = "tau00l", "tau10l", "tau00h", "tau10h"
    else:
        # If baseline has only country-0 keys, mirror them
        p_iil = float(row.get("p11l", row["p00l"]))
        p_jil = float(row.get("p01l", row.get("p10l", 0.0)))
        p_iih = float(row.get("p11h", row["p00h"]))
        p_jih = float(row.get("p01h", row.get("p10h", 0.0)))
        tau_ll_key, tau_lm_key, tau_hl_key, tau_hm_key = "tau11l", "tau01l", "tau11h", "tau01h"

    pi = [p_iil, p_jil, p_iih, p_jih]
    pop_sum = p_iil + p_jil + p_iih + p_jih

    # ----- Wages: wage(j, i, a, pi)
    i = int(country)
    j = 1 - i
    w_native_l = float(wage(i, i, "l", pi))
    w_immig_l  = float(wage(j, i, "l", pi))
    w_native_h = float(wage(i, i, "h", pi))
    w_immig_h  = float(wage(j, i, "h", pi))

    # ----- Average tax rates from CSV if present
    tau_native_l = float(row[tau_ll_key]) if tau_ll_key in row else float("nan")
    tau_immig_l  = float(row[tau_lm_key]) if tau_lm_key in row else float("nan")
    tau_native_h = float(row[tau_hl_key]) if tau_hl_key in row else float("nan")
    tau_immig_h  = float(row[tau_hm_key]) if tau_hm_key in row else float("nan")

    # ----- Aggregates and output
    N_il = float(_Nia(p_iil, p_jil))
    N_ih = float(_Nia(p_iih, p_jih))
    F_i  = float(_Fi(i, pi))

    # ----- Display symbols (keep i/j notation)
    sym_ll  = r"$\pi_{ii}^l$"; sym_lm  = r"$\pi_{ji}^l$"
    sym_hl  = r"$\pi_{ii}^h$"; sym_hm  = r"$\pi_{ji}^h$"
    sym_wll = r"$w_{ii}^l$";   sym_wlm = r"$w_{ji}^l$"
    sym_whl = r"$w_{ii}^h$";   sym_whm = r"$w_{ji}^h$"
    sym_tll = r"$\tau_{ii}^l$"; sym_tlm = r"$\tau_{ji}^l$"
    sym_thl = r"$\tau_{ii}^h$"; sym_thm = r"$\tau_{ji}^h$"
    Nil_label, Nih_label, Fi_label = r"$N_i^l$", r"$N_i^h$", r"$F_i$"

    r2 = lambda x: f"{x:.2f}"
    rows = [
        (sym_ll,   "Share of Native Low-skill",        r2(p_iil)),
        (sym_lm,   "Share of Immigrant Low-skill",     r2(p_jil)),
        (sym_hl,   "Share of Native High-skill",       r2(p_iih)),
        (sym_hm,   "Share of Immigrant High-skill",    r2(p_jih)),
        (r"$\sum\pi$", "Population",                    r2(pop_sum)),
        (sym_wll,  "Wage of Native Low-skill",         r2(w_native_l)),
        (sym_wlm,  "Wage of Imm. Low-skill",           r2(w_immig_l)),
        (sym_whl,  "Wage of Native High-skill",        r2(w_native_h)),
        (sym_whm,  "Wage of Imm. High-skill",          r2(w_immig_h)),
        (sym_tll,  "Avg. tax rate of Native Low-skill",r2(tau_native_l)),
        (sym_tlm,  "Avg. tax rate of Imm. Low-skill",  r2(tau_immig_l)),
        (sym_thl,  "Avg. tax rate of Native High-skill",r2(tau_native_h)),
        (sym_thm,  "Avg. tax rate of Imm. High-skill", r2(tau_immig_h)),
        (Nil_label,"Low-skill Labour",                  r2(N_il)),
        (Nih_label,"High-skill Labour",                 r2(N_ih)),
        (Fi_label, "Output",                            r2(F_i)),
    ]

    # ----- LaTeX
    if label is None:
        scenario_tag = Path(baseline_eq_path).stem.replace("_equilibrium", "")
        label = f"tab:{scenario_tag}_baseline_model_vars_i{i}"

    tex_lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\begin{tabular}{lll}",
        r"\toprule",
        r"\textbf{Variable} & \textbf{Description} & \textbf{Value} \\",
        r"\midrule",
        *[f"{v} & {desc} & {val} \\\\" for (v, desc, val) in rows],
        r"\bottomrule",
        r"\end{tabular}",
        fr"\caption{{{caption}}}",
        fr"\label{{{label}}}",
        r"\end{table}",
    ]

    tex_path = os.path.join(tex_dir, f"baseline_model_vars_i{i}.tex")
    with open(tex_path, "w", encoding="utf-8") as f:
        f.write("\n".join(tex_lines))

    return tex_path
