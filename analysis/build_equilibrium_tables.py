# analysis/build_equilibrium_tables.py
import os
from pathlib import Path
import pandas as pd

from config import pref_tax, σ, τl, tau_k, ell_l, ell_h
from model.parameters import κl, κh
from model.equations import ri

# ---------- helpers ----------
def _pct_str(x: float) -> str:
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
    Net income for group `code` ∈ {00l,10l,00h,10h,01l,11l,01h,11h}.

    Preferred: use CSV column 'inc{code}' if present (already labor+land net).
    Fallback:  w_{code}*(1 - tau_{code}) + r_{origin} * ell_{skill} * (1 - tau_k).
    """
    # 1) Use precomputed income if available
    col = f"inc{code}"
    if col in row and pd.notna(row[col]):
        return float(row[col])

    # 2) Compute labor net + origin-based land net
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

def _kappa(code: str) -> float:
    pair, skill = code[:2], code[-1]
    if pair in ("00", "11"):  # stayers’ utility shifter = 1
        return 1.0
    return κl if skill == "l" else κh

def _needs_policy_mirroring(tag: str) -> bool:
    """
    Old SOE-preferential CSVs may lack (01*,11*) columns. If so, mirror for display.
    Two-country files should NOT mirror.
    """
    t = tag.lower()
    return ("pref_tax_" in t) and ("large" not in t) and ("two_country" not in t)

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
    baseline_eq_path=f"output/tables/csv/baseline_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    out_dir="output/tables/csv",
    tex_dir="output/tables/tex",
    caption=None,
    label=None,
):
    """
    Compare BASELINE (B) vs POLICY (P). Outputs CSVs and the two-tabular LaTeX.
    """
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
    #t1.to_csv(os.path.join(out_dir, f"baseline_vs_{scenario_tag}_levels.csv"), index=True)

    # ---------- After-tax incomes (build full map with mirroring rules)
    inc_b = {
        "00l": _inc(b,"00l"), "10l": _inc(b,"10l"),
        "00h": _inc(b,"00h"), "10h": _inc(b,"10h"),
    }
    inc_b["01l"] = inc_b["10l"]; inc_b["11l"] = inc_b["00l"]
    inc_b["01h"] = inc_b["10h"]; inc_b["11h"] = inc_b["00h"]

    inc_p = {
        "00l": _inc(p,"00l"), "10l": _inc(p,"10l"),
        "00h": _inc(p,"00h"), "10h": _inc(p,"10h"),
    }

    if _needs_policy_mirroring(scenario_tag):
        inc_p["01l"] = inc_b["10l"]; inc_p["11l"] = inc_b["00l"]
        inc_p["01h"] = inc_b["10h"]; inc_p["11h"] = inc_b["00h"]
        pi_p = {
            "00l": float(p["p00l"]), "10l": float(p["p10l"]),
            "00h": float(p["p00h"]), "10h": float(p["p10h"]),
        }
    else:
        def need_inc(code):
            v = _inc(p, code)
            if v is None:
                raise ValueError(f"{scenario_tag}: missing {code} (needed for CEV) in {pref_eq_path}")
            return v
        inc_p["01l"] = need_inc("01l"); inc_p["11l"] = need_inc("11l")
        inc_p["01h"] = need_inc("01h"); inc_p["11h"] = need_inc("11h")
        pi_p = {
            "00l": float(p["p00l"]), "10l": float(p["p10l"]),
            "00h": float(p["p00h"]), "10h": float(p["p10h"]),
        }

    pi_b = {
        "00l": float(b["p00l"]), "10l": float(b["p10l"]),
        "00h": float(b["p00h"]), "10h": float(b["p10h"]),
    }

    # ---------- Display rows
    disp_B = {"inc00l": inc_b["00l"], "inc10l": inc_b["10l"], "inc00h": inc_b["00h"], "inc10h": inc_b["10h"]}
    disp_P_pct = {
        "inc00l": _pct_str(inc_p["00l"]/inc_b["00l"] - 1.0),
        "inc10l": _pct_str(inc_p["10l"]/inc_b["10l"] - 1.0),
        "inc00h": _pct_str(inc_p["00h"]/inc_b["00h"] - 1.0),
        "inc10h": _pct_str(inc_p["10h"]/inc_b["10h"] - 1.0),
    }

    # ---------- Taxes table CSV
    tax_cols = ["tau00l","tau10l","tau00h","tau10h"]
    t2 = pd.DataFrame(
        [
            [float(b[c]) for c in tax_cols] + [disp_B["inc00l"], disp_B["inc10l"], disp_B["inc00h"], disp_B["inc10h"]],
            [float(p[c]) for c in tax_cols] + [disp_P_pct["inc00l"], disp_P_pct["inc10l"], disp_P_pct["inc00h"], disp_P_pct["inc10h"]],
        ],
        columns=tax_cols + ["inc00l","inc10l","inc00h","inc10h"], index=["B","P"],
    )
    #t2.to_csv(os.path.join(out_dir, f"baseline_vs_{scenario_tag}_taxes_incomes.csv"), index=True)

    # ---------- CEV_{i,0} with CES (v=c/κ)
    num_b_l_i0 = _ces_num_for_i(inc_b, "l", i=0)
    num_b_l_i1 = _ces_num_for_i(inc_b, "l", i=1)
    num_p_l_i0 = _ces_num_for_i(inc_p, "l", i=0)
    num_p_l_i1 = _ces_num_for_i(inc_p, "l", i=1)

    num_b_h_i0 = _ces_num_for_i(inc_b, "h", i=0)
    num_b_h_i1 = _ces_num_for_i(inc_b, "h", i=1)
    num_p_h_i0 = _ces_num_for_i(inc_p, "h", i=0)
    num_p_h_i1 = _ces_num_for_i(inc_p, "h", i=1)

    ωB_c0l = num_b_l_i0
    ωP_c0l = num_p_l_i0
    ωB_c1l = num_b_l_i1
    ωP_c1l = num_p_l_i1

    ωB_c0h = num_b_h_i0
    ωP_c0h = num_p_h_i0
    ωB_c1h = num_b_h_i1
    ωP_c1h = num_p_h_i1

    cev_vals_pct = [
        _pct_str(ωP_c0l/ωB_c0l - 1.0),
        _pct_str(ωP_c1l/ωB_c1l - 1.0),
        _pct_str(ωP_c0h/ωB_c0h - 1.0),
        _pct_str(ωP_c1h/ωB_c1h - 1.0),
    ]

    # ---------- LaTeX (two tabulars + CEV row)
    hdr_t1 = r"& $\pi_{00}^l$ & $\pi_{10}^l$ & $\pi_{00}^h$ & $\pi_{10}^h$ & $\Sigma\pi$ & $w_{00}^l$ & $w_{10}^l$ & $w_{00}^h$ & $w_{10}^h$ \\"
    hdr_t2 = r"& $\tau_{00}^l$ & $\tau_{10}^l$ & $\tau_{00}^h$ & $\tau_{10}^h$ & $c_{00}^l$ & $c_{10}^l$ & $c_{0}^h$ & $c_{10}^h$ \\"
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
        r2(disp_B["inc00l"]), r2(disp_B["inc10l"]), r2(disp_B["inc00h"]), r2(disp_B["inc10h"])
    ]) + r" \\"
    rowP_t2 = "P & " + " & ".join([
        r2(float(p["tau00l"])), r2(float(p["tau10l"])), r2(float(p["tau00h"])), tau10hP,
        disp_P_pct["inc00l"], disp_P_pct["inc10l"], disp_P_pct["inc00h"], disp_P_pct["inc10h"]
    ]) + r" \\"

    cev_row = r"{CEV} &  &  &  &  & " + " & ".join(cev_vals_pct) + r" \\"

    tex_lines = [
r"\begin{table}[ht]",
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
    baseline_eq_path=f"output/tables/csv/baseline_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    pref_eq_path=f"output/tables/csv/soe_pref_tax_{int(pref_tax*100)}_sig_{int(σ)}_taul_{τl}_equilibrium.csv",
    out_dir="output/tables/csv",
    tex_dir="output/tables/tex",
    caption=None,
    label=None,
):
    """
    Country 1 perspective. Baseline CSV has only Country 0 variables (symmetric),
    so we mirror:
        00 -> 11 (natives in 1),    10 -> 01 (immigrants into 1 from 0)
    Policy CSV already contains 11/01, so we read as-is.
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

    highlight_tau10h = _needs_policy_mirroring(scenario_tag)  # kept if you use it elsewhere

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
    # Build complete maps (both countries) so _ces_num_for_i can pick i=0/1.
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

    cev_row = r"{CEV} &  &  &  &  & " + " & ".join(cev_vals_pct) + r" \\"

    tex_lines = [
        r"\begin{table}[ht]",
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