# analysis/build_paper_tables.py
"""
Write the tables of the paper draft (text/main_final.tex, text/appendix_final.tex).

Each table environment of the draft reads its tabular block(s) with
\\input{tables/<name>}; this script writes those files. Captions, labels and notes stay
in the draft. A table with two panels has one file per panel.

Sources: only existing outputs, written by run_paper.py:
- output/tables/{csv,tex}/ (main run) and output/variants/{eps20,gamma005,eu_born,
  sigma_theta_grid}/;
- the calibration (config.py, through a subprocess for the epsilon = 20 specification);
- data/eurostat/ (data tables);
- the dynamic exercise (Appendix M): summarize_cops_run.py applied to the saved .mat files
  of the closed-budget runs; skipped with a message when they are absent.

Destination: $MIG_PAPER_TABLES; default ../../text/tables relative to the package root if
that directory exists, else output/paper_tables/.

Usage (from the package root): python -m analysis.build_paper_tables
"""
import glob
import importlib.util
import json
import os
import re
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PKG, "output")
TAUL = 0.218            # progressivity of the paper's runs (config.τl)
PREF_TAG = 11           # int(100 * pref_tax), pref_tax = 0.114
SUBS = ((r"_{10}", r"_{F0}"), (r"_{01}", r"_{F1}"), (r"_{ji}", r"_{Fi}"))


# ---------------------------------------------------------------- helpers
def destination() -> str:
    d = os.environ.get("MIG_PAPER_TABLES")
    if d:
        return d
    text_tables = os.path.normpath(os.path.join(PKG, "..", "..", "text", "tables"))
    return text_tables if os.path.isdir(text_tables) else os.path.join(OUT, "paper_tables")


def _fmt(x) -> str:
    # filename tag of a parameter value (config._fmt)
    return str(int(x)) if float(x).is_integer() else str(x).replace(".", "p")


def _round(x, nd: int) -> str:
    """Half-up rounding of the decimal value as written (1.275 -> 1.28), trailing zeros dropped."""
    q = Decimal(str(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP)
    s = format(q.normalize(), "f")
    return s


def _shape(pred) -> dict:
    with open(os.path.join(PKG, "model", "shape_params.json"), encoding="utf-8") as f:
        params = json.load(f)
    hits = [v for k, v in params.items() if pred(k)]
    assert len(hits) == 1, "shape_params.json: expected one matching specification"
    return hits[0]


class Spec:
    """Output folder and file tags of one specification."""

    def __init__(self, root: str, shape_pred):
        sh = _shape(shape_pred)
        self.root, self.sigma, self.theta = root, sh["sigma"], sh["theta"]
        self.sig = f"sig_{_fmt(self.sigma)}"
        self.dd = f"sig_{_fmt(self.sigma)}_theta_{_fmt(self.theta)}"

    def tex(self, name: str) -> str:
        return os.path.join(self.root, "tables", "tex", name.format(S=self.sig, D=self.dd, t=TAUL))

    def csv(self, name: str) -> pd.Series:
        p = os.path.join(self.root, "tables", "csv", name.format(S=self.sig, D=self.dd, t=TAUL))
        return pd.read_csv(p).iloc[0]


MAIN = Spec(OUT, lambda k: k == "eps=7,gamma=0,rho=3,phi=0.72")
EPS20 = Spec(os.path.join(OUT, "variants", "eps20"), lambda k: k.startswith("eps=20,"))
LAND = Spec(os.path.join(OUT, "variants", "gamma005"), lambda k: "gamma=0.05" in k)
EUB = Spec(os.path.join(OUT, "variants", "eu_born"), lambda k: "targets=eu_born" in k)

CMP = "baseline_vs_{kind}_" + f"{PREF_TAG}" + "_{tag}_taul_{t}_combined{suffix}.tex"


def cmp_path(spec: Spec, kind: str, dd: bool = False, cty1: bool = False) -> str:
    tag = "{D}" if dd else "{S}"
    return spec.tex(CMP.format(kind=kind, tag=tag, t="{t}", suffix="_cty1" if cty1 else ""))


def _read(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


def _tabulars(text: str) -> list:
    return re.findall(r"\\begin\{tabular\}.*?\\end\{tabular\}", text, re.S)


def _rows(tabular: str) -> dict:
    """Data rows (between \\midrule and \\bottomrule): label -> list of cell strings."""
    body = tabular.split(r"\midrule", 1)[1].split(r"\bottomrule", 1)[0]
    out = {}
    for line in body.strip().splitlines():
        cells = [c.strip() for c in line.strip().rstrip("\\").strip().split("&")]
        out[cells[0]] = cells[1:]
    return out


def _sub(s: str) -> str:
    for a, b in SUBS:
        s = s.replace(a, b)
    return s


def _pct(cell: str) -> str:
    """'+0.1\\%' -> '+0.1'."""
    return cell.replace(r"\%", "")


def _cmp_tabulars(path: str) -> list:
    """The two tabulars (allocation; taxes, consumption, CEV) of a comparison table."""
    tabs = _tabulars(_read(path))
    assert len(tabs) == 2, path
    return [_sub(t) for t in tabs]


def comparison(path: str, vspace: str) -> str:
    a, b = _cmp_tabulars(path)
    return f"{a}\n\n\\vspace{{{vspace}}}\n\n{b}\n"


def _line(label: str, cells) -> str:
    return " & ".join([label] + list(cells)) + r" \\"


# ---------------------------------------------------------------- calibration tables
def _calibration(eps: int) -> dict:
    code = (
        "import json, config as c\n"
        "print(json.dumps({'rho': c.ρ, 'eps': c.ϵ, 'phi': c.φ, 'S': c.S, 'kl': c.κl_cal,"
        " 'kh': c.κh_cal, 'sigma': c.σ, 'theta': c.θ, 'tau': c.τl, 'lam': c.CALIBRATION.lam,"
        " 'G': c.G0, 'g_share': c.TARGETS.g_share, 'wage_ratio': c.WAGE_RATIO_TARGET}))\n"
    )
    env = {k: v for k, v in os.environ.items() if not k.startswith("MIG_")}
    env["MIG_EPS"] = str(eps)
    r = subprocess.run([sys.executable, "-c", code], cwd=PKG, env=env, check=True,
                       capture_output=True, text=True)
    return json.loads(r.stdout.strip().splitlines()[-1])


def tab_model_param() -> str:
    c = _calibration(7)
    rows = [
        (r"$\rho$", "Elasticity of substitution between skills", f"{c['rho']:.1f}", "Card (2009)"),
        (r"$\epsilon$", "Nat-immig. elasticity", f"{c['eps']:.1f}", "Manacorda et al. (2012)"),
        (r"$\varphi$", "Immig. productivity", f"{c['phi']:g}",
         f"Imm./native wage ratio {c['wage_ratio']:.1f}"),
        (r"$S$", "High skill productivity", f"{c['S']:.2f}", "College wage premium"),
        (r"$\kappa^l$", "Migration cost: low skill", f"{c['kl']:.2f}", "Immigrant shares"),
        (r"$\kappa^h$", "Migration cost: high skill", f"{c['kh']:.2f}", "Immigrant shares"),
        (r"$\sigma$", r"Fr\'echet shape (stay/leave)", f"{c['sigma']:g}", "Kleven et al. (2013)"),
        (r"$\theta$", r"Fr\'echet shape (destinations)", f"{c['theta']:g}", "Kleven et al. (2014)"),
        (r"$\tau$", "Income tax progressivity", f"{c['tau']:g}", "Holter et al. (2023)"),
        (r"$\lambda$", "Income tax level", f"{c['lam']:.2f}", "Budget clearing"),
        (r"$G$", "Government consumption", f"{c['G']:.2f}",
         f"{100 * c['g_share']:.0f}\\% of Output"),
    ]
    lines = [r"\begin{tabular}{l|l|c|l}", r"\hline",
             r"\textbf{Symbol} & \textbf{Description} & \textbf{Value} & \textbf{Source/Target} \\",
             r"\hline"]
    for r in rows:
        lines += [_line(r[0], r[1:]), r"\hline"]
    return "\n".join(lines + [r"\end{tabular}"]) + "\n"


def tab_model_param_eps20() -> str:
    c = _calibration(20)
    rows = [
        (r"$\rho$", "Elasticity of substitution between skills", f"{c['rho']:.1f}", "Card (2009)"),
        (r"$\epsilon$", "Nat-Immig. elasticity", f"{c['eps']:.1f}", "Card (2009)"),
        (r"$\varphi$", "Immig. productivity", f"{c['phi']:g}",
         f"Imm./native wage ratio {c['wage_ratio']:.1f}"),
        (r"$S$", "High skill productivity", f"{c['S']:.2f}", "College wage premium"),
        (r"$\kappa^l$", "Migration cost: low skill", f"{c['kl']:.2f}", r"Table \ref{tab:population}"),
        (r"$\kappa^h$", "Migration cost: high skill", f"{c['kh']:.2f}", r"Table \ref{tab:education}"),
        (r"$\sigma$", r"Fr\'echet disp. (stay/leave)", f"{c['sigma']:g}", "Kleven et al. (2013)"),
        (r"$\theta$", r"Fr\'echet disp. (destinations)", f"{c['theta']:g}", "Kleven et al. (2014)"),
        (r"$\tau$", "Income tax progressivity", f"{c['tau']:g}", "Holter et al. (2023)"),
        (r"$\lambda$", "Income tax level", f"{c['lam']:.2f}", "Budget clearing"),
        (r"$G$", "Government consumption", f"{c['G']:.2f}",
         f"{100 * c['g_share']:.0f}\\% of Output"),
    ]
    lines = [r"\begin{tabular}{ p{0.5cm}|p{5cm}|p{1.2cm}|p{7cm} }", r"\hline",
             r" & Description & Value & Source/Target \\", r"\hline"]
    for r in rows:
        lines += [_line(r[0], r[1:]), r"\hline"]
    return "\n".join(lines + [r"\end{tabular}"]) + "\n"


def tab_baseline_vars(spec: Spec) -> str:
    return _sub(_tabulars(_read(spec.tex("baseline_model_vars_i0.tex")))[0]) + "\n"


# ---------------------------------------------------------------- data tables
def _data_path(*parts) -> str:
    for base in (os.path.join(PKG, "data"), os.path.join(PKG, "..", "..", "data")):
        p = os.path.join(base, *parts)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(os.path.join("data", *parts))


def tab_education() -> str:
    import config
    d = pd.read_csv(_data_path("eurostat", "edat_lfs_9915_2024_paper_release.csv"),
                    comment="#").set_index("c_birth")
    assert abs(d.loc["NAT", "ED5-8"] - 100 * config.TARGETS_BASE.tertiary_native) < 1e-9
    assert abs(d.loc["FOR", "ED5-8"] - 100 * config.TARGETS_BASE.tertiary_foreign) < 1e-9
    lines = [r"\begin{tabular}{ p{3cm} | p{2cm} | p{2cm} | p{2cm} }", r"\hline",
             r" & Levels 0-2 & Levels 3-4 & Levels 5-8\\", r"\hline"]
    for _, r in d.iterrows():
        lines.append(_line(r["group"], [f"{r[c]:.1f}" for c in ("ED0-2", "ED3_4", "ED5-8")]))
    return "\n".join(lines + [r"\hline", r"\end{tabular}"]) + "\n"


def _jsonstat(path: str, **sel) -> float:
    d = json.load(open(path, encoding="utf-8"))
    idx, stride = 0, 1
    for dim, size in reversed(list(zip(d["id"], d["size"]))):
        cat = d["dimension"][dim]["category"]["index"]
        pos = cat[sel[dim]] if dim in sel else 0
        assert dim in sel or size == 1, f"{dim} not selected"
        idx += pos * stride
        stride *= size
    return d["value"][str(idx)]


def tab_population() -> str:
    import config
    p = _data_path("eurostat", "eu_born", "migr_pop3ctb_EU27_2024.json")
    nat_all = _jsonstat(p, c_birth="NAT", age="TOTAL")
    for_all = _jsonstat(p, c_birth="FOR", age="TOTAL")
    ages = [f"Y{a}-{a + 4}" for a in range(25, 65, 5)]
    nat = sum(_jsonstat(p, c_birth="NAT", age=a) for a in ages)
    fgn = sum(_jsonstat(p, c_birth="FOR", age=a) for a in ages)
    assert nat + fgn == config.POP_25_64 and fgn == config.FOREIGN_25_64
    sp = lambda n: f"{int(n):,}".replace(",", " ")
    tot = nat + fgn
    lines = [r"\begin{tabular}{ p{3cm} | p{2.5cm} | p{2.5cm} | p{2cm} }", r"\hline",
             r" & Total & 25-64 & Share \\", r"\hline",
             _line("Total", [sp(nat_all + for_all), sp(tot), f"{1:.2f}"]),
             _line("Native", [sp(nat_all), sp(nat), f"{nat / tot:.2f}"]),
             _line("Foreign", [sp(for_all), sp(fgn), f"{fgn / tot:.2f}"]),
             r"\hline", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def tab_earnings() -> str:
    import config
    e, n = config.SES22_EARNINGS, config.SES22_EMPLOYEES
    th = lambda x: f"{int(round(x)):,}".replace(",", r"\,")
    lines = [r"\begin{tabular}{ l | c | c | c | c }", r"\hline",
             r" & Levels 0--2 & Levels 3--4 & Non-college & College \\", r"\hline\hline",
             _line(r"Mean annual gross earnings (\texteuro)",
                   [th(e["ED0-2"]), th(e["ED3_4"]), th(config._noncollege_earnings), th(e["ED5-8"])]),
             _line("Employees (millions)",
                   [f"{n['ED0-2'] / 1e6:.1f}", f"{n['ED3_4'] / 1e6:.1f}",
                    f"{(n['ED0-2'] + n['ED3_4']) / 1e6:.1f}", f"{n['ED5-8'] / 1e6:.1f}"]),
             r"\hline", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def tab_nested_elasticities() -> str:
    # formulas only (no numbers)
    return "\n".join([
        r"\begin{tabular}{>{\raggedright\arraybackslash}p{9cm} c c}", r"\toprule",
        r"Response & Elasticity & Parameter \\", r"\midrule",
        r"Natives staying at home, to the home net wage & $\sigma\,(1-P^a_{j,j})$ & $\sigma$ \\",
        r"Natives leaving, to the home net wage & $-\sigma\, P^a_{j,j}$ & $\sigma$ \\",
        r"Natives leaving, to a proportional change in all net wages abroad & $\sigma\, P^a_{j,j}$ & $\sigma$ \\",
        r"Immigration into a small country, to its net wage & $\theta$ & $\theta$ \\",
        r"\quad independent shocks across destinations, $\theta = \sigma$ & $\sigma$ & \\",
        r"\quad perfect substitutes, $\theta \to \infty$ & $\infty$ & \\",
        r"Immigration into a large country (two-country economy), to its net wage & $\sigma\,(1-P^a_{j,i})$ & $\sigma$ \\",
        r"Welfare of workers from other countries, to the net wage of a small country & $0$ & \\",
        r"Welfare of natives, to the home net wage & $P^a_{j,j}$ & $\sigma$ \\",
        r"\bottomrule", r"\end{tabular}"]) + "\n"


# ---------------------------------------------------------------- model tables
def tab_theta_sweep(spec: Spec) -> str:
    return _sub(_tabulars(_read(spec.tex("theta_sweep_pref_tax_" + f"{PREF_TAG}" + "_{S}_taul_{t}.tex")))[0]) + "\n"


def tab_ps_robust() -> str:
    """Perfect-substitutes limit: constant spending (P) and rival share 0.62 (P')."""
    a, b = _cmp_tabulars(cmp_path(MAIN, "soe_pref_tax"))
    ga, gb = _cmp_tabulars(cmp_path(MAIN, "grival_0p62_soe_pref_tax"))
    ra, rb, rga, rgb = _rows(a), _rows(b), _rows(ga), _rows(gb)
    assert ra["B"] == rga["B"] and rb["B"] == rgb["B"]
    head_a = a.split(r"\midrule")[0]
    head_b = b.split(r"\midrule")[0]
    ta = [head_a + r"\midrule", _line("B", ra["B"]), _line("P", ra["P"]),
          _line("P$'$", rga["P"]), r"\bottomrule", r"\end{tabular}"]
    tb = [head_b + r"\midrule", _line("B", rb["B"]), _line("P", rb["P"]),
          _line("{CEV}", rb["{CEV}"]), _line("P$'$", rgb["P"]),
          _line("{CEV$'$}", rgb["{CEV}"]), r"\bottomrule", r"\end{tabular}"]
    return "\n".join(ta) + "\n\n\\vspace{0.3cm}\n\n" + "\n".join(tb) + "\n"


def _cev(path: str, cty_block: int = 1) -> list:
    """CEV row of a comparison table: [c00l, c10l, c00h, c10h] cells without \\%."""
    cells = _rows(_tabulars(_read(path))[cty_block])["{CEV}"]
    return [_pct(c) for c in cells[4:]]


def _welfare_table(rows) -> str:
    lines = [r"\begin{tabular}{l c c c c}", r"\toprule",
             r"\textbf{Scenario} & \multicolumn{4}{c}{\textbf{CEV (\%)}} \\", r"\cmidrule(lr){2-5}",
             r"& Low-skill & Low-skill & High-skill & High-skill \\",
             r"& natives & immigrants & natives & immigrants \\", r"\midrule"]
    lines += [_line(label, cells) for label, cells in rows]
    return "\n".join(lines + [r"\bottomrule", r"\end{tabular}"]) + "\n"


def tab_welfare_summary() -> str:
    ts = _rows(_tabulars(_read(MAIN.tex("theta_sweep_pref_tax_" + f"{PREF_TAG}" + "_{S}_taul_{t}.tex")))[0])
    rows = [
        (f"Small country, $\\theta = {_round(MAIN.theta, 1)}$", _cev(cmp_path(MAIN, "diff_dest_pref_tax", dd=True))),
        (r"Small country, $\theta = 5$", [_pct(c) for c in ts["5"][2:]]),
        (r"Small country, $\theta = 10$", [_pct(c) for c in ts["10"][2:]]),
        ("Small country, perfect-substitutes limit", _cev(cmp_path(MAIN, "soe_pref_tax"))),
        ("Large adopting country", _cev(cmp_path(MAIN, "two_country_pref"))),
        ("Non-adopting country", _cev(cmp_path(MAIN, "two_country_pref", cty1=True))),
        ("Global (multilateral)", _cev(cmp_path(MAIN, "global_pref"))),
    ]
    return _welfare_table(rows)


def tab_welfare_summary_land() -> str:
    rows = [
        (f"Small country, $\\theta = {_round(LAND.theta, 2)}$", _cev(cmp_path(LAND, "diff_dest_pref_tax", dd=True))),
        ("Small country, perfect-substitutes limit", _cev(cmp_path(LAND, "soe_pref_tax"))),
        ("Large country (unilateral)", _cev(cmp_path(LAND, "two_country_pref"))),
        ("Global (multilateral)", _cev(cmp_path(LAND, "global_pref"))),
    ]
    return _welfare_table(rows)


def tab_sigma(name: str) -> str:
    """Sensitivity to σ and θ: σ labels at two decimals, the baseline θ row at one."""
    t = _tabulars(_read(os.path.join(OUT, "variants", "sigma_theta_grid", f"{name}.tex")))[0]
    t = re.sub(r"\\sigma = ([0-9.]+)\$", lambda m: f"\\sigma = {_round(m.group(1), 2)}$", t)
    base = f"{MAIN.theta:g}"
    t = re.sub(rf"(?m)^{re.escape(base)} &", f"{_round(MAIN.theta, 1)} &", t)
    return t + "\n"


def tab_temporary_tax() -> str:
    """Permanent preferential tax and its temporary (5 and 10 of 25 years) equivalents."""
    base = MAIN.csv("baseline_{S}_taul_{t}_equilibrium.csv")
    rates = [0.114] + [s * 0.114 + (1 - s) * base["tau10h"] for s in (5 / 25, 10 / 25)]
    tags = [int(r * 100) for r in rates]
    cases = [("diff_dest_pref_tax_{tag}_{D}", True, "diff_dest_pref_tax"),
             ("soe_pref_tax_{tag}_{S}", False, "soe_pref_tax"),
             ("two_country_pref_{tag}_{S}", False, "two_country_pref")]
    cols = []   # per case, per rate: (rate, pi_h, pi_l, sum, cev l, cev h, cev imm h)
    for eq, dd, kind in cases:
        for tag in tags:
            r = MAIN.csv(eq.replace("{tag}", str(tag)) + "_taul_{t}_equilibrium.csv")
            p = os.path.join(MAIN.root, "tables", "tex", f"baseline_vs_{kind}_{tag}_"
                             f"{MAIN.dd if dd else MAIN.sig}_taul_{TAUL}_combined.tex")
            a, b = (_rows(x) for x in _tabulars(_read(p)))
            cev = b["{CEV}"]
            cols.append((f"{r['tau10h']:.3f}", a["P"][3], a["P"][1], a["P"][4],
                         _pct(cev[4]), _pct(cev[6]), _pct(cev[7])))
    labels = [r"$\tau_{F0}^h$ (pref.\ rate)", r"$\pi_{F0}^h$ (h-skill imm.)",
              r"$\pi_{F0}^l$ (l-skill imm.)", r"$\Sigma\pi$ (population)",
              r"CEV native l-skill (\%)", r"CEV native h-skill (\%)", r"CEV h-skill imm.\ (\%)"]
    lines = [r"\begin{tabular}{l|ccc|ccc|ccc}", r"\hline",
             f" & \\multicolumn{{3}}{{c|}}{{\\textbf{{Small, $\\theta = {_round(MAIN.theta, 1)}$}}}} & "
             r"\multicolumn{3}{c|}{\textbf{Small, $\theta = \infty$}} & \multicolumn{3}{c}{\textbf{Large}} \\",
             r" & Perm. & $T=5$ & $T=10$ & Perm. & $T=5$ & $T=10$ & Perm. & $T=5$ & $T=10$ \\",
             r"\hline"]
    for i, lab in enumerate(labels):
        cells = [c[i] for c in cols]
        if i >= 4:
            cells = [f"${c}$" for c in cells]
        lines.append(_line(lab, cells))
    return "\n".join(lines + [r"\hline", r"\end{tabular}"]) + "\n"


def tab_eu_born() -> str:
    def col(spec: Spec):
        b = spec.csv("baseline_{S}_taul_{t}_equilibrium.csv")
        eqs = [spec.csv("diff_dest_pref_tax_" + f"{PREF_TAG}" + "_{D}_taul_{t}_equilibrium.csv"),
               spec.csv("two_country_pref_" + f"{PREF_TAG}" + "_{S}_taul_{t}_equilibrium.csv"),
               spec.csv("global_pref_" + f"{PREF_TAG}" + "_{S}_taul_{t}_equilibrium.csv"),
               spec.csv("soe_pref_tax_" + f"{PREF_TAG}" + "_{S}_taul_{t}_equilibrium.csv")]
        def dpi(r):
            x = 100 * (r["p10h"] / b["p10h"] - 1)
            return f"{x:.0f}\\%" if abs(x) >= 100 else f"{x:.1f}\\%"
        alloc = [("", f"{1 - b['λ0']:.3f}")] + [(dpi(r), f"{1 - r['λ0']:.3f}") for r in eqs]
        cev = []
        for kind, dd in (("diff_dest_pref_tax", True), ("two_country_pref", False),
                         ("global_pref", False), ("soe_pref_tax", False)):
            c = _cev(cmp_path(spec, kind, dd=dd))
            cev.append(f"${c[0]}$ / ${c[2]}$")
        star = spec.csv("theta_star_pref_tax_" + f"{PREF_TAG}" + "_{S}_taul_{t}.csv")["theta"]
        return alloc, cev, f"{star:.2f}"

    (a0, c0, s0), (a1, c1, s1) = col(MAIN), col(EUB)
    names = ["Small country", "Two countries, unilateral", "Two countries, multilateral",
             "Perfect substitutes"]
    mc = lambda x: f"\\multicolumn{{2}}{{c}}{{{x}}}"
    lines = [r"\begin{tabular}{l cc c cc}", r"\toprule",
             r" & \multicolumn{2}{c}{Baseline} & & \multicolumn{2}{c}{EU-born} \\",
             r"\cmidrule{2-3}\cmidrule{5-6}",
             r" & $\Delta\pi^h_{10}$ & $1-\lambda_0$ & & $\Delta\pi^h_{10}$ & $1-\lambda_0$ \\",
             r"\midrule"]
    for lab, x0, x1 in zip(["No preferential tax"] + names, a0, a1):
        lines.append(_line(lab, [x0[0], x0[1], "", x1[0], x1[1]]))
    lines += [r"\midrule", r" & \multicolumn{5}{c}{CEV of natives, low skill / high skill (\%)} \\",
              r"\midrule"]
    for lab, x0, x1 in zip(names, c0, c1):
        lines.append(_line(lab, [mc(x0), "", mc(x1)]))
    lines += [r"\midrule", _line(r"Self-financing $\theta^*$", [mc(s0), "", mc(s1)]),
              r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- dynamic exercise (Appendix M)
COPS_CANDIDATES = (   # (folder with summarize_cops_run.py, run subfolder)
    (os.path.join(PKG, "..", "..", "COPS_replication", "Replication_JPE_Matlab_files",
                  "Direct_preferential_tax_PT"), "recalib_v8"),
    (os.path.join(PKG, "cops"), "runs"),
)
COPS_PREFIX, COPS_TAUP = "v8", 0.114


def _cops_summary():
    for folder, runs in COPS_CANDIDATES:
        script = os.path.join(folder, "summarize_cops_run.py")
        d = os.path.join(folder, runs)
        if os.path.exists(script) and glob.glob(os.path.join(d, f"cev_paths_{COPS_PREFIX}_*.mat")):
            spec = importlib.util.spec_from_file_location("summarize_cops_run", script)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            from pathlib import Path
            return {cc: mod.summary(Path(d), COPS_PREFIX, cc, COPS_TAUP) for cc in mod.CTRY}
    return None


def _s(x: float, nd: int) -> str:
    v = round(x, nd)
    return f"${v:+.{nd}f}$" if v != 0 else f"${0:.{nd}f}$"


def tab_cops_transition(out: dict) -> str:
    cc = ("pt", "es", "de")
    rows = [
        (r"High-skill stock (\% dev.)", [_s(out[c]["h_stock"], 1) for c in cc]),
        (r"\quad preferential pool (\% of high-skill)", [f"${out[c]['pool']:.1f}$" for c in cc]),
        (r"Low-skill stock (\% dev.)", [_s(out[c]["l_stock"], 1) for c in cc]),
        (r"High-skill wage (\% dev.)", [_s(out[c]["w_h"], 1) for c in cc]),
        (r"Low-skill wage (\% dev.)", [_s(out[c]["w_l"], 1) for c in cc]),
        ("Standard high-skill average rate",
         [f"${out[c]['tau_h0']:.4f} \\to {out[c]['tau_h1']:.4f}$" for c in cc]),
    ]
    lines = [r"\begin{tabular}{ l | c c c }", r"\toprule", r" & Portugal & Spain & Germany \\",
             r"\midrule"] + [_line(a, b) for a, b in rows] + [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def tab_cops_cev(out: dict) -> str:
    ex = lambda x: 100 * np.expm1(x / 100)     # log points -> exact percent
    rows = [("Preferential high-skill, EU-15", "pref_eu15", 1),
            ("Preferential high-skill, NMS-8", "pref_nms8", 1),
            ("Preferential high-skill, rest of world", "pref_oth", 1),
            ("Standard-schedule high-skill", "std_eu15", 1),
            ("Low-skill", "low_eu15", 3)]
    lines = [r"\begin{tabular}{ l | c c | c c | c c }", r"\toprule",
             r" & \multicolumn{2}{c|}{Portugal} & \multicolumn{2}{c|}{Spain} & \multicolumn{2}{c}{Germany} \\",
             r" & announc. & end & announc. & end & announc. & end \\", r"\midrule"]
    for lab, key, nd in rows:
        lines.append(_line(lab, [_s(ex(out[c][key][k]), nd) for c in ("pt", "es", "de") for k in (0, 1)]))
    return "\n".join(lines + [r"\bottomrule", r"\end{tabular}"]) + "\n"


# ---------------------------------------------------------------- registry
def tables() -> dict:
    """File name -> function returning its content. Comments: draft label."""
    c = lambda spec, kind, dd=False, cty1=False, v="1cm": (
        lambda: comparison(cmp_path(spec, kind, dd=dd, cty1=cty1), v))
    return {
        # main text
        "tab_model_param": tab_model_param,                                      # tab:model-param
        "tab_diff_dest": c(MAIN, "diff_dest_pref_tax", dd=True, v="0.3cm"),      # tab:diff_dest_pref_tax_13_...
        "tab_two_country_a": c(MAIN, "two_country_pref", v="0.3cm"),             # tab:two_country_pref_16_..., panel A
        "tab_two_country_b": c(MAIN, "two_country_pref", cty1=True, v="0.3cm"),  # panel B
        "tab_multilateral": c(MAIN, "global_pref", v="0.3cm"),                   # tab:pref_16_sig_2_large_gbc_vs_baseline
        "tab_welfare_summary": tab_welfare_summary,                              # tab:baseline_welfare_comparison
        "sigma_theta_grid": lambda: tab_sigma("sigma_theta_grid"),               # tab:sigma_theta_grid
        "sigma_summary": lambda: tab_sigma("sigma_summary"),                     # tab:sigma_summary
        "tab_ps_robust": tab_ps_robust,                                          # tab:ps-robust
        # appendix
        "tab_nested_elasticities": tab_nested_elasticities,                      # tab:nested-elasticities
        "tab_education": tab_education,                                          # tab:education
        "tab_population": tab_population,                                        # tab:population
        "tab_earnings": tab_earnings,                                            # tab:earnings
        "tab_theta_sweep": lambda: tab_theta_sweep(MAIN),                        # tab:theta_sweep_sig_1p323
        "tab_case1_soe": c(MAIN, "case1_soe_pref_tax"),                          # tab:case1_soe_...
        "tab_diff_dest_taul": c(MAIN, "diff_dest_taul_pref_tax", dd=True),       # tab:diff_dest_taul_theta_2
        "tab_two_country_taul": c(MAIN, "two_country_taul_case2_pref"),          # tab:two_country_taul_case2_...
        "tab_multilateral_taul": c(MAIN, "two_country_taul_case3_pref"),         # tab:global_..._taul
        "tab_gpc_diff_dest": c(MAIN, "diff_dest_gpc_pref_tax", dd=True),         # tab:g-per-capita
        "tab_gpc_ps": c(MAIN, "gpc_soe_pref_tax"),                               # tab:g-per-capita-ps
        "tab_gpc_two_country_a": c(MAIN, "gpc_two_country_pref", v="0.5cm"),     # tab:g-per-capita-lc, panel A
        "tab_gpc_two_country_b": c(MAIN, "gpc_two_country_pref", cty1=True, v="0.5cm"),
        "tab_grival_diff_dest": c(MAIN, "diff_dest_grival_0p62_pref_tax", dd=True),  # tab:g-rival
        "tab_grival_two_country_a": c(MAIN, "grival_0p62_two_country_pref", v="0.5cm"),  # tab:g-rival-lc
        "tab_grival_two_country_b": c(MAIN, "grival_0p62_two_country_pref", cty1=True, v="0.5cm"),
        "tab_model_param_eps20": tab_model_param_eps20,                          # tab:model-param-eps-20
        "tab_baseline_vars_eps20": lambda: tab_baseline_vars(EPS20),             # tab:baseline_..._model_vars_i0
        "tab_two_country_eps20": c(EPS20, "two_country_pref"),                   # tab:two_country_pref_13_..._vs_baseline
        "tab_multilateral_eps20": c(EPS20, "global_pref"),                       # tab:global_..._eps20
        "tab_diff_dest_eps20": c(EPS20, "diff_dest_pref_tax", dd=True),          # tab:diff_dest_eps20
        "tab_theta_sweep_eps20": lambda: tab_theta_sweep(EPS20),                 # tab:theta_sweep_eps20
        "tab_diff_dest_land": c(LAND, "diff_dest_pref_tax", dd=True),            # tab:diff_dest_land
        "tab_ps_land": c(LAND, "soe_pref_tax"),                                  # tab:soe_..._gamma
        "tab_two_country_land": c(LAND, "two_country_pref"),                     # tab:two_country_..._gamma
        "tab_multilateral_land": c(LAND, "global_pref"),                         # tab:global_..._gamma
        "tab_welfare_summary_land": tab_welfare_summary_land,                    # tab:land_welfare_comparison
        "tab_temporary_tax": tab_temporary_tax,                                  # tab:temporary-tax
        "tab_eu_born": tab_eu_born,                                              # tab:eu-born
    }


def build(dest: str = None) -> dict:
    """Write every table file; returns name -> content. COPS tables only if their runs exist."""
    dest = dest or destination()
    os.makedirs(dest, exist_ok=True)
    print(f"build_paper_tables: writing to {dest}")
    written = {}
    for name, fn in tables().items():
        written[name] = fn()
    out = _cops_summary()
    if out is None:
        print("build_paper_tables: dynamic-exercise runs (.mat) not found; "
              "tab_cops_transition and tab_cops_cev not written")
    else:
        written["tab_cops_transition"] = tab_cops_transition(out)
        written["tab_cops_cev"] = tab_cops_cev(out)
    for name, text in written.items():
        # the closing % keeps the end of the file from adding a space after the tabular
        written[name] = text = text.rstrip("\n") + "%\n"
        with open(os.path.join(dest, f"{name}.tex"), "w", encoding="utf-8") as f:
            f.write(text)
    print(f"build_paper_tables: {len(written)} files")
    return written


if __name__ == "__main__":
    os.chdir(PKG)
    build()
