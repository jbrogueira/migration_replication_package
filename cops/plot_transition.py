"""Plot the Approach-2 preferential-tax transition (tau_P = 0.13, closed budget).

Reads the saved run outputs in this directory:
  pt_stocks_r4_0p13.mat / pt_stocks_r4_0p19.mat / pt_stocks_r4_0p16.mat /
  pt_stocks_treat2_closed.mat   (dump_pt_stocks: PT stocks by skill/sub-pool)
  r4_taup_0p13.mat              (solution; info.dg per-period paths: taxes, lambda,
                                 wedges, PT wage deviations)
  cev_paths_r4_taup13.mat       (per-period CEV from value functions alone)

Writes PNG figures to figures_transition/. Options select another run: --dir (folder of
the run outputs), --prefix (e.g. v6_pt reads v6_pt_taup_0p13.mat, pt_stocks_v6_pt_*.mat,
cev_paths_v6_pt_taup13.mat), --tau-h / --tau-l (baseline average rates drawn in the tax
figure), --out (output folder). Defaults reproduce the July r4 run. Reform is announced and implemented,
unanticipated, at period 2; the terminal condition pins the value hats at T, so
paths flatten mechanically near the horizon end.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

HERE = Path(__file__).parent
OUT = HERE / "figures_transition"

# Okabe-Ito subset, validated (dataviz six checks). Fixed assignment by entity:
BLOC = {"eu15": "#0072B2", "nms8": "#D55E00", "oth": "#009E73"}  # origin groups
BLOC_NAME = {"eu15": "EU-15", "nms8": "NMS-8", "oth": "RoW"}  # display names
HIGH, LOW = "#0072B2", "#D55E00"  # high-/low-skill when groups are pooled
ACCENT = "#E69F00"  # low contrast on white: use only with a direct label
BASE_GRAY = "#666666"
SEQ_BLUES = ["#9ecae1", "#6baed6", "#3182bd", "#08519c"]  # shallow -> deep cut

TP = 0.13
RUN = {"dir": HERE, "prefix": None, "tau_h": 0.24, "tau_l": 0.18}
T_SHOW = 199  # slot 200 is the pinned terminal slice


def load(name: str):
    return scipy.io.loadmat(RUN["dir"] / name, squeeze_me=True, struct_as_record=False)


@dataclass
class Data:
    st: object            # PT stocks at tau_P = 0.13
    dg: object            # per-period paths at tau_P = 0.13
    cev: object           # CEV paths
    st_homotopy: dict     # tau_P -> stocks


def load_all() -> Data:
    p = RUN["prefix"]
    if p is None:  # July r4 run
        st = load("pt_stocks_r4_0p13.mat")["St"]
        dg = load("r4_taup_0p13.mat")["info"].dg
        cev = load("cev_paths_r4_taup13.mat")["out"]
        st_h = {0.2127: load("pt_stocks_treat2_closed.mat")["St"],
                0.19: load("pt_stocks_r4_0p19.mat")["St"],
                0.16: load("pt_stocks_r4_0p16.mat")["St"], 0.13: st}
    else:
        st = load(f"pt_stocks_{p}_0p13.mat")["St"]
        dg = load(f"{p}_taup_0p13.mat")["info"].dg
        cev = load(f"cev_paths_{p}_taup13.mat")["out"]
        st_h = {0.2127: load(f"pt_stocks_{p}_0p2127.mat")["St"],
                0.19: load(f"pt_stocks_{p}_0p19.mat")["St"],
                0.16: load(f"pt_stocks_{p}_0p16.mat")["St"], 0.13: st}
    return Data(st, dg, cev, st_h)


def style_axis(ax, ylabel: str = "") -> None:
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axvline(2, color=BASE_GRAY, linewidth=0.8, linestyle=":", alpha=0.8)
    ax.set_xlabel("year")
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_xlim(1, T_SHOW)


def save(fig, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=150)
    plt.close(fig)
    print(f"wrote {OUT / name}")


def t_of(x) -> np.ndarray:
    return np.arange(1, len(x) + 1)


def trim(x) -> np.ndarray:
    """Drop the pinned terminal slice (period 200) — its stock entry is degenerate."""
    return np.asarray(x)[:T_SHOW]


def fig1_stock_levels(d: Data) -> None:
    st = d.st
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, cf, base, lab, col in (
        (axes[0], trim(st.h_tot), trim(st.h_base), "high-skill", HIGH),
        (axes[1], trim(st.l_tot), trim(st.l_base), "low-skill", LOW),
    ):
        tt = t_of(cf)
        ax.plot(tt, cf, color=col, linewidth=1.8, label="reform path")
        ax.plot(tt, base, color=BASE_GRAY, linewidth=1.6, linestyle="--",
                label="baseline")
        ax.set_title(f"PT {lab} stock (level)")
        style_axis(ax, "workers (population units)")
        ax.legend(frameon=False, fontsize=9)
    save(fig, "fig1_stocks_levels.png")


def fig2_response_composition(d: Data) -> None:
    st = d.st
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    hr, lr = trim(st.h_ratio), trim(st.l_ratio)
    tt = t_of(hr)
    axes[0].plot(tt, 100 * (hr - 1), color=HIGH, linewidth=1.8,
                 label="high-skill")
    axes[0].plot(tt, 100 * (lr - 1), color=LOW, linewidth=1.8,
                 label="low-skill")
    axes[0].axhline(0, color=BASE_GRAY, linewidth=0.8)
    axes[0].set_title("stocks vs baseline (%)")
    style_axis(axes[0], "% deviation from baseline")
    axes[0].legend(frameon=False, fontsize=9)

    axes[1].plot(tt, 100 * trim(st.pref_share_of_h), color=HIGH, linewidth=1.8)
    axes[1].set_title("preferential pool, share of high-skill (%)")
    style_axis(axes[1], "%")
    save(fig, "fig2_response_composition.png")


def fig3_wages(d: Data) -> None:
    dg = d.dg
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    s = slice(1, T_SHOW)  # slots 2..199 (paths start moving at the reform)
    tt = np.arange(2, T_SHOW + 1)
    ax.plot(tt, 100 * (dg.What_h_path[s] - 1), color=HIGH, linewidth=1.8,
            label="high-skill")
    ax.plot(tt, 100 * (dg.What_l_path[s] - 1), color=LOW, linewidth=1.8,
            label="low-skill")
    ax.axhline(0, color=BASE_GRAY, linewidth=0.8)
    ax.set_title("PT wages vs baseline (%)")
    style_axis(ax, "% deviation from baseline")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "fig3_wages.png")


def fig4_taxes_budget(d: Data) -> None:
    dg = d.dg
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    s = slice(1, T_SHOW)
    tt = np.arange(2, T_SHOW + 1)
    ax = axes[0]
    ax.plot(tt, dg.tau_h_t_path[s], color=HIGH, linewidth=1.8,
            label="standard high-skill")
    ax.plot(tt, dg.tau_l_t_path[s], color=LOW, linewidth=1.8, label="low-skill")
    ax.axhline(TP, color=ACCENT, linewidth=1.8)
    ax.text(T_SHOW * 0.98, TP + 0.003, "preferential rate 0.13", ha="right",
            fontsize=9, color="#8a5a00")
    ax.axhline(RUN["tau_h"], color=BASE_GRAY, linewidth=0.8, linestyle="--")
    ax.axhline(RUN["tau_l"], color=BASE_GRAY, linewidth=0.8, linestyle="--")
    ax.text(T_SHOW * 0.98, RUN["tau_h"] + 0.001, "baseline rates", ha="right", fontsize=8,
            color=BASE_GRAY)
    ax.set_title("average tax rates")
    style_axis(ax, "average rate")
    ax.legend(frameon=False, fontsize=9, loc="center right")

    ax = axes[1]
    ax.plot(tt, 100 * (dg.lambda_hat_path[s] - 1), color=HIGH, linewidth=1.8)
    ax.axhline(0.0, color=BASE_GRAY, linewidth=0.8, linestyle="--")
    ax.set_title("budget-clearing schedule level $\\hat\\lambda_{t}$")
    style_axis(ax, "% deviation from baseline")
    save(fig, "fig4_taxes_budget.png")


def fig5_wedges(d: Data) -> None:
    dg = d.dg
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    s = slice(1, T_SHOW)
    tt = np.arange(2, T_SHOW + 1)

    ax = axes[0]
    ax.plot(tt, dg.f_pref_path[s], color=HIGH, linewidth=1.8, label="preferential")
    ax.plot(tt, dg.f_std_path[s], color=LOW, linewidth=1.8, label="standard h")
    ax.plot(tt, dg.f_l_path[s], color=BLOC["oth"], linewidth=1.8, label="low-skill")
    ax.set_title("tax wedge: retention-rate ratio $(1-\\tau')/(1-\\bar\\tau)$")
    style_axis(ax, "share of income kept, reform / baseline")
    ax.legend(frameon=False, fontsize=9)

    # after-tax income = pre-tax wage change x retention-rate ratio
    ax = axes[1]
    ax.plot(tt, 100 * (dg.What_h_path[s] * dg.f_pref_path[s] - 1), color=HIGH,
            linewidth=1.8, label="preferential")
    ax.plot(tt, 100 * (dg.What_h_path[s] * dg.f_std_path[s] - 1), color=LOW,
            linewidth=1.8, label="standard h")
    ax.plot(tt, 100 * (dg.What_l_path[s] * dg.f_l_path[s] - 1),
            color=BLOC["oth"], linewidth=1.8, label="low-skill")
    ax.axhline(0, color=BASE_GRAY, linewidth=0.8)
    ax.set_title("after-tax income vs baseline (wage $\\times$ retention)")
    style_axis(ax, "% deviation from baseline")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "fig5_wedges.png")


def fig6_welfare(d: Data) -> None:
    cev = d.cev
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    tt = np.arange(1, T_SHOW + 1)
    ax = axes[0]
    for b in BLOC:
        ax.plot(tt, getattr(cev, f"cev_{b}_h_pref_16")[:T_SHOW], color=BLOC[b],
                linewidth=1.8, label=f"{BLOC_NAME[b]}, employed")
    ax.plot(tt, getattr(cev, "cev_eu15_h_pref_33")[:T_SHOW], color=BLOC["eu15"],
            linewidth=1.4, linestyle="--", label="EU-15, non-employed")
    ax.axhline(0, color=BASE_GRAY, linewidth=0.8)
    ax.set_title("CEV: preferential high-skill")
    style_axis(ax, "CEV (% consumption)")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for b in BLOC:
        ax.plot(tt, getattr(cev, f"cev_{b}_h_std_16")[:T_SHOW], color=BLOC[b],
                linewidth=1.8, label=f"{BLOC_NAME[b]}, standard high-skill")
    ax.plot(tt, getattr(cev, "cev_eu15_l_16")[:T_SHOW], color=BASE_GRAY,
            linewidth=1.4, linestyle="--", label="EU-15, low-skill")
    ax.axhline(0, color=BASE_GRAY, linewidth=0.8)
    ax.set_title("CEV: standard-schedule and low-skill")
    style_axis(ax, "CEV (% consumption)")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "fig6_welfare_cev.png")


def fig7_homotopy(d: Data) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for col, (tp, st) in zip(SEQ_BLUES, sorted(d.st_homotopy.items(),
                                               reverse=True)):
        hr = trim(st.h_ratio)
        tt = t_of(hr)
        ax.plot(tt, 100 * (hr - 1), color=col, linewidth=1.8,
                label=f"$\\tau_P$ = {tp:.4g}")
    ax.axhline(0, color=BASE_GRAY, linewidth=0.8)
    ax.set_title("high-skill stock vs baseline (%), by preferential rate")
    style_axis(ax, "% deviation from baseline")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "fig7_homotopy.png")


def main() -> None:
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(HERE))
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--tau-h", type=float, default=0.24)
    ap.add_argument("--tau-l", type=float, default=0.18)
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    RUN.update(dir=Path(a.dir), prefix=a.prefix, tau_h=a.tau_h, tau_l=a.tau_l)
    OUT = Path(a.out)
    OUT.mkdir(exist_ok=True)
    d = load_all()
    fig1_stock_levels(d)
    fig2_response_composition(d)
    fig3_wages(d)
    fig4_taxes_budget(d)
    fig5_wedges(d)
    fig6_welfare(d)
    fig7_homotopy(d)


if __name__ == "__main__":
    main()
