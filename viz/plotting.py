import matplotlib.pyplot as plt


def plot_tax_revenue(df, star_tau, G0, filename=None):
    """
    Plot Tax Revenue vs. 1 - lambda.

    Parameters:
    - df: DataFrame containing 'λ0' and 'TR' columns
    - star_tau: the value of 1 - λ* at the symmetric solution
    - G0: baseline government expenditure to mark horizontal line
    - filename: if provided, save the figure to this path; otherwise show it.
    """
    x = 1 - df['λ0']
    y = df['TR']

    plt.figure(figsize=(8,5))
    # Curve: solid black with hollow markers
    plt.plot(
        x, y,
        linestyle='-',
        marker='o',
        markersize=6,
        markevery=(0, 12),
        markerfacecolor='white',
        markeredgecolor='#0072B2',
        linewidth=1.5,
        color='#0072B2',
        label="Tax revenue"
    )

    # Vertical equilibrium line: dashed
    plt.axvline(
        x=star_tau,
        linestyle='--',
        linewidth=1.2,
        color='#E69F00',
        label=r"baseline $1-\lambda_0$"
    )

    # Horizontal G0 line: dotted
    plt.axhline(
        y=G0,
        linestyle=':',
        linewidth=1.2,
        color='#009E73',
        label=r"$G_0$"
    )

    plt.xlabel(r"$1-\lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    # axis limits
    #plt.xlim(0.15, 0.30)
    #plt.ylim(0.25, 0.45)
    
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


# def plot_comparison(df_sym, df_soe, star_tau, G0, filename=None):
#     """
#     Plot Tax Revenue curves for Multilateral (df_sym) vs. Unilateral/SOE (df_soe).

#     df_sym: DataFrame from the symmetric sweep.
#     df_soe: DataFrame from the SOE sweep.
#     star_tau: scalar = 1 - λ* from the symmetric solution.
#     G0: baseline government spending level.
#     """
#     x_sym = 1 - df_sym['λ0']
#     y_sym = df_sym['TR']

#     x_soe = 1 - df_soe['λ0']
#     y_soe = df_soe['TR']

#     plt.figure(figsize=(8,5))

#     # Multilateral (symmetric): dashed line with square markers
#     plt.plot(
#         x_sym, y_sym,
#         linestyle='-',
#         marker='o',
#         markersize=6,
#         markerfacecolor='white',
#         markeredgecolor='#0072B2',
#         color='#0072B2',
#         linewidth=1.5,
#         label="Multilateral (symmetric)"
#     )

#     # Unilateral (SOE): solid line with circle markers
#     plt.plot(
#         x_soe, y_soe,
#         linestyle='-',
#         marker='s',
#         markersize=6,
#         markerfacecolor='white',
#         markeredgecolor='#D55E00',
#         color='#D55E00',
#         linewidth=1.5,
#         label="Unilateral (SOE)"
#     )

#     # Equilibrium 1 - λ*: vertical dotted line
#     plt.axvline(
#         x=star_tau,
#         linestyle='--',
#         linewidth=1.2,
#         color='#E69F00',
#         label=r"baseline $1-\lambda_0$"
#     )

#     # Baseline G0: horizontal dash-dot line
#     plt.axhline(
#         y=G0,
#         linestyle=':',
#         linewidth=1.2,
#         color='#009E73',
#         label=r"$G_0$"
#     )

#     plt.xlabel(r"$1-\lambda_0$")
#     plt.ylabel("Tax revenue")
#     plt.title("Tax revenue: Unilateral vs. Multilateral")
#     plt.legend(frameon=False)
#     plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
#     # axis limits
#     #plt.xlim(0.15, 0.30)
#     #plt.ylim(0.25, 0.45)
    
#     plt.tight_layout()

#     if filename:
#         plt.savefig(filename)
#         plt.close()
#     else:
#         plt.show()

def plot_unilateral_laffer(df_mult, df_soe, df_dd, star_tau, G0, filename=None,
                           dd_label="Unilateral, differentiated destinations"):
    """
    Tax revenue vs. 1-λ0: multilateral (symmetric) vs the two unilateral cases —
    perfect substitutes (common shock) and differentiated destinations.

    Parameters:
    - df_mult: multilateral (symmetric) sweep ('λ0','TR' columns)
    - df_soe: unilateral sweep, perfect-substitutes case
    - df_dd: unilateral sweep, differentiated destinations
    - star_tau: 1 - λ* at the calibrated baseline
    - G0: government expenditure (horizontal line)
    """
    plt.figure(figsize=(8,5))

    plt.plot(
        1 - df_mult['λ0'], df_mult['TR'],
        linestyle='-',
        marker='o', markersize=6, markevery=(0, 12),
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label="Multilateral (symmetric)"
    )

    plt.plot(
        1 - df_soe['λ0'], df_soe['TR'],
        linestyle='-',
        marker='s', markersize=6, markevery=(4, 12),
        markerfacecolor='white', markeredgecolor='#006400',
        color='#006400', linewidth=1.5,
        label="Unilateral, perfect substitutes"
    )

    plt.plot(
        1 - df_dd['λ0'], df_dd['TR'],
        linestyle='-',
        marker='^', markersize=6, markevery=(8, 12),
        markerfacecolor='white', markeredgecolor='#D55E00',
        color='#D55E00', linewidth=1.5,
        label=dd_label
    )

    plt.axvline(
        x=star_tau, linestyle='--', linewidth=1.2,
        color='#E69F00', label=r"baseline $1-\lambda_0$"
    )

    plt.axhline(
        y=G0, linestyle=':', linewidth=1.2,
        color='#009E73', label=r"$G_0$"
    )

    plt.xlabel(r"$1-\lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')

    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_diff_dest_laffer(df_base, df_pref, star_tau, G0, filename=None,
                          df_soe_pref=None):
    """
    Tax revenue vs. 1-λ0 for the differentiated-destinations small country,
    with and without the preferential tax.

    Parameters:
    - df_base: sweep without the preferential tax ('λ0','TR' columns)
    - df_pref: sweep with the preferential tax
    - star_tau: 1 - λ* at the calibrated baseline
    - G0: government expenditure (horizontal line)
    - df_soe_pref: optional overlay, the SOE preferential-tax sweep
      (perfect-substitutes limit)
    """
    plt.figure(figsize=(8,5))

    plt.plot(
        1 - df_base['λ0'], df_base['TR'],
        linestyle='-',
        marker='o', markersize=6, markevery=(0, 12),
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label="No preferential tax"
    )

    plt.plot(
        1 - df_pref['λ0'], df_pref['TR'],
        linestyle='-',
        marker='s', markersize=6, markevery=(4, 12),
        markerfacecolor='white', markeredgecolor='#D55E00',
        color='#D55E00', linewidth=1.5,
        label="Preferential tax"
    )

    if df_soe_pref is not None:
        plt.plot(
            1 - df_soe_pref['λ0'], df_soe_pref['TR'],
            linestyle='--',
            marker='^', markersize=6, markevery=(8, 12),
            markerfacecolor='white', markeredgecolor='#CC79A7',
            color='#CC79A7', linewidth=1.5,
            label="Preferential tax, perfect-substitutes limit"
        )

    plt.axvline(
        x=star_tau, linestyle='--', linewidth=1.2,
        color='#E69F00', label=r"baseline $1-\lambda_0$"
    )

    plt.axhline(
        y=G0, linestyle=':', linewidth=1.2,
        color='#009E73', label=r"$G_0$"
    )

    plt.xlabel(r"$1-\lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')

    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_theta_laffer(dfs_pref, df_base, star_tau, G0, filename=None):
    """
    Tax revenue vs. 1-λ0 in the small country with the preferential tax, for several
    values of θ (choice among foreign destinations); the no-preferential-tax curve is
    drawn for reference.

    Parameters:
    - dfs_pref: dict {label: sweep with the preferential tax ('λ0','TR' columns)},
      in the order to be drawn
    - df_base: sweep without the preferential tax
    - star_tau: 1 - λ* at the calibrated baseline
    - G0: government expenditure (horizontal line)
    """
    plt.figure(figsize=(8,5))

    plt.plot(
        1 - df_base['λ0'], df_base['TR'],
        linestyle='-', color='0.35', linewidth=1.5,
        marker='o', markersize=6, markevery=(0, 12),
        markerfacecolor='white', markeredgecolor='0.35',
        label="No preferential tax"
    )

    styles = [('#0072B2', 's'), ('#D55E00', '^'), ('#CC79A7', 'D'),
              ('#006400', 'v'), ('#56B4E9', 'P')]
    for k, (lab, df) in enumerate(dfs_pref.items()):
        col, mk = styles[k % len(styles)]
        plt.plot(
            1 - df['λ0'], df['TR'],
            linestyle='-', color=col, linewidth=1.5,
            marker=mk, markersize=6, markevery=(2 + 2*k, 12),
            markerfacecolor='white', markeredgecolor=col,
            label=lab
        )

    plt.axvline(
        x=star_tau, linestyle='--', linewidth=1.2,
        color='#E69F00', label=r"baseline $1-\lambda_0$"
    )

    plt.axhline(
        y=G0, linestyle=':', linewidth=1.2,
        color='#009E73', label=r"$G_0$"
    )

    plt.xlabel(r"$1-\lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')

    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_g_per_capita_laffer(df_base, df_pref, star_tau, G0, eq_fixed, eq_gpc,
                             filename=None, share=1.0):
    """
    Government spending proportional to the resident population N of country 0 (small
    country with differentiated destinations or in the perfect-substitutes limit, or the
    adopting large country). Left: tax revenue and required spending G0·N vs. 1-λ0, with
    and without the preferential tax. Right: tax revenue per resident vs. G0.

    Parameters:
    - df_base, df_pref: sweeps without / with the preferential tax ('λ0','TR',
      'p00l','p10l','p00h','p10h' columns)
    - star_tau: 1 - λ* at the calibrated baseline
    - G0: government expenditure at the baseline (N = 1)
    - eq_fixed, eq_gpc: equilibria with the preferential tax, (λ0, N), under G = G0
      and G = G0·N; None skips the marker (no budget-clearing equilibrium)
    - share: share s of spending proportional to the population; with s < 1 the required
      spending is G0·[(1 − s) + s·N], drawn dashed in both panels
    """
    scale = lambda N: (1.0 - share) + share * N
    g_tex = (r"$G_0 N$" if share == 1.0
             else rf"$G_0[{1 - share:.2f} + {share:.2f} N]$")
    stocks = ['p00l', 'p10l', 'p00h', 'p10h']
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, per_resident in zip(axes, (False, True)):
        for df, col, mk, mev, lab in (
                (df_base, '#0072B2', 'o', (0, 12), "No preferential tax"),
                (df_pref, '#D55E00', 's', (4, 12), "Preferential tax")):
            N = df[stocks].sum(axis=1)
            ax.plot(
                1 - df['λ0'], df['TR'] / N if per_resident else df['TR'],
                linestyle='-', color=col, linewidth=1.5,
                marker=mk, markersize=6, markevery=mev,
                markerfacecolor='white', markeredgecolor=col, label=lab
            )
            if not per_resident:
                ax.plot(1 - df['λ0'], G0 * scale(N), linestyle='--', color=col, linewidth=1.1,
                        label=g_tex + ", " + lab[0].lower() + lab[1:])
            elif share != 1.0:
                ax.plot(1 - df['λ0'], G0 * scale(N) / N, linestyle='--', color=col,
                        linewidth=1.1, label=g_tex + r"$/N$, " + lab[0].lower() + lab[1:])

        ax.axvline(x=star_tau, linestyle='--', linewidth=1.2,
                   color='#E69F00', label=r"baseline $1-\lambda_0$")
        ax.axhline(y=G0, linestyle=':', linewidth=1.2, color='#009E73', label=r"$G_0$")

        for (λ0, N), mk, lab in ((eq_fixed or (None, None), 'o', r"equilibrium, $G = G_0$"),
                                 (eq_gpc or (None, None), 'D', r"equilibrium, $G = $" + g_tex)):
            if (λ0, N) == (None, None):
                continue
            G = G0 * scale(N) if mk == 'D' else G0
            ax.plot(1 - λ0, G / N if per_resident else G, mk, markersize=9,
                    markerfacecolor='none', markeredgecolor='k', label=lab)

        ax.set_xlabel(r"$1-\lambda_0$")
        ax.set_ylabel("Tax revenue per resident" if per_resident else "Tax revenue")
        ax.legend(frameon=False, fontsize=9)
        ax.grid(True, linestyle=':', linewidth=0.6, color='0.75')

    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_comparison(df_sym, df_soe, star_tau, G0, filename=None,
                    mark_equilibrium=False,  # new
                    eq_kwargs=None):         # optional style overrides
    """
    Plot Tax Revenue curves for Multilateral (df_sym) vs. Unilateral/SOE (df_soe).

    mark_equilibrium: if True, annotate the intersection of SOE curve with y=G0.
    eq_kwargs: optional dict with keys:
        - text (str): label text (default: "Equilibrium")
        - circle_edge (str): circle edge color (default: '#B22222')
        - circle_lw (float): circle linewidth (default: 2.0)
        - arrow_color (str): arrow color (default: 'black')
        - xytext (tuple): text position; if None, auto-offset
    """
    import numpy as np
    import matplotlib.pyplot as plt

    x_sym = 1 - df_sym['λ0']
    y_sym = df_sym['TR']

    x_soe = 1 - df_soe['λ0']
    y_soe = df_soe['TR']

    plt.figure(figsize=(8,5))

    # Unilateral, no preferential tax
    plt.plot(
        x_sym, y_sym,
        linestyle='-',
        marker='o', markersize=6, markevery=(0, 12),
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label="No preferential tax"
    )

    # Unilateral (SOE)
    plt.plot(
        x_soe, y_soe,
        linestyle='-',
        marker='s', markersize=6, markevery=(6, 12),
        markerfacecolor='white', markeredgecolor='#006400',
        color='#006400', linewidth=1.5,
        label="Preferential tax"
    )

    # Equilibrium 1 - λ*: vertical line
    plt.axvline(
        x=star_tau, linestyle='--', linewidth=1.2,
        color='#E69F00', label=r"baseline $1-\lambda_0$"
    )

    # Baseline G0: horizontal line
    plt.axhline(
        y=G0, linestyle=':', linewidth=1.2,
        color='#009E73', label=r"$G_0$"
    )

    # ----- Optional equilibrium annotation (SOE curve ∩ G0) -----
    if mark_equilibrium:
        # defaults + user overrides
        defaults = dict(
            text="Equilibrium",
            circle_edge='#B22222',
            circle_lw=2.0,
            arrow_color='black',
            xytext=None
        )
        if eq_kwargs:
            defaults.update(eq_kwargs)

        # Find intersection by linear interpolation
        ydiff = y_soe - G0
        # indices where sign changes
        cross_idx = np.where(np.sign(ydiff[:-1]) * np.sign(ydiff[1:]) <= 0)[0]
        if len(cross_idx) > 0:
            i = int(cross_idx[0])
            x0, x1 = float(x_soe.iloc[i]), float(x_soe.iloc[i+1])
            y0, y1 = float(y_soe.iloc[i]), float(y_soe.iloc[i+1])
            # handle flat segment edge-case
            if y1 != y0:
                t = (G0 - y0) / (y1 - y0)
            else:
                t = 0.0
            x_star = x0 + t * (x1 - x0)
        else:
            # fallback: closest point
            j = int(np.argmin(np.abs(ydiff.values)))
            x_star = float(x_soe.iloc[j])

        y_star = float(G0)

        ax = plt.gca()
        # auto offset if not provided
        if defaults["xytext"] is None:
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            dx = 0.05 * (xmax - xmin)
            dy = 0.05 * (ymax - ymin)
            xytext = (x_star + dx, y_star + dy)
        else:
            xytext = defaults["xytext"]

        # draw an open circle
        ax.scatter([x_star], [y_star],
                   facecolors='none',
                   edgecolors=defaults["circle_edge"],
                   linewidths=defaults["circle_lw"],
                   s=300, zorder=5)

        # annotation with arrow
        ax.annotate(
            defaults["text"],
            xy=(x_star, y_star),
            xytext=(x_star - 0.03, y_star + 0.01),  # only a small step left/up
            textcoords='data',
            arrowprops=dict(arrowstyle='->', lw=1.2, color=defaults["arrow_color"]),
            ha='right', va='bottom'
        )

    # ------------------------------------------------------------

    plt.xlabel(r"$1-\lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()



def plot_tax_vs_tau10h(df, G0, tau10h_star=None, filename=None):
    """
    Plot Tax Revenue vs. high-skill immigrant tax rate τ10h.

    Parameters:
    - df: DataFrame with columns 'tau10h' and 'TR'
    - G0: baseline government expenditure
    - tau10h_star: optional equilibrium vertical line value
    - filename: path to save the figure
    """
    x = df['tau10h']
    y = df['TR']

    plt.figure(figsize=(8,5))
    # Curve: solid black with hollow markers
    plt.plot(
        x, y,
        linestyle='-',
        marker='o',
        markersize=6,
        markerfacecolor='white',
        markeredgecolor='black',
        color='black',
        linewidth=1.5,
        label='Tax revenue'
    )

    # Baseline G0: horizontal dotted line
    plt.axhline(
        y=G0,
        linestyle=':',
        linewidth=1.2,
        color='#009E73',
        label=r'$G_0$'
    )

    # Optional equilibrium τ10h*: vertical dashed line
    if tau10h_star is not None:
        plt.axvline(
            x=tau10h_star,
            linestyle='--',
            linewidth=1.2,
            color='#E69F00',
            label=fr'$\bar{{\tau}}_{{10}}^h = {tau10h_star:.2f}$'
        )

    plt.xlabel(r'$\tau_{10}^h$ (high-skill immigrant tax)')
    plt.ylabel('Total tax revenue')
    plt.title('Tax revenue vs. high-skill immigrant tax')
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()
        
        
def plot_two_country_comparison(df_pe, df_sym, pi_eq, G0, filename=None):
   """
   Plot Tax Revenue curves for Country 0, symmetric multilateral, and Country 1.

   Parameters
   ----------
   df_pe : DataFrame
       Output from sweep_equilibrium_2country in PE mode, with columns 'λ0', 'TR0', 'TR1'.
   df_sym : DataFrame
       Output from sweep_equilibrium for the symmetric (multilateral) case, with columns 'λ0', 'TR'.
   pi_eq : array-like
       The full symmetric equilibrium solution from one-country GE: pi_ge_eq.
   G0 : float
       Baseline government expenditure marker.
   filename : str, optional
       Path to save the figure. If None, shows the plot.
   """
   x = 1 - df_pe['λ0']
   rev0 = df_pe['TR0']
   rev1 = df_pe['TR1']
   if 'TR' in df_sym.columns:
        rev_sym = df_sym['TR']
   elif 'TR0' in df_sym.columns:
        rev_sym = df_sym['TR0']

   # equilibrium reference
   tau_star = 1 - pi_eq[4]

   plt.figure(figsize=(8,5))
   
    # Multilateral: dashed + squares
   plt.plot(
       x, rev_sym,
       linestyle='None',
       marker='+', markersize=15, markevery=(0, 12),
       markerfacecolor='white', markeredgecolor='#CC79A7',
       color='#CC79A7', linewidth=1.5,
       label="Country 0, unilateral adoption"
   )
   
   # Country 0 (PE): solid + circles
   plt.plot(
       x, rev0,
       linestyle='--',
       marker='o', markersize=6, markevery=(4, 12),
       markerfacecolor='white', markeredgecolor='#0072B2',
       color='#0072B2', linewidth=1.5,
       label="Country 0"
   )
    
    # Country 1 (PE): dash-dot + triangles
   plt.plot(
       x, rev1,
       linestyle='-.',
       marker='s', markersize=6, markevery=(8, 12),
       markerfacecolor='white', markeredgecolor='#D55E00',
       color='#D55E00', linewidth=1.5,
       label="Country 1"
   )
    
    # Reference lines (monochrome)
   plt.axvline(x=tau_star, linestyle='--', linewidth=1.2, color='#E69F00', label=r"baseline $1-\lambda_0$")
   plt.axhline(y=G0,       linestyle=':', linewidth=1.2, color='#009E73', label=r'$G_0$')
    
   plt.xlabel(r"$1 - \lambda_0$")
   plt.ylabel("Tax revenue")
   plt.legend(frameon=False)
   plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
   # axis limits
   #plt.xlim(0.15, 0.30)
   #plt.ylim(0.25, 0.45)
    
   plt.tight_layout()
    
   if filename:
       plt.savefig(filename)
       plt.close()
   else:
       plt.show()
    
       
    
def plot_two_country_unilateral(df, pi_eq, G0, filename=None):
    """
    Plot Country 0 and Country 1 unilateral tax‐revenue curves.

    df must have columns 'λ0', 'TR0', 'TR1'.
    pi_eq is the one‐country GE solution array (so we can mark the equilibrium τ).
    """
    x = 1 - df['λ0']
    rev0 = df['TR0']
    rev1 = df['TR1']
    tau_star = 1 - pi_eq[4]

    plt.figure(figsize=(8,5))
    # Country 0: solid + circles
    plt.plot(
        x, rev0,
        linestyle='-',
        marker='o', markersize=6, markevery=(0, 12),
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label="Country 0"
    )

    # Country 1: dashed + squares
    plt.plot(
        x, rev1,
        linestyle='--',
        marker='s', markersize=6, markevery=(6, 12),
        markerfacecolor='white', markeredgecolor='#D55E00',
        color='#D55E00', linewidth=1.5,
        label="Country 1"
    )

    # Reference lines (black only)
    plt.axvline(x=tau_star, linestyle='--', linewidth=1.2, color='#E69F00', label=r"baseline $1-\lambda_0$")
    plt.axhline(y=G0,       linestyle=':', linewidth=1.2, color='#009E73', label=r'$G_0$')

    plt.xlabel(r"$1 - \lambda_0$")
    plt.ylabel("Tax revenue")
    plt.legend(frameon=False)
    plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    # axis limits
    #plt.xlim(0.15, 0.30)
    #plt.ylim(0.25, 0.45)
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()




def plot_tax_rates_vs_lambda(df, pi_eq, filename=None):
    """
    Plot tax rates vs 1 - lambda_0, with vertical line at equilibrium.
    
    Parameters:
    - df: DataFrame with tax rates columns and 'λ0'
    - pi_eq: equilibrium solution array (e.g. pi_ge_eq)
    - filename: optional save path
    """
    lambda_vals = df['λ0']
    x_val = 1 - pi_eq[4]  # equilibrium location
    
    plt.figure(figsize=(8,5))
    plt.plot(1 - lambda_vals, df['tau00l'], label=r'$\tau_{00}^l$', marker='o')
    plt.plot(1 - lambda_vals, df['tau10l'], label=r'$\tau_{10}^l$', marker='s')
    plt.plot(1 - lambda_vals, df['tau00h'], label=r'$\tau_{00}^h$', marker='^')
    plt.plot(1 - lambda_vals, df['tau10h'], label=r'$\tau_{10}^h$', marker='d')
    plt.axvline(x=x_val, color='red', linestyle='--', label=fr'$\lambda_0^* = {x_val:.2f}$')
    plt.xlabel(r'$1-\lambda_0$')
    plt.ylabel('Tax Rates')
    plt.title(r'Tax Rates vs $1-\lambda_0$')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_net_incomes_vs_lambda(df, pi_eq, filename=None):
    """
    Plot net incomes vs 1 - lambda_0, with vertical line at equilibrium.
    
    Parameters:
    - df: DataFrame with income columns and 'λ0'
    - pi_eq: equilibrium solution array (e.g. pi_ge_eq)
    - filename: optional save path
    """
    lambda_vals = df['λ0']
    x_val = 1 - pi_eq[4]  # equilibrium location
    
    plt.figure(figsize=(8,5))
    plt.plot(1 - lambda_vals, df['inc00l'], label=r'$y_{00}^l$', marker='o')
    plt.plot(1 - lambda_vals, df['inc10l'], label=r'$y_{10}^l$', marker='s')
    plt.plot(1 - lambda_vals, df['inc00h'], label=r'$y_{00}^h$', marker='^')
    plt.plot(1 - lambda_vals, df['inc10h'], label=r'$y_{10}^h$', marker='d')
    plt.axvline(x=x_val, color='red', linestyle='--', label=fr'$\lambda_0^* = {x_val:.2f}$')
    plt.xlabel(r'$1-\lambda_0$')
    plt.ylabel('Net incomes')
    plt.title(r'Net incomes vs $1-\lambda_0$')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.close()
    else:
        plt.show()


def plot_taul_mlt_vs_soe(df_mlt, df_soe, G0, tau_l_base, filename=None):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,5))
    plt.plot(df_mlt['tau_l'], df_mlt['TR'],
             linestyle='--', marker='o', markersize=6,
             markerfacecolor='white', markeredgecolor='#0072B2',
             color='#0072B2', linewidth=1.5, label='Multilateral')
    plt.plot(df_soe['tau_l'], df_soe['TR'],
             linestyle='-', marker='s', markersize=6,
             markerfacecolor='white', markeredgecolor='#006400',
             color='#006400', linewidth=1.5, label='Unilateral (SOE)')
    plt.axhline(G0, linestyle=':', linewidth=1.2, color='#009E73', label=r'$G_0$')
    plt.axvline(tau_l_base, linestyle='--', linewidth=1.2, color='#E69F00', label=r'baseline $\tau_l$')
    plt.xlabel(r'$\tau_l$ (progressivity)'); plt.ylabel('Tax revenue')
    plt.title(r'Tax revenue vs. $\tau_l$ (λ fixed at baseline)')
    plt.legend(frameon=False); plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    plt.tight_layout()
    if filename: plt.savefig(filename, dpi=300); plt.close()
    else: plt.show()

def plot_taul_soe_pref(df_soe, df_pref, G0, tau_l_base, pref_tax, filename=None):
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,5))
    plt.plot(df_soe['tau_l'], df_soe['TR'],
             linestyle='--', marker='o', markersize=6,
             markerfacecolor='white', markeredgecolor='#0072B2',
             color='#0072B2', linewidth=1.5, label='SOE (no pref tax)')
    plt.plot(df_pref['tau_l'], df_pref['TR'],
             linestyle='-', marker='s', markersize=6,
             markerfacecolor='white', markeredgecolor='#D55E00',
             color='#D55E00', linewidth=1.5, label=fr'SOE ($\tau_{{10}}^h={pref_tax:.2f}$)')
    plt.axhline(G0, linestyle=':', linewidth=1.2, color='#009E73', label=r'$G_0$')
    plt.axvline(tau_l_base, linestyle='--', linewidth=1.2, color='#E69F00', label=r'baseline $\tau_l$')
    plt.xlabel(r'$\tau_l$ (progressivity)'); plt.ylabel('Tax revenue')
    plt.title(r'SOE Laffer vs $\tau_l$ (λ fixed) — Preferential tax')
    plt.legend(frameon=False); plt.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    plt.tight_layout()
    if filename: plt.savefig(filename, dpi=300); plt.close()
    else: plt.show()


# --- τ_l sweeps in two large economies ---


def plot_two_country_taul_case2(df, tauL_baseline, G_line, filename=None, show=False):
    """
    Case 2: Country 0 uses preferential tax on high-skill immigrants; Country 1 does not.
    x-axis = tauL (progressivity), λ fixed. Plots TR0 and TR1.
    """
    x = df['tauL0']
    fig, ax = plt.subplots(figsize=(8,5))
    # Country 0: solid + circles (blue)
    ax.plot(
        x, df['TR0'],
        linestyle='-',
        marker='o', markersize=6,
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label='Country 0 (pref on $h$-immigrants)'
    )
    # Country 1: dashed + squares (orange)
    ax.plot(
        x, df['TR1'],
        linestyle='--',
        marker='s', markersize=6,
        markerfacecolor='white', markeredgecolor='#D55E00',
        color='#D55E00', linewidth=1.5,
        label='Country 1 (no pref)'
    )
    # Reference lines
    ax.axvline(x=tauL_baseline, linestyle='--', linewidth=1.2, color='#E69F00', label=r'baseline $\tau_l$')
    ax.axhline(y=G_line,       linestyle=':',  linewidth=1.2, color='#009E73', label=r'$G_0$')

    ax.set_xlabel(r'$\tau_l$ (progressivity)')
    ax.set_ylabel('Tax revenue')
    ax.set_title(r'Tax revenue vs. $\tau_l$ (λ fixed; pref only in Country 0)')
    ax.legend(frameon=False)
    ax.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    fig.tight_layout()
    if filename: fig.savefig(filename, dpi=300)
    if show: plt.show()
    plt.close(fig)


def plot_two_country_taul_case3(df, tauL_baseline, G_line, filename=None, show=False):
    """
    Case 3: Both countries adopt preferential tax on high-skill immigrants.
    x-axis = tauL (progressivity), λ fixed. Plots TR0 and TR1.
    """
    x = df['tauL0']
    fig, ax = plt.subplots(figsize=(8,5))
    # Country 0: solid + circles (blue)
    ax.plot(
        x, df['TR0'],
        linestyle='-',
        marker='o', markersize=6,
        markerfacecolor='white', markeredgecolor='#0072B2',
        color='#0072B2', linewidth=1.5,
        label='Country 0 (bilateral pref)'
    )
    # Country 1: dashed + squares (orange)
    ax.plot(
        x, df['TR1'],
        linestyle='--',
        marker='s', markersize=6,
        markerfacecolor='white', markeredgecolor='#D55E00',
        color='#D55E00', linewidth=1.5,
        label='Country 1 (bilateral pref)'
    )
    # Reference lines
    ax.axvline(x=tauL_baseline, linestyle='--', linewidth=1.2, color='#E69F00', label=r'baseline $\tau_l$')
    ax.axhline(y=G_line,       linestyle=':',  linewidth=1.2, color='#009E73', label=r'$G_0$')

    ax.set_xlabel(r'$\tau_l$ (progressivity)')
    ax.set_ylabel('Tax revenue')
    ax.set_title(r'Tax revenue vs. $\tau_l$ (λ fixed; bilateral pref)')
    ax.legend(frameon=False)
    ax.grid(True, linestyle=':', linewidth=0.6, color='0.75')
    fig.tight_layout()
    if filename: fig.savefig(filename, dpi=300)
    if show: plt.show()
    plt.close(fig)


