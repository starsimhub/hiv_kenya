# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font, get_y

location = 'kenya'


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile or individual plots of HIV epi dynamics """
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(18, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]

    pn = 0
    x = dfplot.index

    # Population size
    ax = axes[pn]
    resname = 'n_alive'
    ax.scatter(hiv_data.year, hiv_data[resname], color='k', label='Data')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='Modeled')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Population size')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # PLHIV
    ax = axes[pn]
    resname = 'hiv_n_infected'
    ax.scatter(hiv_data.year, hiv_data[resname], label='Data', color='k')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='PLHIV')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    resname = 'hiv_prevalence_15_49'
    ax.scatter(hiv_data.year, hiv_data[resname] * 100, label='Data', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y*100, label='Prevalence')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Infections
    ax = axes[pn]
    resname = 'hiv_new_infections'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV infections')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv_new_deaths'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV deaths')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # 90-90-90
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv_n_infected'], color='k')  # label='UNAIDS',
    resnames = {'PLHIV': 'hiv_n_infected', 'Dx': 'hiv_n_diagnosed', 'Treated': 'hiv_n_on_art'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y, label=rlabel)
        # if which == 'multi':
        #     for idx, percentile_pair in enumerate(percentile_pairs):
        #         yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
        #         yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
        #         ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Diagnosed and treated')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig

