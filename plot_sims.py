import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font, get_y

location = 'kenya'


def _load_data(start_year, end_year):
    calib = pd.read_csv(f'data/{location}_hiv_calib.csv')
    calib = calib.loc[(calib.time >= start_year) & (calib.time <= end_year)]
    n_art = pd.read_csv('data/n_art.csv')
    n_art = n_art.loc[(n_art.year >= start_year) & (n_art.year <= end_year)]
    return calib, n_art


def _plot_panel(ax, x_data, y_data, x_model, dfplot, mod_col, which,
                percentile_pairs, alphas, data_label='Data', model_label='Model',
                mult=1.0, ribbon_slice=None):
    ax.scatter(x_data, y_data, color='k', label=data_label)
    y = get_y(dfplot, which, mod_col) * mult
    line, = ax.plot(x_model, y, label=model_label)
    if which == 'multi':
        sl = ribbon_slice if ribbon_slice is not None else slice(None)
        for idx, pair in enumerate(percentile_pairs):
            yl = dfplot[(mod_col, f"{pair[0]:.0%}")] * mult
            yu = dfplot[(mod_col, f"{pair[1]:.0%}")] * mult
            ax.fill_between(x_model[sl], yl[sl], yu[sl], alpha=alphas[idx], facecolor=line.get_color())


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single',
                  percentile_pairs=[[.1, .99]], title='hiv_plots'):
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(18, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    calib, n_art = _load_data(start_year, end_year)
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]
    x = dfplot.index

    panels = [
        ('Population size',           calib.time, calib['n_alive'],                    'n_alive',              1.0,  None,        'Data'),
        ('PLHIV',                     calib.time, calib['hiv.n_infected'],             'hiv.n_infected',       1.0,  None,        'UNAIDS'),
        ('HIV prevalence 15-49 (%)',  calib.time, calib['hiv.prevalence_15_49'] * 100, 'hiv.prevalence_15_49', 100.0, None,        'UNAIDS'),
        ('New HIV infections/yr',     calib.time, calib['hiv.new_infections'],         'hiv.new_infections',   1.0,  None,        'UNAIDS'),
        ('HIV-related deaths',        calib.time, calib['hiv.new_deaths'],             'hiv.new_deaths',       1.0,  slice(0,-1), 'UNAIDS'),
        ('On ART',                    n_art.year, n_art['n_art'],                      'hiv.n_on_art',         1.0,  None,        'UNAIDS'),
    ]

    for ax, (ptitle, dx, dy, mcol, mult, sl, dlab) in zip(axes, panels):
        _plot_panel(ax, dx, dy, x, dfplot, mcol, which, percentile_pairs, alphas,
                    data_label=dlab, mult=mult, ribbon_slice=sl)
        ax.set_title(ptitle)
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
    axes[0].legend(frameon=False)

    sc.figlayout()
    sc.savefig(f'figures/{title}.png', dpi=100)
    return fig
