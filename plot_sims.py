# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font

location = 'zimbabwe'


def plot_sims(df, start_year=2000, end_year=2025, percentile_pairs=[[.1, .99]], title='sim_plots'):
    """ Create quantile plots """
    set_font(size=20)
    fig, axes = pl.subplots(2, 5, figsize=(18, 9))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    ng_data = pd.read_csv(f'data/{location}_gonorrhea_data.csv')
    ct_data = pd.read_csv(f'data/{location}_chlamydia_data.csv')
    tv_data = pd.read_csv(f'data/{location}_trichomoniasis_data.csv')
    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    ng_data = ng_data.loc[(ng_data.year >= start_year) & (ng_data.year <= end_year)]
    ct_data = ct_data.loc[(ct_data.year >= start_year) & (ct_data.year <= end_year)]
    tv_data = tv_data.loc[(tv_data.year >= start_year) & (tv_data.year <= end_year)]
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]

    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
    dfplot['year'] = np.floor(np.round(dfplot.index, 1)).astype(int)

    # Population size
    pn = 0
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['n_alive'], color='k')
    resname = 'n_alive'
    x = np.unique(dfplot['year'])
    y = dfplot.groupby(by='year')[resname].mean()[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1])
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot.groupby(by='year')[resname].mean()[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot.groupby(by='year')[resname].mean()[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color(), label=f"{percentile_pair[0]:.0%}" + ' - ' + f"{percentile_pair[1]:.0%}")
    ax.set_title('Population size')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV infections
    pn = 5
    ax = axes[pn]
    resname = 'hiv.new_infections'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = np.unique(dfplot['year'])
    y = dfplot.groupby(by='year')[resname].sum()[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot.groupby(by='year')[resname].sum()[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot.groupby(by='year')[resname].sum()[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv.new_deaths'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = np.unique(dfplot['year'])
    y = dfplot.groupby(by='year')[resname].sum()[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot.groupby(by='year')[resname].sum()[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot.groupby(by='year')[resname].sum()[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # PLHIV
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], color='k')  # label='UNAIDS',
    resnames = {'Total': 'hiv.n_infected', 'Dx': 'hiv.n_diagnosed', 'Treated': 'hiv.n_on_art'}
    for rlabel, rname in resnames.items():
        x = np.unique(dfplot['year'])
        y = dfplot.groupby(by='year')[rname].mean()[(rname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label=rlabel)
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    resname = 'hiv.prevalence'
    ax.scatter(hiv_data.year, hiv_data[resname] * 100, label='Data', color='k')
    x = np.unique(dfplot['year'])
    y = dfplot.groupby(by='year')[resname].mean()[(resname, '50%')]
    line, = ax.plot(x, y * 100, label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot.groupby(by='year')[resname].mean()[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot.groupby(by='year')[resname].mean()[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + ".png", dpi=100)

    return fig


if __name__ == '__main__':

    df_stats = sc.loadobj('results/multi_res_stats.df')
    percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
    sims = plot_sims(df_stats, start_year=1990, percentile_pairs=percentile_pairs)

