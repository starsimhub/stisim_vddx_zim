# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font

location = 'zimbabwe'


def plot_hiv_sims(df, start_year=2000, end_year=2025, percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile plots of HIV """
    set_font(size=20)
    fig, axes = pl.subplots(2, 2, figsize=(8, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]

    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
    dfplot['year'] = np.floor(np.round(dfplot.index, 1)).astype(int)

    # HIV infections
    pn = 0
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
    sc.savefig("figures/" + title + str(start_year)+".png", dpi=100)

    return fig


def plot_sti_sims(df, start_year=2000, end_year=2025, percentile_pairs=[[.1, .99]], title='sti_plots'):
    """ Create quantile plots """
    set_font(size=20)
    fig, axes = pl.subplots(3, 4, figsize=(15, 12))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    ng_data = pd.read_csv(f'data/{location}_gonorrhea_data.csv')
    ct_data = pd.read_csv(f'data/{location}_chlamydia_data.csv')
    tv_data = pd.read_csv(f'data/{location}_trichomoniasis_data.csv')
    ng_data = ng_data.loc[(ng_data.year >= start_year) & (ng_data.year <= end_year)]
    ct_data = ct_data.loc[(ct_data.year >= start_year) & (ct_data.year <= end_year)]
    tv_data = tv_data.loc[(tv_data.year >= start_year) & (tv_data.year <= end_year)]

    df['year'] = np.floor(np.round(df.index, 1)).astype(int)
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    disease_map = {'ng': 'Gonorrhea', 'ct': 'Chlamydia', 'tv': 'Trich', 'vd': 'Other'}
    disease_data = {'ng': ng_data, 'ct': ct_data, 'tv': tv_data, 'vd': None}

    pn = 0

    # Incidence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resname = dname+'.new_infections'
        data = disease_data[dname]
        if data is not None:
            ax.scatter(data.year, data[resname], label='Data', color='k')
        resnames = {'Total': dname+'.new_infections', 'Symptomatic': dname+'.new_symptomatic' , 'Care seekers': dname+'.new_care_seekers'}
        for rlabel, rname in resnames.items():
            x = np.unique(dfplot['year'])
            y = dfplot.groupby(by='year')[resname].sum()[(resname, '50%')]
            line, = ax.plot(x[:-1], y[:-1], label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot.groupby(by='year')[rname].sum()[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot.groupby(by='year')[rname].sum()[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' burden')
        if pn == 7: ax.legend(frameon=False, prop={'size': 15})
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    # Burden
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resname = dname+'.n_infected'
        data = disease_data[dname]
        if data is not None:
            ax.scatter(data.year, data[resname], label='Data', color='k')
        resnames = {'Total': dname+'.n_infected', 'Symptomatic': dname+'.n_symptomatic'}
        for rlabel, rname in resnames.items():
            x = np.unique(dfplot['year'])
            y = dfplot.groupby(by='year')[rname].mean()[(rname, '50%')]
            line, = ax.plot(x[:-1], y[:-1], label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' burden')
        if pn == 7: ax.legend(frameon=False, prop={'size': 15})
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    # Prevalence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resnames = {'Total': dname+'.adult_prevalence', 'Symptomatic': dname+'.symp_adult_prevalence'}
        for rlabel, rname in resnames.items():
            x = np.unique(dfplot['year'])
            y = dfplot.groupby(by='year')[rname].mean()[(rname, '50%')]
            line, = ax.plot(x[:-1], y[:-1], label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' burden')
        if pn == 7: ax.legend(frameon=False, prop={'size': 15})
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + ".png", dpi=100)

    return fig

if __name__ == '__main__':

    df_stats = sc.loadobj('results/multi_res_stats.df')
    percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
    # plot_hiv_sims(df_stats, start_year=1980, percentile_pairs=percentile_pairs)
    plot_sti_sims(df_stats, start_year=1980, percentile_pairs=percentile_pairs)

    # Coinfection stats
    df_stats['year'] = np.floor(np.round(df_stats.index, 1)).astype(int)
    dfplot = df_stats.iloc[(df_stats.index >= 2020) & (df_stats.index <= 2030)]

    # sim.results.coinfection_stats.ng_only[si] + sim.results.coinfection_stats.ct_only[si] + sim.results.coinfection_stats.tv_only[si] + sim.results.coinfection_stats.vd_only[si]
