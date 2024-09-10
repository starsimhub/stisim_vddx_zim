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

    # HIV infections
    pn = 0
    ax = axes[pn]
    resname = 'hiv.new_infections'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv.new_deaths'
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
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
        x = dfplot.index
        y = dfplot[(rname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label=rlabel)
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
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
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x, y * 100, label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year)+".png", dpi=100)

    return fig


def plot_sti_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='sti_plots'):
    """ Create quantile plots """
    set_font(size=30)
    fig, axes = pl.subplots(3, 4, figsize=(25, 12))
    axes = axes.ravel()
    if which == 'multi': alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    ng_data = pd.read_csv(f'data/{location}_gonorrhea_data.csv')
    ct_data = pd.read_csv(f'data/{location}_chlamydia_data.csv')
    tv_data = pd.read_csv(f'data/{location}_trichomoniasis_data.csv')
    ng_data = ng_data.loc[(ng_data.year >= start_year) & (ng_data.year <= end_year)]
    ct_data = ct_data.loc[(ct_data.year >= start_year) & (ct_data.year <= end_year)]
    tv_data = tv_data.loc[(tv_data.year >= start_year) & (tv_data.year <= end_year)]

    # df['year'] = np.floor(np.round(df.index, 1)).astype(int)
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    disease_map = {'ng': 'Gonorrhea', 'ct': 'Chlamydia', 'tv': 'Trich', 'bv': 'Other'}
    disease_data = {'ng': ng_data, 'ct': ct_data, 'tv': tv_data, 'bv': None}

    pn = 0
    def get_y(df, which, rname):
        if which == 'single': y = df[rname]
        elif which == 'multi': y = dfplot[(rname, '50%')]
        return y

    # Incidence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resname = dname+'.new_infections'
        data = disease_data[dname]
        if data is not None:
            ax.scatter(data.year, data[resname], label='Data', color='k')
        resnames = {'Total': dname+'.new_infections', 'Symptomatic': dname+'.new_symptomatic', 'Care seekers': dname+'.new_care_seekers'}
        for rlabel, rname in resnames.items():
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())

        ax.set_title(dlabel+' incidence')
        if pn == 2: ax.legend(frameon=False, prop={'size': 20})
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
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' burden')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    # Prevalence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resnames = {'Total': dname+'.adult_prevalence', 'Symptomatic': dname+'.symp_adult_prevalence'}
        if dname == 'ng':
            data = disease_data[dname]
            ax.scatter(data.year, 100*data['ng.adult_prevalence'], label='Data', color='k')
        for rlabel, rname in resnames.items():
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y*100, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl* 100, yu* 100, alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' prevalence (%)')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig


def plot_sti_tx(df, start_year=2000, end_year=2020):
    set_font(size=24)
    legend_font=16
    fig, axes = pl.subplots(2, 4, figsize=(25, 8))
    axes = axes.ravel()
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    pn = 0

    # Care seeking plot
    ax = axes[pn]
    x = dfplot.index
    y1 = dfplot['total_symptomatic.new_symptoms']
    y2 = dfplot['syndromicmgmt.new_care_seekers']
    ax.plot(x, y1, label='Discharge incidence')
    ax.plot(x, y2, label='Number seeking care')
    ax.set_title('Symptomatic care')
    ax.legend(frameon=False, prop={'size': legend_font})
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Care seeker epidemiology
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['syndromicmgmt.new_sti1'],
        dfplot['syndromicmgmt.new_sti2'],
        dfplot['syndromicmgmt.new_sti3'],
        dfplot['syndromicmgmt.new_sti4'],
    ]
    labels = ["1 STI ", "2 STIs", "3 STIs", "4 STIs"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels)
    ax.set_title('Coinfection among care seekers')
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Care seeker management
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['syndromicmgmt.new_tx0'],
        dfplot['syndromicmgmt.new_tx1'],
        dfplot['syndromicmgmt.new_tx2'],
        dfplot['syndromicmgmt.new_tx3'],
    ]
    labels = ["0", "1", "2", "3"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels)
    ax.set_title('Treatments prescribed')
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Overtreatment
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['ng_tx.new_treated_unnecessary'],
        dfplot['ct_tx.new_treated_unnecessary'],
        dfplot['metronidazole.new_treated_unnecessary'],
    ]
    labels = ["Ceftriaxone", "Doxycycline", "Metronidazole"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels)
    ax.set_title('Treatments overprescribed')
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Treatment stats by disease
    for disease in ['ng', 'ct', 'tv', 'bv']:
        ax = axes[pn]
        x = dfplot.index
        Y = [
            dfplot[disease+'.new_treated_success'],
            dfplot[disease+'.new_treated_failure'],
            dfplot[disease+'.new_treated_unnecessary'],
        ]
        labels = ["Treatment success", "Treatment failure", "Overtreatment"]
        ax.stackplot(x, *Y, baseline='zero', labels=labels)
        ax.set_title('Care seekers treated for ' + disease.upper())
        ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    sc.figlayout()
    sc.savefig("figures/sti_tx.png", dpi=100)
    return


def print_results(df):

    year = 2020
    metric = 'mean'
    dfg = df.groupby(by='year')

    output = ''
    r1 = (dfg[('ng.new_symptomatic', metric)][year] +
          dfg[('ct.new_symptomatic', metric)][year] +
          dfg[('tv.new_symptomatic', metric)][year] +
          dfg[('vd.new_symptomatic', metric)][year] )
    r2 = dfg[('syndromicmgmt.care_seekers', metric)][year]
    output += f'Total newly symptomatic: {int(r1)}\n'
    output += f'Total seeking care: {int(r2)}\n'

    r11 = (dfg[('syndromicmgmt.ng_only', metric)][year])
    r12 = (dfg[('syndromicmgmt.ct_only', metric)][year])
    r13 = (dfg[('syndromicmgmt.tv_only', metric)][year])
    r14 = (dfg[('syndromicmgmt.vd_only', metric)][year])
    output += f'NG only: {int(r11)}\n'
    output += f'CT only: {int(r12)}\n'
    output += f'TV only: {int(r13)}\n'
    output += f'VD only: {int(r14)}\n'

    ng1 = (dfg[('syndromicmgmt.ng_ct', metric)][year] +
           dfg[('syndromicmgmt.ng_tv', metric)][year] +
           dfg[('syndromicmgmt.ng_vd', metric)][year] )
    ng2 = (dfg[('syndromicmgmt.ng_ct_tv', metric)][year] +
           dfg[('syndromicmgmt.ng_ct_vd', metric)][year] +
           dfg[('syndromicmgmt.ng_tv_vd', metric)][year] )
    sti4 = dfg[('syndromicmgmt.ng_ct_tv_vd', metric)][year]

    output += f'NG+1: {int(ng1)}\n'
    output += f'NG+2: {int(ng2)}\n'
    output += f'NG+3: {int(sti4)}\n'

    ct1 = (dfg['syndromicmgmt.ng_ct'].sum()[('syndromicmgmt.ng_ct', metric)][year] +
           dfg['syndromicmgmt.ct_tv'].sum()[('syndromicmgmt.ct_tv', metric)][year] +
           dfg['syndromicmgmt.ct_vd'].sum()[('syndromicmgmt.ct_vd', metric)][year] )
    ct2 = (dfg['syndromicmgmt.ng_ct_tv'].sum()[('syndromicmgmt.ng_ct_tv', metric)][year] +
           dfg['syndromicmgmt.ng_ct_vd'].sum()[('syndromicmgmt.ng_ct_vd', metric)][year] +
           dfg['syndromicmgmt.ct_tv_vd'].sum()[('syndromicmgmt.ct_tv_vd', metric)][year] )
    sti4 = dfg['syndromicmgmt.ng_ct_tv_vd'].sum()[('syndromicmgmt.ng_ct_tv_vd', metric)][year]

    output += f'CT+1: {int(ct1)}\n'
    output += f'CT+2: {int(ct2)}\n'
    output += f'CT+3: {int(sti4)}\n'

    tv1 = (dfg['syndromicmgmt.ng_tv'].sum()[('syndromicmgmt.ng_tv', metric)][year] +
           dfg['syndromicmgmt.ct_tv'].sum()[('syndromicmgmt.ct_tv', metric)][year] +
           dfg['syndromicmgmt.tv_vd'].sum()[('syndromicmgmt.tv_vd', metric)][year] )
    tv2 = (dfg['syndromicmgmt.ng_ct_tv'].sum()[('syndromicmgmt.ng_ct_tv', metric)][year] +
           dfg['syndromicmgmt.ng_tv_vd'].sum()[('syndromicmgmt.ng_tv_vd', metric)][year] +
           dfg['syndromicmgmt.ct_tv_vd'].sum()[('syndromicmgmt.ct_tv_vd', metric)][year] )
    sti4 = dfg['syndromicmgmt.ng_ct_tv_vd'].sum()[('syndromicmgmt.ng_ct_tv_vd', metric)][year]

    output += f'TV+1: {int(tv1)}\n'
    output += f'TV+2: {int(tv2)}\n'
    output += f'TV+3: {int(sti4)}\n'

    vd1 = (dfg['syndromicmgmt.ng_vd'].sum()[('syndromicmgmt.ng_vd', metric)][year] +
           dfg['syndromicmgmt.ct_vd'].sum()[('syndromicmgmt.ct_vd', metric)][year] +
           dfg['syndromicmgmt.tv_vd'].sum()[('syndromicmgmt.tv_vd', metric)][year] )
    vd2 = (dfg['syndromicmgmt.ng_ct_vd'].sum()[('syndromicmgmt.ng_ct_vd', metric)][year] +
           dfg['syndromicmgmt.ng_tv_vd'].sum()[('syndromicmgmt.ng_tv_vd', metric)][year] +
           dfg['syndromicmgmt.ct_tv_vd'].sum()[('syndromicmgmt.ct_tv_vd', metric)][year] )
    sti4 = dfg['syndromicmgmt.ng_ct_tv_vd'].sum()[('syndromicmgmt.ng_ct_tv_vd', metric)][year]

    output += f'VD+1: {int(vd1)}\n'
    output += f'VD+2: {int(vd2)}\n'
    output += f'VD+3: {int(sti4)}\n'

    # Number with 2 STIs
    sti2 = (dfg['syndromicmgmt.ng_ct'].sum()[('syndromicmgmt.ng_ct', metric)][year] +
            dfg['syndromicmgmt.ng_tv'].sum()[('syndromicmgmt.ng_tv', metric)][year] +
            dfg['syndromicmgmt.ng_vd'].sum()[('syndromicmgmt.ng_vd', metric)][year] +
            dfg['syndromicmgmt.ct_tv'].sum()[('syndromicmgmt.ct_tv', metric)][year] +
            dfg['syndromicmgmt.ct_vd'].sum()[('syndromicmgmt.ct_vd', metric)][year] +
            dfg['syndromicmgmt.tv_vd'].sum()[('syndromicmgmt.tv_vd', metric)][year] )

    # Number with 3 STIs
    sti3 = (dfg['syndromicmgmt.ng_ct_tv'].sum()[('syndromicmgmt.ng_ct_tv', metric)][year] +
            dfg['syndromicmgmt.ng_ct_vd'].sum()[('syndromicmgmt.ng_ct_vd', metric)][year] +
            dfg['syndromicmgmt.ng_tv_vd'].sum()[('syndromicmgmt.ng_tv_vd', metric)][year] +
            dfg['syndromicmgmt.ct_tv_vd'].sum()[('syndromicmgmt.ct_tv_vd', metric)][year] )

    output += f'Number with 1 STI: {int(r11+r12+r13+r14)}\n'
    output += f'Number with 2 STIs: {int(sti2)}\n'
    output += f'Number with 3 STIs: {int(sti3)}\n'

    print(output)

    return


def plot_tx(df, start_year=2000, percentile_pairs=None, title='tx_plots'):
    set_font(size=20)
    fig, axes = pl.subplots(2, 2, figsize=(8, 7))

    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    pn = 0
    # AMR
    ax = axes[pn]
    rname = 'ng.rel_treat'
    x = np.unique(dfplot['year'])
    y = dfplot.groupby(by='year')[rname].mean()[(rname, '50%')]
    line, = ax.plot(x[:-1], y[:-1] * 100)
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot.groupby(by='year')[rname].mean()[(rname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1]* 100, yu[:-1]* 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Susceptible to treatment (%)')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + ".png", dpi=100)

    return fig


if __name__ == '__main__':

    plot_single = False
    plot_multi = True

    if plot_multi:
        df_stats = sc.loadobj('results/multi_res_stats.df')
        percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
        plot_hiv_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs)
        plot_sti_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs, which='multi')

        # # Coinfection stats
        # dfplot = df_stats.iloc[(df_stats.index >= 2000) & (df_stats.index <= 2025)]
        # print_results(dfplot)

        # # Treatment results
        # plot_tx(df_stats, start_year=2000, percentile_pairs=percentile_pairs)

    if plot_single:
        df = sc.loadobj('results/sim.df')
        plot_sti_sims(df, start_year=2000, which='single')
