# %% Imports and settings
import numpy as np
import pandas as pd
import sciris as sc
import matplotlib.pyplot as pl
from utils import set_font, get_y

location = 'zimbabwe'
show = False # Whether to show the plots by default (else just save)


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots', show=show):
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
    resname = 'hiv_prevalence'
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
    if show:
        pl.show()

    return fig


def plot_sti_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='sti_plots', fext='', show=show):
    """ Create quantile or individual sim plots of STIs """
    set_font(size=30)
    fig, axes = pl.subplots(3, 3, figsize=(25, 12))
    axes = axes.ravel()
    if which == 'multi': alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    sti_data = pd.read_csv(f'data/{location}_sti_data.csv')
    sti_data = sti_data.loc[(sti_data.time >= start_year) & (sti_data.time <= end_year)]
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]

    disease_map = {'ng': 'NG', 'ct': 'CT', 'tv': 'TV'}  #, 'bv': 'Other'}
    result_map = {
        # 'prevalence_f_15_25': 'Prevalence F 15-25',
        # 'n_infected_f_15_25': 'Burden F 15-25',
        'prevalence': 'Prevalence',
        'new_infections': 'Infections',
        'n_infected': 'Burden',
    }

    pn = 0
    x = dfplot.index

    # Incidence
    for dname, dlabel in disease_map.items():
        for rname, reslabel in result_map.items():
            ax = axes[pn]

            resname = dname+'_'+rname
            if resname in sti_data.columns:
                ax.scatter(sti_data.time, sti_data[resname], color='k', label='Data')
            y = get_y(dfplot, which, resname)
            line, = ax.plot(x, y, label=reslabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())

            ax.set_title(dlabel+' '+reslabel)
            if pn == 2: ax.legend(frameon=False, prop={'size': 20})
            ax.set_ylim(bottom=0)
            sc.SIticks(ax=ax)

            pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + "_" + fext + ".png", dpi=100)
    if show:
        pl.show()

    return fig


def plot_sti_tx(df, start_year=2000, end_year=2025, fext='', sex=None, show=show):
    set_font(size=24)
    legend_font = 20
    fig, axes = pl.subplots(2, 3, figsize=(20, 8))
    axes = axes.ravel()
    dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
    dfplot = dfplot.set_index('timevec')
    if sex is not None:
        sex = '_' + sex
        sexlabel = ' - ' + sex.upper()
    else:
        sex = ''
        sexlabel = ''

    pn = 0

    # Care seeker epidemiology
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['syndromicmgmt_new_sti1'+sex],
        dfplot['syndromicmgmt_new_sti2'+sex],
        dfplot['syndromicmgmt_new_sti3'+sex],
        dfplot['syndromicmgmt_new_sti4'+sex],
    ]
    labels = ["1 infection", "2 infections", "3 infections", "4 infections"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.vectocolor(4, reverse=True))
    ax.set_title('Coinfection among care seekers'+sexlabel)
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Care seeker management
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['syndromicmgmt_new_tx0'+sex],
        dfplot['syndromicmgmt_new_tx1'+sex],
        dfplot['syndromicmgmt_new_tx2'+sex],
        dfplot['syndromicmgmt_new_tx3'+sex],
    ]
    labels = ["0", "1", "2", "3"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.vectocolor(4, reverse=True))
    ax.set_title('Treatments prescribed'+sexlabel)
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Overtreatment
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['ng_tx_new_treated_unnecessary'+sex],
        dfplot['ct_tx_new_treated_unnecessary'+sex],
        dfplot['metronidazole_new_treated_unnecessary'+sex],
    ]
    labels = ["Ceftriaxone", "Doxycycline", "Metronidazole"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.gridcolors(4))
    ax.set_title('Treatments overprescribed')
    ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Treatment stats by disease
    for disease in ['ng', 'ct', 'tv']:  #, 'bv']:
        ax = axes[pn]
        x = dfplot.index
        Y = [
            dfplot[disease+'_new_treated_success'+sex],
            dfplot[disease+'_new_treated_failure'+sex],
            dfplot[disease+'_new_treated_unnecessary'+sex],
        ]
        labels = ["Treatment success", "Treatment failure", "Overtreatment"]
        ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.gridcolors(4))
        ax.set_title('Care seekers treated for ' + disease.upper())
        if pn == 5: ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    sc.figlayout()
    sc.savefig(f"figures/sti_tx{fext}.png", dpi=100)
    if show:
        pl.show()

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


def plot_ng_sim(df, start_year=2000, end_year=2025, which='single', title='ng_plots', show=show):
    """ Create quantile or individual sim plots of NG """
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()

    ng_data = pd.read_csv(f'data/zimbabwe_ng_data.csv')
    ng_data = ng_data.loc[(ng_data.year >= start_year) & (ng_data.year <= end_year)]
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    pn = 0
    dname = 'ng'

    # Incidence
    ax = axes[pn]
    resname = 'ng.new_infections'
    ax.scatter(ng_data.year, ng_data[resname], label='Data', color='k')
    resnames = {'Total': dname+'.new_infections', 'Symptomatic': dname+'.new_symptomatic'}  #, 'Care seekers': dname+'.new_care_seekers'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = dfplot[rname]
        ax.plot(x, y, label=rlabel)
        ax.set_title('Infections')
        ax.legend(frameon=False, prop={'size': 20})
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
    pn += 1

    # Burden
    ax = axes[pn]
    resname = dname+'.n_infected'
    ax.scatter(ng_data.year, ng_data[resname], label='Data', color='k')
    resnames = {'Total': dname+'.n_infected', 'Symptomatic': dname+'.n_symptomatic'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = dfplot[rname]
        ax.plot(x, y, label=rlabel)
        ax.set_title('Burden')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
    pn += 1

    # Prevalence
    ax = axes[pn]
    resnames = {'Total': dname+'.adult_prevalence', 'Symptomatic': dname+'.symp_adult_prevalence'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = dfplot[rname]
        ax.plot(x, y*100, label=rlabel)
        ax.set_title('Prevalence (%)')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
    pn += 1

    # Testing outcomes
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['ng.new_true_pos'],
        dfplot['ng.new_false_pos'],
        dfplot['ng.new_true_neg'],
        dfplot['ng.new_false_neg'],
    ]
    labels = ["True positive", "False positives", "True negatives", "False negatives"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels)
    ax.set_title('Testing outcomes')
    ax.legend(frameon=True, prop={'size': 14}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Treatment stats by disease
    ax = axes[pn]
    x = dfplot.index
    Y = [
        dfplot['ng.new_treated_unnecessary'],
        dfplot['ng.new_treated_success'],
        dfplot['ng.new_treated_failure'],
    ]
    labels = ["Overtreatment", "Success", "Failure"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels)
    ax.set_title('Care seekers treated for NG')
    ax.legend(frameon=True, prop={'size': 14}, loc='lower left')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    ax = axes[pn]
    rname = 'ng.rel_treat'
    x = dfplot.index
    y = dfplot[rname]
    ax.plot(x, y)
    ax.set_title('Susceptible to treatment (%)')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)
    if show:
        pl.show()

    return fig


if __name__ == '__main__':

    show = True
    plot_single = False
    plot_multi = True

    if plot_multi:
        df_stats = sc.loadobj('results/multi_res_stats.df')
        percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]
        plot_hiv_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs, show=show)
        plot_sti_sims(df_stats, start_year=2000, percentile_pairs=percentile_pairs, which='multi', show=show)

    if plot_single:
        df = sc.loadobj('results/sim.df')
        plot_sti_sims(df, start_year=2000, which='single', show=show)
