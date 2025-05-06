"""
Plot degree of overtreatment by scenario
Generate files by running run_syndromic_scens.py on VMs -- takes about 10min
"""

# %% Imports and settings
import pandas as pd
import sciris as sc
import matplotlib.pyplot as pl
import seaborn as sns
import utils as ut


def plot_overtx(odf):

    # Plot settings
    ut.set_font(size=30)
    fig, axes = pl.subplots(1, 2, figsize=(12, 5))
    axes = axes.ravel()
    legendfont = 16

    clist = sc.gridcolors(3)
    colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

    # First row: reduction in overtreatment over time, by treatment type
    t = odf.timevec.unique()
    si = sc.findfirst(t, 2010)
    ei = sc.findfirst(t, 2040)

    for pn, (txname, txlabel) in enumerate(ut.tx_labels.items()):
        if txname != 'metronidazole':
            ax = axes[pn]
            for scenario in ut.scenarios:
                socdf = odf.loc[(odf.scenario == scenario) & (odf.treatment == txname) & (odf.variable == txname+'.new_treated_unnecessary_f')]
                socy = socdf['value'][si:ei]
                socy = socy.rolling(3, min_periods=1).mean()
                ax.plot(t[si:ei], socy, label=ut.txscenlabels[scenario], color=colors[scenario])
            for scenario in ut.scenarios:
                pocdf = odf.loc[(odf.scenario == (scenario+'poc')) & (odf.treatment == txname) & (odf.variable == txname+'.new_treated_unnecessary_f')]
                pocy = pocdf['value'][si:ei]
                pocy = pocy.rolling(3, min_periods=1).mean()
                ax.plot(t[si:ei], pocy, label=ut.txscenlabels[scenario], color=colors[scenario], ls='--')

            if pn == 0:
                h,l = ax.get_legend_handles_labels()  # #Get the legend handles and labels
                l1 = ax.legend(h[:3], l[:3], loc='upper left', frameon=False, prop={'size': legendfont})
                from matplotlib.lines import Line2D
                myHandle = [Line2D([], [], ls='-', color='k'), Line2D([], [], ls='--', color='k')]
                l2 = ax.legend(handles=myHandle, labels=['SOC', 'POC'], loc='lower left', bbox_to_anchor=(0, 0), frameon=False, prop={'size': legendfont})
                ax.add_artist(l1)

            ax.set_title(f'{txlabel} overtreatment')
            ax.set_ylim(bottom=0)
            sc.SIticks(ax)

    fig.tight_layout()
    pl.savefig(f"figures/fig4_slide_poctx.png", dpi=100)
    if show:
        pl.show()
    return


def plot_health(hdf):
    # Plot settings
    ut.set_font(size=20)
    fig, ax = pl.subplots(1, 1, figsize=(7, 4))
    legendfont = 14

    clist = sc.gridcolors(2)
    hdf.reset_index(inplace=True)
    toploth = hdf.loc[hdf.scenario != 'Treat-all']
    sns.boxplot(data=toploth, x="disease", y="infections", hue="scenario", palette=clist, ax=ax)
    ax.legend(frameon=False, prop={'size': legendfont})
    ax.set_title('% reduction in infections, 2027-2040')
    ax.set_ylim(0, 100)
    ax.set_xlabel('')
    ax.set_ylabel('')

    fig.tight_layout()
    pl.savefig(f"figures/fig4_impact_slide.png", dpi=100)
    if show:
        pl.show()
    return


def plot_fig4(odf, hdf, tdf, ddf):
    # Plot settings
    ut.set_font(size=30)
    fig = pl.figure(figsize=(30, 16))
    legendfont = 25

    gs1 = pl.GridSpec(1, 4, left=0.05, right=0.99, bottom=0.55, top=0.92, wspace=0.2)
    gs2 = pl.GridSpec(1, 3, left=0.05, right=0.99, bottom=0.05, top=0.45, wspace=0.2)

    clist = sc.gridcolors(3)
    colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

    # First row: reduction in overtreatment over time, by treatment type
    t = odf.timevec.unique()
    si = sc.findfirst(t, 2020)
    ei = sc.findfirst(t, 2040)

    # tx_dict = {'ng_tx': 'NG', 'ct_tx': 'CT'}  #ut.tx_labels

    # for pn, (txname, txlabel) in enumerate(tx_dict.items()):
    for pn, disease in enumerate(['ng', 'ct', 'tv']):
        ax = fig.add_subplot(gs1[pn])
        for scenario in ut.scenarios:
            socdf = odf.loc[(odf.scenario == scenario) & (odf.treatment == disease) & (odf.variable == disease+'.new_treated_unnecessary_f')]
            # socdf = odf.loc[(odf.scenario == scenario) & (odf.treatment == txname) & (odf.variable == txname+'.new_treated_unnecessary_f')]
            socy = socdf['value'][si:ei]
            socy = socy.rolling(3, min_periods=1).mean()
            ax.plot(t[si:ei], socy, label=ut.txscenlabels[scenario], color=colors[scenario])
        for scenario in ut.scenarios:
            pocdf = odf.loc[(odf.scenario == (scenario+'poc')) & (odf.treatment == disease) & (odf.variable == disease+'.new_treated_unnecessary_f')]
            # pocdf = odf.loc[(odf.scenario == (scenario+'poc')) & (odf.treatment == txname) & (odf.variable == txname+'.new_treated_unnecessary_f')]
            pocy = pocdf['value'][si:ei]
            pocy = pocy.rolling(3, min_periods=1).mean()
            ax.plot(t[si:ei], pocy, label=ut.txscenlabels[scenario], color=colors[scenario], ls='--')

        if pn == 0:
            h,l = ax.get_legend_handles_labels()  # #Get the legend handles and labels
            l1 = ax.legend(h[:3], l[:3], loc='upper left', frameon=False, prop={'size': legendfont})
            from matplotlib.lines import Line2D
            myHandle = [Line2D([], [], ls='-', color='k'), Line2D([], [], ls='--', color='k')]
            l2 = ax.legend(handles=myHandle, labels=['SOC', 'POC'], loc='lower left', bbox_to_anchor=(0, .1), frameon=False, prop={'size': legendfont})
            ax.add_artist(l1)

        ax.set_title(f'{disease.upper()} overtreatment')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)

    # First row summary plot: Cumulative reduction in overtreatment
    ax = fig.add_subplot(gs1[pn+1])
    tdf.reset_index(inplace=True)
    # tdf_plot = tdf.loc[tdf.treatment != 'MTNZ']
    tdf_plot = tdf
    sns.boxplot(data=tdf_plot, x="treatment", y="overtreatments", hue="scenario", palette=clist, ax=ax)
    ax.set_title('% reduction in\novertreatment 2027-2040')
    ax.set_ylim(0, 100)
    ax.get_legend().set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Second row, plot 1: Cumulative reduction in missed diagnoses
    pn = 0
    ax = fig.add_subplot(gs2[pn])
    hdf.reset_index(inplace=True)
    # sns.boxplot(data=hdf, x="disease", y="infections", hue="scenario", palette=clist, ax=ax)
    sns.boxplot(data=hdf, x="disease", y="new_false_neg", hue="scenario", palette=clist, ax=ax)
    ax.legend(frameon=False, prop={'size': legendfont}, loc='upper right')
    # ax.legend(frameon=False, prop={'size': legendfont})
    # ax.set_title('% reduction in infections, 2027-2040')
    ax.set_title('% reduction in\nundertreatment 2027-2040')
    # ax.set_ylim(-15, 50)
    ax.set_xlabel('')
    ax.set_ylabel('')
    pn += 1

    # Second row, plot 2: Cumulative reduction in infections
    ax = fig.add_subplot(gs2[pn])
    hdf.reset_index(inplace=True)
    sns.boxplot(data=hdf, x="disease", y="n_infected_f", hue="scenario", palette=clist, ax=ax)
    # ax.legend(frameon=False, prop={'size': legendfont})
    ax.get_legend().set_visible(False)
    # ax.set_title('% reduction in infections, 2027-2040')
    ax.set_title('% reduction in number infected')
    # ax.set_ylim(-15, 50)
    ax.set_xlabel('')
    ax.set_ylabel('')
    pn += 1

    # Reduction in time to treatment
    gs01 = gs2[pn].subgridspec(1, 3, wspace=0)

    for wn, disease in enumerate(['ng', 'ct', 'tv']):
        ax = fig.add_subplot(gs01[wn])
        toplot = ddf.loc[(ddf.disease == disease.upper())]
        sns.pointplot(data=toplot, x="poc", y="prop_fast", hue="scenario", palette=clist, ax=ax)

        # Labeling
        ax.get_legend().set_visible(False)
        if wn == 1: ax.set_title(f'% of cases resolved within 3m')
        if wn != 0: ax.set_yticklabels([])
        # if wn == 2: ax.legend(frameon=False, prop={'size': legendfont})
        ax.text(0.5, 95, disease.upper(), va="center", ha="center")

        ax.set_ylim(0, 100)
        ax.set_xlabel('')
        ax.set_ylabel('')
    pn += 1

    fig.tight_layout()

    pl.figtext(0.04, 0.95, 'A', fontsize=40, ha='center', va='center')
    pl.figtext(0.28, 0.95, 'B', fontsize=40, ha='center', va='center')
    pl.figtext(0.52, 0.95, 'C', fontsize=40, ha='center', va='center')
    pl.figtext(0.77, 0.95, 'D', fontsize=40, ha='center', va='center')
    pl.figtext(0.04, 0.48, 'E', fontsize=40, ha='center', va='center')
    pl.figtext(0.37, 0.48, 'F', fontsize=40, ha='center', va='center')
    pl.figtext(0.70, 0.48, 'G', fontsize=40, ha='center', va='center')

    pl.savefig(f"figures/fig4_poctx.png", dpi=100)
    if show:
        pl.show()
    return


if __name__ == '__main__':

    show = False

    # Load data
    odf = sc.loadobj(f'results/overtx.obj')
    hdf = sc.loadobj('results/synd_health.obj')
    tdf = sc.loadobj('results/synd_treat.obj')
    ddf = sc.loadobj('results/dur_df.obj')

    # Process durations
    process_durs = False
    if process_durs:
        mt = 4
        dfs = sc.autolist()
        idx = 0
        for rn, disease in enumerate(['ng', 'ct', 'tv']):
            dur_df = sc.loadobj(f'results/dur_df_{disease}.obj')
            for scenario, scenlabel in ut.scenlabels.items():
                for pocflag, poclabel in {'': 'SOC', 'poc': 'POC'}.items():
                    durs = dur_df.loc[dur_df.scenario.isin([scenario+pocflag])].groupby('months')['dur_inf'].mean()
                    res = durs.loc[durs.index <= mt].sum()
                    dd = dict(
                        scenario=scenlabel,
                        disease=disease.upper(),
                        poc=poclabel,
                        prop_fast=res*100,
                    )
                    dfs += pd.DataFrame(dd, index=[idx])
                    idx += 1
        ddf = pd.concat(dfs)
        sc.saveobj('results/dur_df.obj', ddf)

    # Condensed versions for slides
    make_main_fig = True
    if make_main_fig:
        plot_fig4(odf, hdf, tdf, ddf)

    # Make durations plot
    # ut.set_font(size=30)
    make_dur_plot = False
    if make_dur_plot:
        fig, axes = pl.subplots(3, 3, figsize=(15, 12))
        width = 0.45
        multiplier = 0
        for rn, disease in enumerate(['ng', 'ct', 'tv']):
            dur_df = sc.loadobj(f'results/dur_df_{disease}.obj')
            for cn, scenario in enumerate(['treat50', 'treat80', 'treat100']):
                ax = axes[rn, cn]
                for pocflag, poclabel in {'':'SOC', 'poc':'POC'}.items():
                    offset = width*multiplier-width/2
                    toplot = dur_df.loc[dur_df.scenario.isin([scenario+pocflag])].groupby('months')['dur_inf'].mean()
                    ax.bar(toplot.index.values+offset, toplot.values, width=width, label=poclabel)
                    multiplier += 1
                ax.set_title(f'{disease.upper()} duration of infection')
                if rn == 1 and cn == 0:
                    ax.legend()
        pl.show()

    # Condensed versions for slides
    make_slide_figs = False
    if make_slide_figs:
        # Plot overtreatment
        plot_overtx(odf)
        # Plot health impact
        plot_health(hdf)

    print('Done!')





