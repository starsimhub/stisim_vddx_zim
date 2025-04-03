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

if __name__ == '__main__':

    show = True

    # Load data
    odf = sc.loadobj(f'results/overtx.obj')
    hdf = sc.loadobj('results/synd_health.obj')
    tdf = sc.loadobj('results/synd_treat.obj')

    # Plot settings
    ut.set_font(size=30)
    fig = pl.figure(figsize=(20, 16))
    legendfont = 25

    gs1 = pl.GridSpec(1, 3, left=0.05, right=0.95, bottom=0.55, top=0.95, wspace=0.3)
    gs2 = pl.GridSpec(1, 2, left=0.05, right=0.95, bottom=0.05, top=0.45, wspace=0.1)

    clist = sc.gridcolors(3)
    colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

    # First row: reduction in overtreatment over time, by treatment type
    t = odf.timevec.unique()
    si = sc.findfirst(t, 2010)
    ei = sc.findfirst(t, 2040)

    for pn, (txname, txlabel) in enumerate(ut.tx_labels.items()):
        ax = fig.add_subplot(gs1[pn])
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
            l2 = ax.legend(handles=myHandle, labels=['SOC', 'POC'], loc='upper left', bbox_to_anchor=(0, .7), frameon=False, prop={'size': legendfont})
            ax.add_artist(l1)

        ax.set_title(f'{txlabel} overtreatment')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)

    # Second row, plot 1: Cumulative reduction in overtreatment
    ax = fig.add_subplot(gs2[0])
    tdf.reset_index(inplace=True)
    sns.boxplot(data=tdf, x="treatment", y="overtreatments", hue="scenario", palette=clist, ax=ax)
    ax.set_title('% reduction in overtreatment, 2027-2040')
    ax.set_ylim(0, 100)
    ax.get_legend().set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Second row, plot 2: Cumulative reduction in infections
    ax = fig.add_subplot(gs2[1])
    hdf.reset_index(inplace=True)
    sns.boxplot(data=hdf, x="disease", y="infections", hue="scenario", palette=clist, ax=ax)
    ax.legend(frameon=False, prop={'size': legendfont})
    ax.set_title('% reduction in infections, 2027-2040')
    ax.set_ylim(-20, 100)
    ax.set_xlabel('')
    ax.set_ylabel('')

    fig.tight_layout()
    pl.savefig(f"figures/fig4_poctx.png", dpi=100)
    if show:
        pl.show()

    print('Done!')





