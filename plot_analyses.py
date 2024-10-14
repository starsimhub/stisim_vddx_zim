# %% Imports and settings
import sciris as sc
import pylab as pl
import seaborn as sns
import pandas as pd

# From this repo
from utils import set_font, get_y

scenarios = dict(soc='SOC', panel='Panel')


def plot_tx_diffs(df):
    set_font(size=20)
    fig, axes = pl.subplots(4, 4, figsize=(25, 15))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}  #, 'bv': 'BV'}
    test_res = {'new_false_pos': ' - false positives'}
    # test_res = {'new_true_pos': ' - true positives', 'new_false_pos': ' - false positives', 'new_false_neg': ' - false negatives', 'new_true_neg': ' - true negatives'}

    for trname, trlabel in test_res.items():
        for dname, dlabel in labels.items():
            ax = axes[pn]
            sns.lineplot(df, x=df.index, y=f"{dname}.{trname}", hue="scenario", ax=ax)
            ax.set_title(dlabel + trlabel)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sc.SIticks(ax=ax)
            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/sti_tx_analyses.png", dpi=100)
    return


def plot_result1(df, which='multi'):
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    test_res = {'new_false_pos': ' - false positives'}
    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}  #, 'bv': 'BV'}
    # test_res = {'new_true_pos': ' - true positives', 'new_false_pos': ' - false positives', 'new_false_neg': ' - false negatives', 'new_true_neg': ' - true negatives'}

    for trname, trlabel in test_res.items():
        for dname, dlabel in labels.items():
            resname = f"{dname}.{trname}"
            ax = axes[pn]
            for sname, slabel in scenarios.items():
                thisdf = df.loc[[sname]]

                x = thisdf.index.get_level_values('year')
                y = get_y(thisdf, which, resname)
                ax.plot(x, y, label=slabel)
                if which == 'multi':
                    quartiles1 = thisdf[(resname, '10%')]
                    quartiles3 = thisdf[(resname, '90%')]
                    ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)

            if pn == 0: ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
            ax.set_title(dlabel + trlabel)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sc.SIticks(ax=ax)
            ax.axvline(x=intv_year, color='k', ls='--')
            pn += 1

    # Treatment stats by disease
    panel_df = df.loc[['panel']]
    for tname, tlabel in {'ng_tx': 'Ceftriaxone', 'ct_tx': 'Doxycycline', 'metronidazole':'Metronidazole'}.items():  #, 'bv']:
        ax = axes[pn]
        x = panel_df.index.get_level_values('year')
        Y = [
            panel_df[(tname+'.new_treated_success', '50%')],
            panel_df[(tname+'.new_treated_failure', '50%')],
            panel_df[(tname+'.new_treated_unnecessary', '50%')],
        ]
        labels = ["Treatment success", "Treatment failure", "Overtreatment"]
        ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.gridcolors(3))
        ax.set_title(tlabel)
        if pn == 5: ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
        ax.set_ylim(bottom=0)
        ax.axvline(x=intv_year, color='k', ls='--')
        sc.SIticks(ax=ax)
        pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_false_pos.png", dpi=100)
    return


def plot_result2(df, which='multi'):
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    diseases = {'ng': "NG", 'ct': "CT", 'tv': "TV"}
    results = {'new_symptomatic': ' - new symptomatic infections', 'n_symptomatic': ' - symptomatic burden', }  #'new_true_pos_f': ' - true positives (F)'}

    for rname, rlabel in results.items():
        for dname, dlabel in diseases.items():
            ax = axes[pn]
            resname = f"{dname}.{rname}"

            for sname, slabel in scenarios.items():
                thisdf = df.loc[[sname]]

                x = thisdf.index.get_level_values('year')
                y = get_y(thisdf, which, resname)
                ax.plot(x, y, label=slabel)

                if which == 'multi':
                    quartiles1 = thisdf[(resname, '10%')]
                    quartiles3 = thisdf[(resname, '90%')]
                    ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)

            ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title(dlabel + rlabel)
            sc.SIticks(ax=ax)
            if pn == 0: ax.legend(frameon=False, prop={'size': legend_font}, loc='upper left')
            # else: ax.legend_.remove()
            ax.axvline(x=intv_year, color='k', ls='--')
            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_epi.png", dpi=100)
    return


def plot_result3(df):
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    test_res = {'new_false_pos': ' - false positives averted', 'n_infected': ' - incremental infections'}
    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}  #, 'bv': 'BV'}
    bv_dict = {0.12: '15%', 0.15: '20%', 0.2: '25%'}
    colors = ['#E2D6BF', '#B6985E', '#5F4D2B']  #sc.vectocolor(4, cmap='cividis')
    # test_res = {'new_true_pos': ' - true positives', 'new_false_pos': ' - false positives', 'new_false_neg': ' - false negatives', 'new_true_neg': ' - true negatives'}

    for trname, trlabel in test_res.items():
        for dname, dlabel in labels.items():
            ax = axes[pn]
            cn = 0
            for bv_val, bv_label in bv_dict.items():
                df_base = df[df.index.get_level_values('bv_beta') == bv_val].loc[['soc']]
                df_intv = df[df.index.get_level_values('bv_beta') == bv_val].loc[['panel']]
                x = df_base.index.get_level_values('year')
                y = df_base[(f"{dname}.{trname}", '50%')].values - df_intv[(f"{dname}.{trname}", '50%')].values
                ax.plot(x, y, label=bv_label, color=colors[cn])

                quartiles1 = df_base[(f"{dname}.{trname}", '10%')].values - df_intv[(f"{dname}.{trname}", '10%')].values
                quartiles3 = df_base[(f"{dname}.{trname}", '90%')].values - df_intv[(f"{dname}.{trname}", '90%')].values
                ax.fill_between(x, quartiles1, quartiles3, alpha=0.3, color=colors[cn])

                cn += 1

            ax.set_title(dlabel + trlabel)
            # ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sc.SIticks(ax=ax)
            if pn == 2: ax.legend(frameon=False, prop={'size': legend_font}, loc='upper left')
            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_sens.png", dpi=100)
    return


def plot_result4(df):
    """ Comparison of the direct model and the transmission model """
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(1, 2, figsize=(10, 5))
    axes = axes.ravel()
    pn = 0
    bi = 20

    resdict = {'new_female_symptomatic': ' - cases averted in women'}  #,'new_symptomatic': ' - incremental true positives'}  #,
    ddict = {'ng': "NG", 'ct': "CT"}  #, 'tv': "TV"}
    colors = sc.gridcolors(2)
    dm = pd.read_csv('data/new_tp_zim.csv')
    bv_val = 0.15

    for rname, rlabel in resdict.items():
        for dname, dlabel in ddict.items():
            ax = axes[pn]
            df_base = df[df.index.get_level_values('bv_beta') == bv_val].loc[['soc']]
            df_intv = df[df.index.get_level_values('bv_beta') == bv_val].loc[['panel']]

            # df_base = df.loc[['soc']]
            # df_intv = df.loc[['panel']]
            x = df_base.index.get_level_values('year')[bi:]
            y = (df_base[(f"{dname}.{rname}", '50%')].values - df_intv[(f"{dname}.{rname}", '50%')].values)[bi:]
            ax.plot(x, y, color=colors[0], label='Indirect impact')
            quartiles1 = (df_base[(f"{dname}.{rname}", '10%')].values - df_intv[(f"{dname}.{rname}", '10%')].values)[bi:]
            quartiles3 = (df_base[(f"{dname}.{rname}", '90%')].values - df_intv[(f"{dname}.{rname}", '90%')].values)[bi:]
            ax.fill_between(x, quartiles1, quartiles3, alpha=0.3, color=colors[0])

            # Add direct model
            xdm = dm.year
            ydm = dm[dname+'.new_true_pos']
            ax.plot(xdm, ydm, color=colors[1], label='Direct impact')

            ax.set_title(dlabel + rlabel)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sc.SIticks(ax=ax)
            ax.set_ylim(bottom=0)
            # if pn == 2: ax.legend(frameon=False, prop={'size': legend_font}, loc='upper left')
            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_comparison.png", dpi=100)
    return


def plot_health(df):

    # Health impact
    set_font(size=20)
    fig, axes = pl.subplots(2, 2, figsize=(12, 10))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0
    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV", 'hiv': 'HIV'}
    # labels = {'ng_tx': "Ceftriaxone", 'ct_tx': "Doxycycline", 'metronidazole': "Metronidazole"}

    for tx, txlabel in labels.items():
        ax = axes[pn]
        sns.lineplot(df, x=df.index, y=tx+".new_infections", hue="scenario", ax=ax)
        # sns.lineplot(df, x=df.index, y=f"{tx}.n_infected", hue="scenario", ax=ax)
        ax.set_title(txlabel)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(txlabel + ' - new infections')
        sc.SIticks(ax=ax)
        pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/sti_tx_analyses_health.png", dpi=100)

    return


if __name__ == '__main__':

    make_stats = True

    # Make stats
    if make_stats:
        df = sc.loadobj('results/scens.obj')
        df = df.iloc[(df.index >= 2000)]
        centraldf = df.loc[df.bv_beta == 0.15]

        df['year'] = df.index
        df_stats = df.groupby(by=['scenario', 'year', 'bv_beta']).describe(percentiles=[.1, .5, .9])
        central_df_stats = df_stats[df_stats.index.get_level_values('bv_beta') == 0.15]
        sc.saveobj('results/scen_stats_central.obj', central_df_stats)
        sc.saveobj('results/scen_stats.obj', df_stats)

    central_df_stats = sc.loadobj('results/scen_stats_central.obj')
    df_stats = sc.loadobj('results/scen_stats.obj')

    # Plots
    plot_result1(central_df_stats)
    plot_result2(central_df_stats)
    plot_result3(df_stats)
    plot_result4(df_stats)


