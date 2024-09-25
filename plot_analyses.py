# %% Imports and settings
import sciris as sc
import pylab as pl
import seaborn as sns

# From this repo
from utils import set_font


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


def plot_result1(df):
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
            ax = axes[pn]
            sns.lineplot(df, x=df.index, y=f"{dname}.{trname}", hue="scenario", ax=ax)
            ax.set_title(dlabel + trlabel)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sc.SIticks(ax=ax)
            pn += 1

    # Treatment stats by disease
    panel_df = df.loc[df.scenario == 'panel']
    for tname, tlabel in {'ng_tx': 'Ceftriaxone', 'ct_tx': 'Doxycycline', 'metronidazole':'Metronidazole'}.items():  #, 'bv']:
        ax = axes[pn]
        x = panel_df.index
        Y = [
            panel_df[tname+'.new_treated_success'],
            panel_df[tname+'.new_treated_failure'],
            panel_df[tname+'.new_treated_unnecessary'],
        ]
        labels = ["Treatment success", "Treatment failure", "Overtreatment"]
        ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.gridcolors(3))
        ax.set_title(tlabel)
        if pn == 5: ax.legend(frameon=True, prop={'size': legend_font}, loc='lower left')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_false_pos.png", dpi=100)
    return


def plot_result2(df):
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}
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
        if pn == 0: ax.legend(frameon=False, prop={'size': legend_font}, loc='upper left')
        else: ax.legend_.remove()
        ax.axvline(x=intv_year, color='k', ls='--')
        pn += 1

    test_res = {'new_true_pos': ' - true positives (F)'}
    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}  #, 'bv': 'BV'}
    # test_res = {'new_true_pos': ' - true positives', 'new_false_pos': ' - false positives', 'new_false_neg': ' - false negatives', 'new_true_neg': ' - true negatives'}

    for trname, trlabel in test_res.items():
        for dname, dlabel in labels.items():
            ax = axes[pn]
            sns.lineplot(df, x=df.index, y=f"{dname}.{trname}_f", hue="scenario", ax=ax)
            ax.set_title(dlabel + trlabel)
            ax.set_ylim(bottom=0)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.axvline(x=intv_year, color='k', ls='--')
            ax.legend_.remove()
            sc.SIticks(ax=ax)
            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/analyses_epi.png", dpi=100)
    return


def plot_result3(df):
    set_font(size=20)
    legend_font = 18
    fig, axes = pl.subplots(1, 3, figsize=(15, 4))
    axes = axes.ravel()
    intv_year = 2027
    pn = 0

    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}
    # labels = {'ng_tx': "Ceftriaxone", 'ct_tx': "Doxycycline", 'metronidazole': "Metronidazole"}

    test_res = {'new_false_pos': ' - false positives averted'}  #, 'new_infections': ' - incremental infections'}
    labels = {'ng': "NG", 'ct': "CT", 'tv': "TV"}  #, 'bv': 'BV'}
    bv_dict = {0.12: '15%', 0.15: '20%', 0.2: '25%'}
    colors = ['#E2D6BF', '#B6985E', '#5F4D2B']  #sc.vectocolor(4, cmap='cividis')
    # test_res = {'new_true_pos': ' - true positives', 'new_false_pos': ' - false positives', 'new_false_neg': ' - false negatives', 'new_true_neg': ' - true negatives'}

    for trname, trlabel in test_res.items():
        for dname, dlabel in labels.items():
            ax = axes[pn]
            cn = 0
            for bv_val, bv_label in bv_dict.items():
                df_base = df.loc[(df.scenario == 'soc') & (df.bv_beta == bv_val)]
                df_intv = df.loc[(df.scenario == 'panel') & (df.bv_beta == bv_val)]
                ax.plot(df_base.index, df_base[f"{dname}.{trname}"] - df_intv[f"{dname}.{trname}"], label=bv_label, color=colors[cn])
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

    df = sc.loadobj('results/scens.obj')
    df = df.iloc[(df.index >= 2000)]
    centraldf = df.loc[df.bv_beta == 0.15]

    # Plots
    # plot_tx_diffs(df)
    # plot_health(df)
    plot_result1(centraldf)
    plot_result2(centraldf)
    plot_result3(df)


