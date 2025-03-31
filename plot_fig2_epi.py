"""
Plot parameter estimates of care seeking alongside prevalence
Requires running model.py first with calibrated parameters.
Steps:
    1. Run run_calibration to generate the files 'results/zim_sti_calib_treat80.obj'
    2. Run run_plot_data.py to generate the epi result files:
            epi_df = 'results/epi_df_{scenario}.df'
            sw_df = 'results/sw_df_{scenario}.df'
            hiv_df = 'results/hiv_df_{scenario}.df'

"""

# Import packages
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
import seaborn as sns
from utils import set_font
from matplotlib.gridspec import GridSpec


# %% Plotting functions
def plot_infections_by_sw(df, disease=None, ax=None, start_year=2000, end_year=2019):
    set_font(size=20)
    width = 0.6

    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'F', 'non_client': 'M'}
    colors = ['#d46e9c', '#2f734a', '#d46e9c', '#2f734a']
    alphas = [0.9, 0.9, 0.3, 0.3]

    si = sc.findfirst(df.index, start_year)
    ei = sc.findfirst(df.index, end_year)

    # New infections acquired by sex and sex work status
    bottom = np.zeros(2)
    x = np.array([0.5, 1.5])
    g = 0
    for group, glabel in groups.items():
        vals = np.array([
            df[f'new_infections_{group}_{disease}'][si:ei].mean(),
            df[f'new_transmissions_{group}_{disease}'][si:ei].mean(),
        ])
        p = ax.barh(x, vals, width, label=glabel, left=bottom, color=colors[g], alpha=alphas[g])
        ax.bar_label(p, labels=[glabel, glabel], label_type='center')
        bottom += vals
        g += 1

    ax.set_title(disease.upper()+f" infections\n{start_year}-{end_year} average")
    sc.SIticks(ax, axis='x')
    ax.set_xlim(left=0)
    ax.set_yticks(x, ['Acquired', 'Transmitted'])
    ax.set_xlabel('')
    ax.set_ylabel('')

    total_trans_ng = sw_df[f'new_transmissions_fsw_ng'][si:ei].mean()+sw_df[f'new_transmissions_client_ng'][10:30].mean()+sw_df[f'new_transmissions_non_fsw_ng'][10:30].mean()+sw_df[f'new_transmissions_non_client_ng'][10:30].mean()
    total_trans_ct = sw_df[f'new_transmissions_fsw_ct'][si:ei].mean()+sw_df[f'new_transmissions_client_ct'][10:30].mean()+sw_df[f'new_transmissions_non_fsw_ct'][10:30].mean()+sw_df[f'new_transmissions_non_client_ct'][10:30].mean()
    total_trans_tv = sw_df[f'new_transmissions_fsw_tv'][si:ei].mean()+sw_df[f'new_transmissions_client_tv'][10:30].mean()+sw_df[f'new_transmissions_non_fsw_tv'][10:30].mean()+sw_df[f'new_transmissions_non_client_tv'][10:30].mean()
    print(f'NG SW share: {(sw_df[f"new_transmissions_fsw_ng"][si:ei].mean()+sw_df[f"new_transmissions_client_ng"][10:30].mean())/total_trans_ng}')
    print(f'CT SW share: {(sw_df[f"new_transmissions_fsw_ct"][si:ei].mean()+sw_df[f"new_transmissions_client_ct"][10:30].mean())/total_trans_ct}')
    print(f'TV SW share: {(sw_df[f"new_transmissions_fsw_tv"][si:ei].mean()+sw_df[f"new_transmissions_client_tv"][10:30].mean())/total_trans_tv}')
    # ax.set_ylim(bottom=0)
    return ax


def plot_hiv(hiv_df, ax=None):

    set_font(size=20)
    colors = ['#ee7989', '#ee7989', '#4682b4', '#4682b4']
    linestyles = ['--', '-', '--', '-']
    rdict = {'symp_prev_no_hiv_f': 'HIV- F', 'symp_prev_has_hiv_f': 'HIV+ F', 'symp_prev_no_hiv_m': 'HIV- M', 'symp_prev_has_hiv_m': 'HIV+ M'}

    cn = 0
    bi = 20
    for rname, rlabel in rdict.items():
        x = hiv_df.index[bi:]
        y = hiv_df[rname][bi:]
        # y = y.rolling(10, min_periods=1).mean()
        ax.plot(x, y*100, color=colors[cn], ls=linestyles[cn], label=rlabel)
        cn += 1

    ax.legend()
    ax.set_ylabel('')
    ax.set_xlabel('')
    ax.set_title('Prevalence of discharge (%)')
    ax.set_ylim(bottom=0)

    return ax


# %% Run as a script
if __name__ == '__main__':

    scenario = 'treat80'
    epi_df = sc.loadobj(f'results/epi_df_{scenario}.df')
    sw_df = sc.loadobj(f'results/sw_df_{scenario}.df')
    hiv_df = sc.loadobj(f'results/hiv_df_{scenario}.df')

    # Initialize plot
    set_font(size=20)
    # fig = pl.figure(figsize=(20, 12))
    # gs1 = GridSpec(3, 3, left=0.05, right=0.98, wspace=0.05, hspace=0.05)
    fig, axes = pl.subplots(2, 3, figsize=(20, 8))
    color = sc.vectocolor([.5, .8, 1])[1]
    axes = axes.ravel()
    scolors = ['#ee7989', '#4682b4']

    # Plot transmission and acquisition
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        # ax = fig.add_subplot(gs1[0, ai])
        ax = axes[ai]
        ax = plot_infections_by_sw(sw_df, disease=disease, ax=ax)

    # Plot prevalence by age
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        # ax = fig.add_subplot(gs1[1, ai])
        ax = axes[ai+3]
        thisdf = epi_df.loc[(epi_df.disease == disease) & (epi_df.age != '0-15') & (epi_df.age != '65+')]
        thisdf['prevalence'] = thisdf['prevalence']*100
        sns.barplot(data=thisdf, x="age", y="prevalence", hue="sex", ax=ax, palette=scolors)
        ax.set_title(disease.upper())
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.get_legend().set_visible(False)
        ax.set_title(disease.upper()+' prevalence by age')
        ax.set_ylim(bottom=0)

    # # Plot symptoms by sex and HIV status
    # ax = fig.add_subplot(gs1[2, :])
    # ax = plot_hiv(hiv_df, ax=ax)

    fig.tight_layout()
    pl.savefig(f"figures/fig2_epi.png", dpi=100)

    print('Done.')
