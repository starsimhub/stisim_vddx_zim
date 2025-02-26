"""
Plot parameter estimates of care seeking alongside prevalence 
"""

# Import packages
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
import seaborn as sns
from utils import set_font


# %% Plotting functions
def plot_infections_by_sw(df, disease=None, ax=None, start_year=2000, end_year=2019):
    set_font(size=20)
    years = np.arange(start_year, end_year + 1)
    width = 0.8

    df['year'] = np.floor(np.round(df.index, 1)).astype(int)
    dfplot = df.iloc[(df.index >= (end_year - 1)) & (df.index <= end_year)]

    to_plot = [
        'sw_stats.new_infections_fsw',
        'sw_stats.new_infections_client',
        'sw_stats.new_infections_non_fsw',
        'sw_stats.new_infections_non_client',
        'sw_stats.new_transmissions_fsw',
        'sw_stats.new_transmissions_client',
        'sw_stats.new_transmissions_non_fsw',
        'sw_stats.new_transmissions_non_client',
    ]
    ysums = dfplot[to_plot].sum()
    # ysums = dfplot.groupby()[to_plot].sum()
    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'Other F', 'non_client': 'Other M'}
    colors = ['#d46e9c', '#2f734a', '#d46e9c', '#2f734a']
    alphas = [0.9, 0.9, 0.3, 0.3]

    # New infections acquired by sex and sex work status
    # bottom = np.zeros((n_years))
    bottom = np.zeros((2))
    x = np.arange(2)
    g = 0
    for group, glabel in groups.items():
        # vals = ysums[f'syph_sw_stats.new_infections_{group}'].values
        vals = np.array([ysums[f'syph_sw_stats.new_infections_{group}'], ysums[f'syph_sw_stats.new_transmissions_{group}']])
        p = ax.bar(x, vals, width, label=glabel, bottom=bottom, color=colors[g], alpha=alphas[g])
        ax.bar_label(p, labels=[glabel, glabel], label_type='center')
        bottom += vals
        g += 1

    ax.set_title("Syphilis infections\n2000-2020 average")
    # ax.legend(frameon=False, loc='upper right')
    sc.SIticks(ax)
    ax.set_ylim([0, 60_000])
    ax.set_xticks([0, 1], ['Acquired', 'Transmitted'])

    sc.figlayout()
    sc.savefig("figures/infections_by_sw_2.png", dpi=100)



# %% Run as a script
if __name__ == '__main__':

    scenario = 'treat80'
    epi_df = sc.loadobj('results/epi_df.df')

    # Initialize plot
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(20, 8))
    color = sc.vectocolor([.5, .8, 1])[1]
    axes = axes.ravel()
    scolors = ['#ee7989', '#4682b4']

    # Plot transmission and acquisition
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        ax = axes[ai]
        ax.set_title(disease.upper()+' cases acquired and transmitted')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(bottom=0)

    # Plot prevalence by age
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        ax = axes[ai+3]
        thisdf = epi_df.loc[(epi_df.disease == disease) & (epi_df.age > 1)]
        sns.barplot(data=thisdf, x="age", y="prevalence", hue="sex", ax=ax, palette=scolors)
        ax.set_title(disease.upper())
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.get_legend().set_visible(False)
        ax.set_title(disease.upper()+' prevalence by age')
        ax.set_ylim(bottom=0)

    fig.tight_layout()
    pl.savefig(f"figures/fig2_epi.png", dpi=100)

    print('Done.')
