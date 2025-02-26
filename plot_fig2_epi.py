"""
Plot parameter estimates of care seeking alongside prevalence
Requires running model.py first with calibrated parameters.
Steps:
    1. Run run_calibration to generate the files 'results/zim_sti_calib_treat80.obj'
    2. Run model.py to generate the files 'results/epi_df.df' and 'results/treat80_sw.df'
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
    width = 0.6

    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'F', 'non_client': 'M'}
    colors = ['#d46e9c', '#2f734a', '#d46e9c', '#2f734a']
    alphas = [0.9, 0.9, 0.3, 0.3]

    # New infections acquired by sex and sex work status
    bottom = np.zeros(2)
    x = np.array([0.5, 1.5])
    g = 0
    for group, glabel in groups.items():
        vals = np.array([
            df[f'new_infections_{group}_{disease}'][10:30].mean(),
            df[f'new_transmissions_{group}_{disease}'][10:30].mean(),
        ])
        p = ax.barh(x, vals, width, label=glabel, left=bottom, color=colors[g], alpha=alphas[g])
        ax.bar_label(p, labels=[glabel, glabel], label_type='center')
        bottom += vals
        g += 1

    ax.set_title(disease.upper()+" infections\n2000-2020 average")
    sc.SIticks(ax, axis='x')
    ax.set_xlim(left=0)
    ax.set_yticks(x, ['Acquired', 'Transmitted'])
    return ax

# %% Run as a script
if __name__ == '__main__':

    scenario = 'treat80'
    epi_df = sc.loadobj('results/epi_df.df')
    sw_df = sc.loadobj(f'results/{scenario}_sw.df')

    # Initialize plot
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(20, 8))
    color = sc.vectocolor([.5, .8, 1])[1]
    axes = axes.ravel()
    scolors = ['#ee7989', '#4682b4']

    # Plot transmission and acquisition
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        ax = axes[ai]
        ax = plot_infections_by_sw(sw_df, disease=disease, ax=ax)
        ax.set_title(disease.upper()+' cases acquired and transmitted')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_ylim(bottom=0)

    # Plot prevalence by age
    for ai, disease in enumerate(['ng', 'ct', 'tv']):
        ax = axes[ai+3]
        thisdf = epi_df.loc[(epi_df.disease == disease) & (epi_df.age > 1)]
        thisdf['prevalence'] = thisdf['prevalence']*100
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
