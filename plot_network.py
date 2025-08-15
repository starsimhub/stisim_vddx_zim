"""
Plot network
"""

# Import packages
import sciris as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font


# %% Run as a script
if __name__ == '__main__':

    sim = sc.loadobj('results/treat80.sim')

    # Initialize plot
    set_font(size=25)
    fig = pl.figure(figsize=(24, 10))
    scolors = ['#ee7989', '#4682b4']

    gs1 = pl.GridSpec(2, 1, left=0.05, right=0.30, bottom=0.1, top=0.91, hspace=0.5)
    gs2 = pl.GridSpec(2, 1, left=0.35, right=0.60, bottom=0.1, top=0.91, hspace=0.5)
    gs3 = pl.GridSpec(1, 1, left=0.65, right=0.99, bottom=0.1, top=0.91)

    # Debut age
    a = sim.analyzers.debutage
    data = dict(
        f=dict(bins = [15, 18, 20, 22, 25], props = [0.057, 0.401, 0.66, 0.819, 0.929]),
        m=dict(bins = [15, 18, 20, 22, 25], props = [0.044, 0.244, 0.441, 0.639, 0.816])
    )
    # Females
    ax = fig.add_subplot(gs1[0])
    for row in a.prop_active_f:
        ax.plot(a.bins, row, color=scolors[0], alpha=0.5)
    ax.scatter(data['f']['bins'], data['f']['props'], color='k')
    ax.set_xlabel('Age')
    ax.set_ylabel('Share')
    ax.set_title('Proportion of females\nwho are sexually active')

    ax = fig.add_subplot(gs1[1])
    for row in a.prop_active_m:
        ax.plot(a.bins, row, color=scolors[1], alpha=0.5)
    ax.scatter(data['m']['bins'], data['m']['props'], color='k')
    ax.set_xlabel('Age')
    ax.set_ylabel('Share')
    ax.set_title('Proportion of males\nwho are sexually active')

    # Network degree
    a = sim.analyzers.networkdegree
    relationship_type = 'lifetime_partners'
    for ai, sex in enumerate(['f', 'm']):
        ax = fig.add_subplot(gs2[ai])
        counts = a.results[f'{relationship_type}_{sex}'].values
        bins = a.bins

        total = sum(counts)
        counts = counts / total
        counts[-2] = counts[-2:].sum()
        counts = counts[:-1]

        ax.bar(bins[:-1], counts, color=scolors[ai])
        ax.set_xlabel(f'Number of {relationship_type}')
        ax.set_title(f'Distribution of partners, {sex}')
        ax.set_ylim([0, 1])

        sex_counts = np.array(getattr(a, f'{relationship_type}_{sex}'))
        stats = f"Mean: {np.mean(sex_counts):.1f}\n"
        stats += f"Median: {np.median(sex_counts):.1f}\n"
        stats += f"Std: {np.std(sex_counts):.1f}\n"
        stats += f"%>20: {np.count_nonzero(sex_counts >= 20) / total * 100:.2f}\n"
        ax.text(10, 0.3, stats)

    # Plot age differences
    ax = fig.add_subplot(gs3[0])
    a = sim.analyzers.partner_age_diff
    ax.hist(list(a.age_diffs.values()), label=list(a.age_diffs.keys()), bins=30, edgecolor='black', alpha=0.7)
    ax.legend()
    ax.set_xlabel('Age Difference (years)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Age Differences Between Partners\n in {a.year} (Male Age - Female Age)')

    # Save
    fig.tight_layout()
    pl.savefig(f"figures/figS_network.png", dpi=100)

    print('Done.')
