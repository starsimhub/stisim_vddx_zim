"""
Plot impact of POC test
Need to run run_syndromic_scens.py first to generate the result dataframes
"""

# %% Imports and settings
import sciris as sc
import matplotlib.pyplot as pl
import seaborn as sns
from utils import set_font


if __name__ == '__main__':

    hdf = sc.loadobj('results/synd_health.obj')
    tdf = sc.loadobj('results/synd_treat.obj')
    set_font(size=20)
    fig, axes = pl.subplots(1, 2, figsize=(20, 8))
    axes = axes.ravel()

    # Plot 1: Health
    ax = axes[0]
    sns.boxplot(data=hdf, x="disease", y="infections", hue="scenario", palette='viridis', ax=ax)
    ax.set_title('% reduction in infections, 2027-2040')
    ax.set_ylim(-10, 100)

    # Plot 2: Treatment
    ax = axes[1]
    sns.boxplot(data=tdf, x="treatment", y="overtreatments", hue="scenario", palette='viridis', ax=ax)
    ax.set_title('% reduction in overtreatment, 2027-2040')
    ax.set_ylim(0, 100)

    fig.tight_layout()
    pl.savefig(f"figures/fig5_impact.png", dpi=100)

    print('Done!')





