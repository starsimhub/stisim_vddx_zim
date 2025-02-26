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
    clist = sc.vectocolor([.5, .8, 1])
    clist = [clist[0], clist[1], clist[2]]
    colors = sc.objdict(treat50=clist[0], treat80=clist[1], treat100=clist[2])

    # Plot 1: Health
    ax = axes[0]
    sns.boxplot(data=hdf, x="disease", y="infections", hue="scenario", palette=clist, ax=ax)
    ax.legend(frameon=False)
    ax.set_title('% reduction in infections, 2027-2040')
    ax.set_ylim(-20, 100)

    # Plot 2: Treatment
    ax = axes[1]
    sns.boxplot(data=tdf, x="treatment", y="overtreatments", hue="scenario", palette=clist, ax=ax)
    ax.set_title('% reduction in overtreatment, 2027-2040')
    ax.set_ylim(0, 100)
    ax.get_legend().set_visible(False)

    fig.tight_layout()
    pl.savefig(f"figures/fig5_impact.png", dpi=100)

    print('Done!')





