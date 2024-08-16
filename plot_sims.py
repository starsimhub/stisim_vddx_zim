# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font

location = 'zimbabwe'


def plot_sims(df, start_year=2000, end_year=2025, percentile_pairs=[[.1, .99]], title='sim_plots'):
    """ Create quantile plots """
    set_font(size=20)
    fig, axes = pl.subplots(2, 5, figsize=(18, 9))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
    dfplot['year'] = np.floor(np.round(dfplot.index, 1)).astype(int)

    # Population size
    pn = 0
    ax = axes[pn]
    pn += 1

    # Adult syphilis prevalence
    ax = axes[pn]
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + ".png", dpi=100)

    return fig


if __name__ == '__main__':

    sims = plot_sims()

