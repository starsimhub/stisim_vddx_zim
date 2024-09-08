"""
Plotting utils
"""
import sciris as sc
import starsim as ss
import numpy as np


def set_font(size=None, font='Libertinus Sans'):
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


unneeded_results = [
    'pregnancy', 'deaths', 'structuredsexual', 'maternalnet', 'new_deaths', 'cum_deaths',
    'fsw_testing', 'other_testing', 'low_cd4_testing', 'art', 'vmmc', 'hivdx'
]


def process_results(sim):
    for key in unneeded_results:
        if key in sim.results.keys(): del sim.results[key]

    df = sim.export_df()
    df['year'] = np.floor(np.round(df.index, 1)).astype(int)
    df['seed'] = sim.pars.rand_seed
    df['scenario'] = sim.scenario

    ymean = df.groupby(by=['scenario', 'seed', 'year'])[mean_results].mean()
    ysums = df.groupby(by=['scenario', 'seed', 'year'])[sum_results].sum()

    return df
