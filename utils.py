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


def count(arr): return np.count_nonzero(arr)


def shrink_calib(calib, n_results=100):
    cal = sc.objdict()
    plot_indices = calib.df.iloc[0:n_results, 0].values
    cal.sim_results = [calib.sim_results[i] for i in plot_indices]
    cal.data = calib.data
    cal.df = calib.df.iloc[0:n_results, ]
    return cal


def get_y(df, which, rname):
    if which == 'single': y = df[rname]
    elif which == 'multi': y = df[(rname, '50%')]
    return y

