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
