"""
Plot calibration
"""

# Import packages
import sciris as sc
import pylab as pl
import pandas as pd
from utils import set_font


# %% Plotting functions
def plot_calibration(calib, start_year=2000, end_year=2025):

    set_font(size=30)
    fig, axes = pl.subplots(3, 3, figsize=(20, 12))
    axes = axes.ravel()
    disease_map = {'ng': 'Gonorrhea', 'ct': 'Chlamydia', 'tv': 'Trich'}  #, 'bv': 'Other'}
    pn = 0

    # Incidence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]

        # Pull out model results and data
        res = calib.sim_results
        target_data = calib.data[dname+'.new_infections']
        resname = dname+'.new_infections'

        year = []
        values = []
        seed = []
        for run_num, run in enumerate(res):
            year += list(run['time'])
            values += list(run[resname])
            seed += [run_num]*len(run['time'])

        modeldf = pd.DataFrame({'time': year, resname: values, 'seed': seed})
        model_stats = modeldf.groupby(['time']).describe()
        model_plot = model_stats.iloc[(model_stats.index >= start_year) & (model_stats.index <= end_year)]

        x = model_plot.index
        medians = model_plot[(resname, '50%')]
        quartiles1 = model_plot[(resname, '25%')]
        quartiles3 = model_plot[(resname, '75%')]

        ax.plot(x, medians)
        ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)

        # Plot data
        data_plot = target_data.iloc[(target_data.index >= start_year) & (target_data.index <= end_year)]
        xdata = data_plot.index
        ydata = data_plot.values
        ax.scatter(xdata, ydata, color='k', marker='s', label='Data')

        ax.set_title(dlabel+' incidence')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)

        pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/sti_calib.png", dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    calib = sc.loadobj('results/calib.obj')
    plot_calibration(calib)

    print('Done.')
