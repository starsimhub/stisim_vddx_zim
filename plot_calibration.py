"""
Plot calibration
"""

# Import packages
import sciris as sc
import pylab as pl
import pandas as pd
from utils import set_font


# %% Plotting functions
def plot_calibration(calib, start_year=2000, end_year=2025, res_to_plot=100):

    set_font(size=30)
    fig, axes = pl.subplots(3, 4, figsize=(20, 12))
    axes = axes.ravel()
    disease_map = {'ng': 'NG', 'ct': 'CT', 'tv': 'TV'}  #, 'bv': 'Other'}
    result_map = {
        'prevalence_f_15_25': 'Prevalence F 15-25',
        'prevalence': 'Prevalence',
        'new_infections': 'Infections',
        'n_infected': 'Burden',
    }
    pn = 0

    # Prevalence in women 15-25
    for dname, dlabel in disease_map.items():
        for resname, reslabel in result_map.items():
            ax = axes[pn]

            # Pull out model results and data
            res = calib.sim_results[:res_to_plot]
            resname = dname+'_'+resname

            year = []
            values = []
            seed = []
            for run_num, run in enumerate(res):
                year += list(run['time'])
                values += list(run[resname])
                seed += [run_num]*len(run['time'])

            modeldf = pd.DataFrame({'time': year, resname: values, 'seed': seed})
            model_stats = modeldf.groupby(['time']).describe(percentiles=[.1, .9])
            model_plot = model_stats.iloc[(model_stats.index >= start_year) & (model_stats.index <= end_year)]

            x = model_plot.index
            medians = model_plot[(resname, '50%')]
            quartiles1 = model_plot[(resname, '10%')]
            quartiles3 = model_plot[(resname, '90%')]

            ax.plot(x, medians)
            ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)

            # Plot data
            if 'prevalence_f_15_25' in resname:
                target_data = calib.data[resname]
                data_plot = target_data.iloc[(target_data.index >= start_year) & (target_data.index <= end_year)]
                xdata = data_plot.index
                ydata = data_plot.values
                ax.scatter(xdata, ydata, color='k', marker='s', label='Data')

            ax.set_title(dlabel+' '+reslabel)
            ax.set_ylim(bottom=0)
            sc.SIticks(ax=ax)

            pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/sti_calib.png", dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    scenario = 'treat80'
    calib = sc.loadobj(f'results/zim_sti_calib_{scenario}.obj')
    plot_calibration(calib)

    print('Done.')
