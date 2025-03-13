"""
Plot HIV calibration
"""

# Import packages
import sciris as sc
import pylab as pl
import pandas as pd
from utils import set_font


# %% Plotting functions
def plot_calibration(calib, start_year=2000, end_year=2025):

    set_font(size=30)
    fig, axes = pl.subplots(2, 2, figsize=(20, 12))
    axes = axes.ravel()
    pn = 0
    result_map = {'new_infections': 'Infections', 'new_deaths': 'Deaths', 'n_infected': 'PLHIV', 'prevalence': 'Prevalence'}  #, 'bv': 'Other'}

    # Incidence
    for resname, reslabel in result_map.items():
        ax = axes[pn]

        # Pull out model results and data
        res = calib.sim_results
        target_data = calib.data['hiv_'+resname]
        resname = 'hiv_'+resname

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

        ax.set_title(reslabel)
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)

        pn += 1

    fig.tight_layout()
    pl.savefig(f"figures/hiv_calib.png", dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    n_results = 100

    percentile_pairs = [[.025, .975], [.05, .95], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
    percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]

    make_stats = True
    do_plot = True

    if make_stats:
        cal = sc.loadobj('results/zim_hiv_calib.obj')
        df = cal.resdf
        df_stats = df.groupby(df.time).describe(percentiles=percentiles)
        sc.saveobj(f'results/zim_hiv_calib_stats.df', df_stats)

    if do_plot:
        from plot_sims import plot_hiv_sims
        df_stats = sc.loadobj('results/zim_hiv_calib_stats.df')
        plot_hiv_sims(df_stats, start_year=2000, end_year=2025, which='multi', percentile_pairs=percentile_pairs, title='hiv_calib')

    best_pars = cal.df.iloc[0].to_dict()

    # Get posteriors
    par_stats = cal.df.describe(percentiles=[0.05, 0.95])
    pars = [p for p in par_stats.columns if p not in ['index', 'mismatch']]
    for p in pars:
        print(f'{p}: {par_stats[p]["mean"]:.3f} ({par_stats[p]["5%"]:.3f}–{par_stats[p]["95%"]:.3f})')

    print('Done.')
