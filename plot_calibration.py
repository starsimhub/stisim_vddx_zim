"""
Plot calibration
"""

# Import packages
import sciris as sc
import numpy as np
import pylab as pl
import pandas as pd
import seaborn as sns

from utils import set_font



def plot_sti_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='sti_plots'):
    """ Create quantile plots """
    set_font(size=30)
    fig, axes = pl.subplots(3, 4, figsize=(25, 12))
    axes = axes.ravel()

    # Incidence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resname = dname+'.new_infections'
        data = disease_data[dname]
        if data is not None:
            ax.scatter(data.year, data[resname], label='Data', color='k')
        resnames = {'Total': dname+'.new_infections', 'Symptomatic': dname+'.new_symptomatic', 'Care seekers': dname+'.new_care_seekers'}
        for rlabel, rname in resnames.items():
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())

        ax.set_title(dlabel+' incidence')
        if pn == 2: ax.legend(frameon=False, prop={'size': 20})
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    # Burden
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resname = dname+'.n_infected'
        data = disease_data[dname]
        if data is not None:
            ax.scatter(data.year, data[resname], label='Data', color='k')
        resnames = {'Total': dname+'.n_infected', 'Symptomatic': dname+'.n_symptomatic'}
        for rlabel, rname in resnames.items():
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' burden')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    # Prevalence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]
        resnames = {'Total': dname+'.adult_prevalence', 'Symptomatic': dname+'.symp_adult_prevalence'}
        # if dname == 'ng':
        #     data = disease_data[dname]
        #     ax.scatter(data.year, 100*data['ng.adult_prevalence'], label='Data', color='k')
        for rlabel, rname in resnames.items():
            x = dfplot.index
            y = get_y(dfplot, which, rname)
            line, = ax.plot(x, y*100, label=rlabel)
            if which == 'multi':
                for idx, percentile_pair in enumerate(percentile_pairs):
                    yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                    yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                    ax.fill_between(x, yl* 100, yu* 100, alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title(dlabel+' prevalence (%)')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax=ax)
        pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig





# %% Plotting functions
def plot_calibration(calib):

    set_font(size=30)
    fig, axes = pl.subplots(3, 4, figsize=(25, 12))
    axes = axes.ravel()
    disease_map = {'ng': 'Gonorrhea', 'ct': 'Chlamydia', 'tv': 'Trich'}  #, 'bv': 'Other'}
    pn = 0

    # Incidence
    for dname, dlabel in disease_map.items():
        ax = axes[pn]

        # Pull out model results and data
        res = calib.sim_results
        target_data = calib.target_data[dname+'.new_infections']
        resname = dname+'.new_infections'

        year = []
        values = []
        seed = []
        for run_num, run in enumerate(res):
            year += target_data.index.values.tolist()
            values += list(run[resname])
            seed += [run_num]*len(target_data.index.values)

        modeldf = pd.DataFrame({'year': year, 'values': values, 'seed': seed})
        model_stats = modeldf.groupby(['year']).describe()
        x = modeldf.year
        medians = model_stats[(resname, '50%')]
        quartiles1 = model_stats[('signal', '25%')]
        quartiles3 = model_stats[('signal', '75%')]

        ax.plot(x, medians)
        ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)

        # Plot data
        xdata = target_data.index
        ydata = target_data.values
        ax.scatter(xdata, ydata, color='k', marker='s', label='Data')

    fig.tight_layout()
    pl.savefig(f"figures/sti_calib.png", dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    calib = sc.loadobj('results/calib.obj')
    plot_calibration(calib)

    print('Done.')
