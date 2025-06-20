"""
Plot degree of overtreatment by scenario
"""

# %% Imports and settings
import pandas as pd
import sciris as sc
import matplotlib.pyplot as pl
import utils as ut
import numpy as np


if __name__ == '__main__':

    do_process = True
    start_year = 2000
    end_year = 2025

    if do_process:

        # Load scenarios
        results = sc.objdict()
        scenarios = ut.scenarios
        for scenario in scenarios:
            df = sc.loadobj(f'results/{scenario}_sim.df')
            dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
            dfplot = dfplot.set_index('timevec')
            x = dfplot.index
            Y = [
                dfplot['syndromic_vds_new_tx0_f'],
                dfplot['syndromic_vds_new_tx1_f'],
                dfplot['syndromic_vds_new_tx2_f'],
                dfplot['syndromic_vds_new_tx3_f'],
            ]
            results[scenario] = sc.objdict(
                x=x,
                y=Y,
                labels=["0", "1", "2", "3"],
                colors=sc.vectocolor(4, reverse=True),
            )
        # Save results
        sc.saveobj('results/tx.obj', results)

    else:
        # Load results
        results = sc.loadobj('results/tx.obj')

    # Make plot
    ut.set_font(size=25)
    legend_font = 20
    n_plots = len(results)
    fig, axes = pl.subplots(1, 2, figsize=(20, 5))

    # Coinfection plot
    ax = axes[0]
    df = sc.loadobj(f'results/treat100_sim.df')
    dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
    dfplot = dfplot.set_index('timevec')
    x = dfplot.index
    Y = [
        dfplot['syndromic_vds_new_sti1_f'],
        dfplot['syndromic_vds_new_sti2_f'],
        dfplot['syndromic_vds_new_sti3_f'],
        dfplot['syndromic_vds_new_sti4_f'],
    ]
    labels = ["1 infection", "2 infections", "3 infections", "4 infections"]
    ax.stackplot(x, *Y, baseline='zero', labels=labels, colors=sc.vectocolor(4, reverse=True, cmap='plasma'))
    ax.set_title('Coinfection among VDS care seekers')
    ax.legend(frameon=False, prop={'size': legend_font}, loc='upper left', ncol=2)
    ax.set_ylim([0, 320_000])
    sc.SIticks(ax=ax)

    # Treatment plot
    ax = axes[1]
    colors = sc.vectocolor(4, reverse=True)
    ally = np.array([[np.sum(res.y[i])/np.sum(res.y) for i in range(4)] for res in results.values()])
    x = [0, 1, 2, 3, 4]  # x positions for the bars

    bottom = np.zeros(len(x))
    for tn in range(4):
        ax.bar(x, np.append(ally[:, tn], 0), bottom=bottom, color=colors[tn], label=tn)
        bottom += np.append(ally[:, tn], 0)
    ax.set_title('Proportion of care-seeking women\nreceiving multiple treatments, 2005-2025')
    ax.set_xticks([0,1,2,3,4], ['Treat-few', 'Treat-half', 'Treat-most', 'Treat-all', ''])
    ax.legend(frameon=False, prop={'size': legend_font}, loc='lower right', title='# treatments\nprescribed', title_fontsize=legend_font)

    fig.tight_layout()
    fig_name = 'figures/tx.png'
    sc.savefig(fig_name, dpi=100)

    print('Done!')





