"""
Plot HIV calibration
"""

# Import packages
import sciris as sc
from plot_sims import plot_hiv_sims, plot_sti_sims
from utils import percentile_pairs


# %% Run as a script
if __name__ == '__main__':

    which = ['hiv', 'sti'][1]
    scenario = 'treat100'

    # Load files - should all be committed to the repository
    df_filename = f'results/zim_{which}_calib_stats' + (f'_{scenario}' if which == 'sti' else '') + '.df'
    par_filename = f'results/zim_{which}_par_stats' + (f'_{scenario}' if which == 'sti' else '') + '.df'
    df_stats = sc.loadobj(df_filename)
    par_stats = sc.loadobj(par_filename)

    # Plot settings
    plot_kwargs = dict(
        start_year=2000,
        end_year=2025,
        which='multi',
        fext=scenario,
        percentile_pairs=percentile_pairs,
        title=f'{which}_calib',
    )

    # Plot
    if which == 'hiv':
        plot_hiv_sims(df_stats, **plot_kwargs)
    elif which == 'sti':
        plot_sti_sims(df_stats, **plot_kwargs)

    # Print posteriors
    pars = [p for p in par_stats.columns if p not in ['index', 'mismatch']]
    for p in pars:
        print(f'{p}: {par_stats[p]["mean"]:.3f} ({par_stats[p]["5%"]:.3f}–{par_stats[p]["95%"]:.3f})')

    print('Done.')
