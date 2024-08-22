# %% Imports and settings
import starsim as ss
import sciris as sc
import pandas as pd
from model import make_sim


def get_results_to_drop():
    """
    Define metrics to save
    """
    to_drop = [
        'pregnancy.pregnancies', 'pregnancy.births', 'pregnancy.cbr',
        'deaths.new', 'deaths.cumulative', 'deaths.cmr',
        'structuredsexual.share_active', 'structuredsexual.partners_f_mean', 'structuredsexual.partners_m_mean',
        'n_alive', 'new_deaths', 'cum_deaths'
    ]

    return to_drop


def run_sims(start=1990, seed=None, n_runs=None):
    sims = [make_sim(start=start, seed=seed + i, verbose=0.01) for i in range(n_runs)]
    sims = ss.parallel(sims).sims
    return sims


def process_output(sim):
    df_res = sim.export_df()
    output = sc.dataframe.from_dict(df_res)
    return output


def post_process_sims(sims, percentiles, to_drop, save_all=False):
    dfs = []
    for s, sim in enumerate(sims):
        output = process_output(sim)
        output['seed'] = sim.pars.rand_seed
        dfs += [output]
    bigdf = pd.concat(dfs)
    if save_all:
        sc.saveobj(f'results/multi_res.df', bigdf)  # NB this is a big file (~7MB), don't commit to repo!!!
    else:
        alldf = bigdf.drop(columns=to_drop)
        sc.saveobj(f'results/multi_res.df', alldf)
        df_stats = alldf.groupby(['yearvec']).describe(percentiles=percentiles)
        sc.saveobj(f'results/multi_res_stats.df', df_stats)


if __name__ == '__main__':

    # SETTINGS
    debug = False
    n_runs = [100, 2][debug]  # Number of runs when using multisim
    seed = 1
    sims = run_sims(start=1980, seed=seed, n_runs=n_runs)

    percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
    percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]
    to_drop = get_results_to_drop()
    post_process_sims(sims, percentiles, to_drop, save_all=False)
