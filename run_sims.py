# %% Imports and settings
import starsim as ss
import sciris as sc
import pandas as pd
from model import make_sim


def get_results_to_save():
    """
    Define metrics to save
    """
    to_save = [
        'seed',
        'yearvec',
        'ng.n_infected',
        'ng.new_infections',
        'ct.n_infected',
        'ct.new_infections',
        'tv.n_infected',
        'tv.new_infections',
        'vd.n_infected',
        'vd.new_infections',
        'hiv.new_infections',
        'hiv.new_deaths',
        'hiv.n_infected',
        'hiv.n_diagnosed',
        'hiv.n_on_art',
        'hiv.prevalence',
        'n_alive',
       ]

    return to_save


def run_sims(seed=None, n_runs=None):
    sims = [make_sim(seed=seed + i, verbose=0.01) for i in range(n_runs)]
    sims = ss.parallel(sims).sims
    return sims


def process_output(sim):
    df_res = sim.export_df()
    output = sc.dataframe.from_dict(df_res)
    return output


def post_process_sims(sims, percentiles, to_save, save_all=False):
    dfs = []
    for s, sim in enumerate(sims):
        output = process_output(sim)
        output['seed'] = sim.pars.rand_seed
        dfs += [output]
    bigdf = pd.concat(dfs)
    if save_all:
        sc.saveobj(f'results/multi_res.df', bigdf)  # NB this is a big file (~7MB), don't commit to repo!!!
    else:
        to_drop = [c for c in bigdf.columns if c not in to_save]
        alldf = bigdf.drop(columns=to_drop)
        sc.saveobj(f'results/multi_res.df', alldf)
        df_stats = alldf.groupby(['yearvec']).describe(percentiles=percentiles)
        sc.saveobj(f'results/multi_res_stats.df', df_stats)


if __name__ == '__main__':

    # SETTINGS
    debug = False
    n_runs = [100, 2][debug]  # Number of runs when using multisim
    seed = 1
    sims = run_sims(seed=seed, n_runs=n_runs)

    percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
    percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]
    to_save = get_results_to_save()
    post_process_sims(sims, percentiles, to_save, save_all=False)
