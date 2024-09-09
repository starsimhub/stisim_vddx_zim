# %% Imports and settings
import starsim as ss
import sciris as sc
import pandas as pd
from model import make_sim
import stisim as sti
from utils import unneeded_results


def run_sims(start=1990, seed=None, n_runs=None):
    sims = [make_sim(start=start, seed=seed + i, verbose=0.01, end=2040) for i in range(n_runs)]
    sims = ss.parallel(sims).sims

    print("Processing... ")
    percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
    percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        sdf['seed'] = sim.pars.rand_seed
        sdf['scenario'] = sim.scenario
        dfs += [sdf]
    df = pd.concat(dfs)
    df_stats = df.describe(percentiles=percentiles)

    return sims, df, df_stats


if __name__ == '__main__':

    # SETTINGS
    debug = False
    n_runs = [100, 2][debug]  # Number of runs when using multisim
    seed = 1
    sims, df, df_stats = run_sims(start=1980, seed=seed, n_runs=n_runs)
    sc.saveobj(f'results/multi_res_stats.df', df_stats)

