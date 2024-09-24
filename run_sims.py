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
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        sdf['seed'] = sim.pars.rand_seed
        sdf['scenario'] = sim.scenario
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


def run_sens_sims(start=1980, seed=None, n_runs=None, bv_range=None, end=2025):

    sims = sc.autolist()
    for bv_beta_m2f in bv_range:
        for i in range(n_runs):
            print(f"Making sim: {bv_beta_m2f=}, seed={seed + i}")
            sim = make_sim(start=start, seed=seed + i, bv_beta_m2f=bv_beta_m2f, verbose=0.01, end=end)
            sim.label = f"BV beta {bv_beta_m2f}: {str(i)}"
            sims += sim

    sims = ss.parallel(sims).sims

    print("Processing... ")
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        sdf['seed'] = sim.pars.rand_seed
        sdf['bv_beta'] = sim.bv_beta_m2f
        sdf['vd_prev'] = sim.results.total_symptomatic.symp_prev_f[-1]
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


if __name__ == '__main__':

    # SETTINGS
    debug = False
    n_runs = [1, 1][debug]  # Number of runs when using multisim
    seed = 1
    sims, df = run_sens_sims(start=1980, seed=seed, n_runs=n_runs, bv_range=[0.12, 0.15, 0.2])

    output = ''
    for disease in ['ng', 'ct', 'tv']:
        output += disease.upper()+": \n"
        for sim in sims:
            output += "New sim \n"

            a = sim.results[disease].new_treated_success_symp[-1]/sim.results[disease].new_treated[-1]
            b = sim.results[disease].new_treated_success_asymp[-1]/sim.results[disease].new_treated[-1]
            c = sim.results[disease].new_treated_failure[-1]/sim.results[disease].new_treated[-1]
            d = sim.results[disease].new_treated_unnecessary[-1]/sim.results[disease].new_treated[-1]
            output += str(a) + '\n'
            output += str(b) + '\n'
            output += str(c) + '\n'
            output += str(d) + '\n'




    # percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
    # percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]
    # df_stats = df.groupby(df.index).describe(percentiles=percentiles)
    # sc.saveobj(f'results/multi_res_stats.df', df_stats)

