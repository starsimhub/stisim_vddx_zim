# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

# From this repo
from model import make_sim
from utils import unneeded_results


def run_analyses(scenarios, end=2040, bv_range=None, parallel=True):
    """
    Run analyses
    """
    sims = sc.autolist()
    for scen in scenarios:
        for bv_beta_m2f in bv_range:
            for i in range(n_scen_runs):
                print(f"Making sim: {scen=}, seed={seed + i}")
                sim = make_sim(seed=seed + i, bv_beta_m2f=bv_beta_m2f, scenario=scen, verbose=0.01, end=end)
                sim.label = scen + str(i)
                sims += sim

    if parallel:
        sims = ss.parallel(sims).sims
    else:
        for sim in sims:
            sim.run()

    print("Processing... ")
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        sdf['seed'] = sim.pars.rand_seed
        sdf['scenario'] = sim.scenario
        sdf['bv_beta'] = sim.bv_beta_m2f
        sdf['vd_prev'] = sim.results.total_symptomatic.symp_prev_f[-1]
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [1, 1][debug]  # Number of seeds per scenarios
    scenarios = ['soc', 'panel']  #, 'panel']

    # Run analyses
    sims, df = run_analyses(scenarios, bv_range=[0.12, 0.15, 0.2], parallel=True, end=2040)
    sc.saveobj('results/scens.obj', df)






