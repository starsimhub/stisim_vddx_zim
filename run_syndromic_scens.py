# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

# From this repo
from model import make_sim, make_scenpars


def run_syndromic_scens(scenarios, stop=2040, parallel=True):
    """
    Run analyses
    """
    sims = sc.autolist()
    for scen in scenarios:
        for i in range(n_scen_runs):
            print(f"Making sim: {scen=}, seed={seed + i}")
            for pocstr in ['', 'poc']:
                scenname = scen + pocstr
                scenpars = make_scenpars(scenname)
                sim = make_sim(seed=seed + i, **scenpars, scenario=scenname, verbose=0.01, stop=stop)
                sim.label = scenname + str(i)
                sims += sim

    if parallel:
        sims = ss.parallel(sims).sims
    else:
        for sim in sims:
            sim.run()

    print("Processing... ")
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sim.to_df(resample='year', use_years=True, sep='.')
        sdf['seed'] = sim.pars.rand_seed
        sdf['scenario'] = sim.scenario
        sdf['poc'] = 1 if 'poc' in sim.scenario else 0
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [10, 1][debug]  # Number of seeds per scenarios
    scenarios = ['treat100', 'treat80', 'treat50']  #, 'panel']

    # Run analyses
    sims, df = run_syndromic_scens(scenarios, parallel=True, stop=2040)
    sc.saveobj('results/synd_scens.obj', df)






