# %% Imports and settings
import pandas as pd
import numpy as np
import sciris as sc
import starsim as ss
import pylab as pl
import stisim as sti
import seaborn as sns

# From this repo
from model import make_sim
from utils import unneeded_results

def run_analyses(scenarios, end=2040, parallel=True):
    """
    Run analyses
    """
    sims = sc.autolist()
    for scen in scenarios:
        for i in range(n_scen_runs):
            print(f"Making sim: {scen=}, seed={seed + i}")
            sim = make_sim(seed=seed + i, scenario=scen, verbose=0.01, end=end)
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
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [1, 1][debug]  # Number of seeds per scenarios
    scenarios = ['soc', 'panel']

    if True:
        # Run analyses
        sims, df = run_analyses(scenarios, parallel=False, end=2010)
        sc.saveobj('results/scens.obj', df)

    df = sc.loadobj('results/scens.obj')
    fig, axes = pl.subplots(2, 2, figsize=(15, 10))
    axes = axes.ravel()
    intv_year = 2000

    for pn, dis in enumerate(['ng', 'ct', 'tv', 'bv']):
        ax = axes[pn]
        # sns.lineplot(df, x=df.index, y=dis+".new_infections", hue="scenario", ax=ax)
        sns.lineplot(df, x=df.index, y=f"{dis}_tx.new_treated", hue="scenario", ax=ax)
    pl.show()







