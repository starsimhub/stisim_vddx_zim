"""
Run syndromic scenarios

Run on VMs - should take approx 10 min
"""

# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss

# From this repo
from model import make_sim
import utils as ut

def run_syndromic_scens(scenarios, stop=2040, parallel=True):
    """
    Run analyses
    """
    sc.heading("Making sims... ")

    sims = sc.autolist()
    for scenario in scenarios:
        for pocstr in ['', 'poc']:
            scenname = scenario + pocstr
            poc = True if pocstr == 'poc' else False
            for i in range(n_scen_runs):
                print(f"Making sim: {scenname=}, param set {i+1}/{n_scen_runs}")
                sim = make_sim(use_calib=True, par_idx=i, scenario=scenname, poc=poc, verbose=1/120, stop=stop)
                sim.label = scenname + str(i)
                sim.parset = i
                sims += sim

    sc.heading(f"Running {n_scen_runs} sims... ")
    if parallel:
        sims = ss.parallel(sims).sims
    else:
        for sim in sims:
            sim.run()

    sc.heading(f"Processing sims... ")
    print("WARNING, this will take a while...")
    dfs = []
    for s, sim in enumerate(sims):
        sdf = sim.to_df(resample='year', use_years=True, sep='.')
        sdf['parset'] = sim.parset
        sdf['scenario'] = sim.scenario
        sdf['poc'] = 1 if 'poc' in sim.scenario else 0
        dfs += [sdf]
    df = pd.concat(dfs)
    print("Finished processing sims.")

    return sims, df


def process_results(df):

    sc.heading(f"Processing results... ")
    healthdfs = sc.autolist()
    treatdfs = sc.autolist()

    from utils import treatments, tx_labels
    from utils import txscenlabels as scen_labels

    for scen in ut.scenarios:
        for parset in df.parset.unique():
            thisdf = df.loc[(df.parset == parset) & (df.scenario.str.contains(scen))]

            for dis in ['ng', 'ct', 'tv']:
                hres = pd.DataFrame()
                hres['scenario'] = [scen_labels[scen]]
                hres['parset'] = [parset]
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.timevec > 2027)][dis+'.new_infections_f'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.timevec > 2027)][dis+'.new_infections_f'].sum()
                hres['disease'] = [dis.upper()]
                hres['infections'] = [(soc - poc)/soc*100]
                healthdfs += hres

            for tx in ['ng_tx', 'ct_tx', 'metronidazole']:
                tres = pd.DataFrame()
                tres['scenario'] = [scen_labels[scen]]
                tres['parset'] = [parset]
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.timevec > 2027)][tx+'.new_treated_unnecessary_f'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.timevec > 2027)][tx+'.new_treated_unnecessary_f'].sum()
                tres['treatment'] = tx_labels[tx]
                tres['overtreatments'] = [(soc - poc)/soc*100]
                treatdfs += tres

    healthdf = pd.concat(healthdfs)
    treatdf = pd.concat(treatdfs)

    # Overtreatment stats - proportion reduction
    results = [tx+'.new_treated_unnecessary_f' for tx in treatments]
    results += [tx+'.new_treated_f' for tx in treatments]
    results += ['parset', 'scenario', 'timevec', 'poc']

    odf = df.loc[:, df.columns.isin(['timevec']+results)]
    for tx in treatments:
        odf[tx+'.overtx'] = odf[tx+'.new_treated_unnecessary_f']/df[tx+'.new_treated_f']

    # Melt dataframe to long form
    dfm = odf.melt(id_vars=['timevec', 'scenario'], var_name='variable', value_name='value')
    dfm['treatment'] = dfm['variable'].apply(lambda x: x.split('.')[0])

    return healthdf, treatdf, dfm


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [50, 1][debug]  # Number of parameter sets to run per scenario
    to_run = [
        'run_syndromic_scens',
        'process_results',
    ]

    # Imports
    from utils import scenarios

    if 'run_syndromic_scens' in to_run:
        # Run analyses
        sims, df = run_syndromic_scens(scenarios, parallel=True, stop=2041)
        sc.saveobj('results/synd_scens.obj', df)  # Don't commit to repo

    if 'process_results' in to_run:
        # Load
        df = sc.loadobj('results/synd_scens.obj')
        healthdf, treatdf, overdf = process_results(df)
        sc.saveobj('results/synd_health.obj', healthdf)  # In repo
        sc.saveobj('results/synd_treat.obj', treatdf)  # In repo
        sc.saveobj(f'results/overtx.obj', overdf)  # In repo

    print('Done!')


