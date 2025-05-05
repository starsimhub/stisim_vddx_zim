"""
Run syndromic scenarios

Run on VMs - should take approx 10 min
"""
import numpy as np

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
    results = ['new_infections', 'new_infections_f', 'new_false_neg', 'n_infected', 'n_infected_f']

    for s, sim in enumerate(sims):
        print(f"Processing sim {s+1}/{len(sims)}")
        sdfs = sc.autolist()
        for res in results:
            for disease in ['ng', 'ct', 'tv']:
                colname = f'{disease}.{res}'
                thisdf = sim.results[disease][res].to_df(resample='year', use_years=True, col_names=colname)
                sdfs += thisdf
        sdf = pd.concat(sdfs, axis=1)
        # sdf = sim.to_df(resample='year', use_years=True, sep='.')
        sdf['parset'] = sim.parset
        sdf['scenario'] = sim.scenario
        sdf['poc'] = 1 if 'poc' in sim.scenario else 0
        dfs += [sdf]
    df = pd.concat(dfs)

    # Process durations
    max_dur_dict = {'ng': 18, 'ct': 24, 'tv': 12}
    for disease in ['ng', 'ct', 'tv']:
        dur_dfs = sc.autolist()
        for s, sim in enumerate(sims):
            di = (sim.people[disease].dur_inf.notnan & sim.people.female).uids
            dur_inf = sim.people[disease].dur_inf[di]
            dur_hist = np.histogram(dur_inf, bins=np.arange(max_dur_dict[disease] + 1), density=True)
            n = len(dur_hist[0])
            dd = dict(
                dur_inf=dur_hist[0],
                months=dur_hist[1][:-1],
                disease=[disease]*n,
                parset=[sim.parset]*n,
                scenario=[sim.scenario]*n,
            )
            dur_dfs += pd.DataFrame(dd)
        dur_df = pd.concat(dur_dfs)
        sc.saveobj(f'results/dur_df_{disease}.obj', dur_df)

    print("Finished processing sims.")

    return sims, df


def process_results(df):

    sc.heading(f"Processing results... ")
    healthdfs = sc.autolist()
    treatdfs = sc.autolist()

    from utils import treatments, tx_labels
    from utils import txscenlabels as scen_labels
    flow_results = ['new_infections', 'new_infections_f', 'new_false_neg']
    stock_results = ['n_infected', 'n_infected_f']

    for scen in ut.scenarios:
        for parset in df.parset.unique():
            thisdf = df.loc[(df.parset == parset) & (df.scenario.str.contains(scen))]

            for dis in ['ng', 'ct', 'tv']:
                hres = pd.DataFrame()
                hres['scenario'] = [scen_labels[scen]]
                hres['parset'] = [parset]
                for fr in flow_results:
                    soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.index > 2027)][dis+'.'+fr].sum()
                    poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.index > 2027)][dis+'.'+fr].sum()
                    hres[fr] = [(soc - poc)/soc*100]
                for sr in stock_results:
                    soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.index == 2040)][dis+'.'+sr].values[0]
                    poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.index == 2040)][dis+'.'+sr].values[0]
                    hres[sr] = [(soc - poc)/soc*100]
                hres['disease'] = [dis.upper()]
                healthdfs += hres

            for tx in ['ng_tx', 'ct_tx', 'metronidazole']:
                tres = pd.DataFrame()
                tres['scenario'] = [scen_labels[scen]]
                tres['parset'] = [parset]
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.index > 2027)][tx+'.new_treated_unnecessary_f'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.index > 2027)][tx+'.new_treated_unnecessary_f'].sum()
                tres['treatment'] = tx_labels[tx]
                tres['overtreatments'] = [(soc - poc)/soc*100]
                treatdfs += tres

    healthdf = pd.concat(healthdfs)
    treatdf = pd.concat(treatdfs)

    # Overtreatment stats - proportion reduction
    results = [tx+'.new_treated_unnecessary_f' for tx in treatments]
    results += [tx+'.new_treated_f' for tx in treatments]
    results += ['parset', 'scenario', 'timevec', 'poc']

    odf = df.loc[:, df.columns.isin(results)]
    for tx in treatments:
        odf[tx+'.overtx'] = odf[tx+'.new_treated_unnecessary_f']/df[tx+'.new_treated_f']

    # Melt dataframe to long form
    dfm = odf.melt(id_vars=['scenario'], var_name='variable', value_name='value')
    dfm['treatment'] = dfm['variable'].apply(lambda x: x.split('.')[0])

    return healthdf, treatdf, dfm


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [2, 1][debug]  # Number of parameter sets to run per scenario
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


