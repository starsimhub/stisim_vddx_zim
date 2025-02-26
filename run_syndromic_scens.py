# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss

# From this repo
from model import make_sim, load_calib_pars


def run_syndromic_scens(scenarios, stop=2040, parallel=True):
    """
    Run analyses
    """
    sims = sc.autolist()
    for scenario in scenarios:
        for pocstr in ['', 'poc']:
            scenname = scenario + pocstr
            calib = sc.loadobj(f'results/zim_sti_calib_{scenario}.obj')
            for i in range(n_scen_runs):
                print(f"Making sim: {scenname=}, param set {i+1}/{n_scen_runs}")
                scenpars = load_calib_pars(scenario=scenname, calib=calib, i=i)
                sim = make_sim(**scenpars, scenario=scenname, verbose=-1, stop=stop)
                sim.label = scenname + str(i)
                sim.parset = i
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
        sdf['parset'] = sim.parset
        sdf['scenario'] = sim.scenario
        sdf['poc'] = 1 if 'poc' in sim.scenario else 0
        dfs += [sdf]
    df = pd.concat(dfs)

    return sims, df


def process_results(df):
    healthdfs = sc.autolist()
    treatdfs = sc.autolist()

    tx_labels = {'ng_tx':'NG', 'ct_tx':'CT', 'metronidazole':'MTNZ'}
    scen_labels = {'treat50':'Poor', 'treat80':'Imperfect', 'treat100':'Perfect'}

    for scen in ['treat50', 'treat80', 'treat100']:
        for parset in df.parset.unique():
            thisdf = df.loc[(df.parset == parset) & (df.scenario.str.contains(scen))]

            for dis in ['ng', 'ct', 'tv']:
                hres = pd.DataFrame()
                hres['scenario'] = [scen_labels[scen]]
                hres['parset'] = [parset]
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
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

    return healthdf, treatdf


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [100, 1][debug]  # Number of parameter sets to run per scenario
    scenarios = ['treat50', 'treat80', 'treat100']  #, 'panel']
    to_run = [
        # 'run_syndromic_scens',
        'process_results',
        # 'plot_results',
    ]

    if 'run_syndromic_scens' in to_run:
        # Run analyses
        sims, df = run_syndromic_scens(scenarios, parallel=True, stop=2040)
        sc.saveobj('results/synd_scens.obj', df)

    if 'process_results' in to_run:
        # Load
        df = sc.loadobj('results/synd_scens.obj')
        healthdf, treatdf = process_results(df)
        sc.saveobj('results/synd_health.obj', healthdf)
        sc.saveobj('results/synd_treat.obj', treatdf)


    print('Done!')





