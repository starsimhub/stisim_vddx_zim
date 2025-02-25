# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti

# From this repo
from model import make_sim, make_scenpars


def load_calib_pars(scenario=None, calib=None, i=None):
    scenpars = make_scenpars(scenario)
    raw_calib_pars = calib.df.iloc[i].to_dict()

    # Overwrite
    diseases = ['ng', 'ct', 'tv']
    for disease in diseases:
        scenpars['p_symp'][disease] = raw_calib_pars[f'{disease}_p_symp']
        scenpars['p_symp_care'][disease] = raw_calib_pars['p_symp_care']
        scenpars['stipars'][disease]['beta_m2f'] = raw_calib_pars[f'{disease}_beta_m2f']
    return


def run_syndromic_scens(scenarios, stop=2040, parallel=True):
    """
    Run analyses
    """
    sims = sc.autolist()
    for scenario in scenarios:
        for pocstr in ['', 'poc']:
            scenname = scenario + pocstr
            calib = sc.loadobj(f'results/zim_sti_calib_{scenname}.obj')
            for i in range(n_scen_runs):
                print(f"Making sim: {scenname=}, param set {i+1}/{n_scen_runs}")
                scenpars = load_calib_pars(scenario=scenario, calib=calib, i=i)
                sim = make_sim(**scenpars, scenario=scenname, verbose=-1, stop=stop)
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


def process_results(df):
    res = pd.DataFrame()
    for scen in ['treat100', 'treat80', 'treat50']:
        for seed in df.seed.unique():
            thisdf = df.loc[(df.seed == seed) & (df.scenario.str.contains(scen))]
            for dis in ['ng', 'ct', 'tv']:
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
                res[dis+'_infections'] = (soc - poc)/soc

    return res


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    n_scen_runs = [10, 1][debug]  # Number of parameter sets to run per scenario
    scenarios = ['treat50', 'treat80', 'treat100']  #, 'panel']

    # Run analyses
    sims, df = run_syndromic_scens(scenarios, parallel=True, stop=2040)
    sc.saveobj('results/synd_scens.obj', df)
    print('Done!')

    # Load
    # df = sc.loadobj('results/synd_scens.obj')






