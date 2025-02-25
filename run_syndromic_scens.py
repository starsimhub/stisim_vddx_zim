# %% Imports and settings
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti
import matplotlib.pyplot as pl
import seaborn as sns

# From this repo
from model import make_sim, make_scenpars
from utils import set_font


def load_calib_pars(scenario=None, calib=None, i=None):
    scenpars = make_scenpars(scenario)
    raw_calib_pars = calib.df.iloc[i].to_dict()

    # Overwrite
    diseases = ['ng', 'ct', 'tv']
    for disease in diseases:
        scenpars['p_symp'][disease] = raw_calib_pars[f'{disease}_p_symp']
        scenpars['p_symp_care'][disease] = raw_calib_pars['p_symp_care']
        scenpars['stipars'][disease]['beta_m2f'] = raw_calib_pars[f'{disease}_beta_m2f']
    return scenpars


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

    for scen in ['treat50', 'treat80', 'treat100']:
        for parset in df.parset.unique():
            thisdf = df.loc[(df.parset == parset) & (df.scenario.str.contains(scen))]

            for dis in ['ng', 'ct', 'tv']:
                hres = pd.DataFrame()
                hres['scenario'] = [scen]
                hres['parset'] = [parset]
                soc = thisdf.loc[(thisdf.poc == 0) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
                poc = thisdf.loc[(thisdf.poc == 1) & (thisdf.timevec > 2027)][dis+'.new_infections'].sum()
                hres['disease'] = [dis.upper()]
                hres['infections'] = [(soc - poc)/soc*100]
                healthdfs += hres

            for tx in ['ng_tx', 'ct_tx', 'metronidazole']:
                tres = pd.DataFrame()
                tres['scenario'] = [scen]
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
        'run_syndromic_scens',
        # 'process_results',
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

    if 'plot_results' in to_run:
        hdf = sc.loadobj('results/synd_health.obj')
        tdf = sc.loadobj('results/synd_treat.obj')
        set_font(size=20)
        fig, axes = pl.subplots(1, 2, figsize=(20, 8))
        axes = axes.ravel()

        # Plot 1: Health
        ax = axes[0]
        sns.boxplot(data=hdf, x="disease", y="infections", hue="scenario", palette='viridis', ax=ax)
        ax.set_title('% reduction in infections, 2027-2040')
        ax.set_ylim(-10, 100)

        # Plot 2: Treatment
        ax = axes[1]
        sns.boxplot(data=tdf, x="treatment", y="overtreatments", hue="scenario", palette='viridis', ax=ax)
        ax.set_title('% reduction in overtreatment, 2027-2040')
        ax.set_ylim(0, 100)

        fig.tight_layout()
        pl.savefig(f"figures/fig5_impact.png", dpi=100)

    print('Done!')





