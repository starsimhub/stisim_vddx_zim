"""
Run calibration for the STI model
"""

# Additions to handle numpy multithreading
import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# %% Imports and settings
import sciris as sc
import stisim as sti
import starsim as ss
import pandas as pd
from model import make_sim, make_scenpars


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [3000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None
do_shrink = True  # Whether to shrink the calibration results
make_stats = True  # Whether to make stats


def build_sim(sim, calib_pars):

   # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        v = pars['value']
        if 'beta' in k :
            sim.diseases[k[:2]].pars[k[3:]] = v
        elif 'p_symp' in k and k != 'p_symp_care':
            sim.diseases[k[:2]].pars[k[3:]][0] = v
        elif 'p_symp_care' in k:
            for dis in ['ng', 'ct', 'tv']:
                sim.diseases[dis].pars[k][0] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def run_calibration(scenario, n_trials=None, n_workers=None):

    # Define the calibration parameters
    calib_pars = dict(
        ng_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        ct_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        tv_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        ng_p_symp=dict(low=0.1, high=0.2, guess=0.15),
        ct_p_symp=dict(low=0.2, high=0.3, guess=0.25),
        tv_p_symp=dict(low=0.15, high=0.75, guess=0.45),
        p_symp_care=dict(low=0.25, high=0.75, guess=0.5, step=0.01),
    )

    # Extra results to save
    sres = sc.autolist()
    for dis in ['ng', 'ct', 'tv']:
        for res in ['prevalence', 'new_infections', 'n_infected']:
            for sk in ['', '_f', '_m']:
                sres += dis+'_'+res+sk

    # Make the sim
    scenpars = make_scenpars(scenario)
    sim = make_sim(scenario=scenario, **scenpars, start=1990, stop=2040, n_agents=10e3, verbose=-1, seed=1)
    data = pd.read_csv('data/zimbabwe_sti_data.csv')

    weights = dict(
        ng_new_infections=2,
        ct_new_infections=2,
        tv_new_infections=0.5,
    )

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        extra_results=sres,
        build_fn=build_sim,
        weights=weights,
        sim=sim,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=False, storage=storage, save_results=True,
    )

    calib.calibrate(load=True)
    sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return sim, calib


if __name__ == '__main__':

    # Settings
    sc.heading('Running STI calibration')
    sc.tic()
    scenarios = ['treat50', 'treat80', 'treat100']
    if len(scenarios) == 3:
        run_all = True

    # Run the calibration
    for scenario in scenarios:
        sim, calib = run_calibration(scenario, n_trials=n_trials, n_workers=n_workers)
        print(f'Best pars are {calib.best_pars}')

        # Save the results
        print('Shrinking and saving...')
        if do_shrink:
            calib = calib.shrink(n_results=500)
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
        else:
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)

        if make_stats:
            print('Making stats...')
            from utils import percentiles
            df = calib.resdf
            df_stats = df.groupby(df.time).describe(percentiles=percentiles)
            sc.saveobj(f'results/zim_sti_calib_stats_{scenario}.df', df_stats)
            par_stats = calib.df.describe(percentiles=[0.05, 0.95])
            sc.saveobj(f'results/zim_sti_par_stats_{scenario}.df', par_stats)

    # Make big dataframe
    if run_all:
        scenlabels = {'treat50': 'Treat-half', 'treat80':'Treat-most', 'treat100':'Treat-all'}
        dfs = sc.autolist()
        cs_dfs = sc.autolist()  # care seeking for VDS - not by disease
        results = sc.objdict()
        for scenario in scenarios:
            calib = sc.loadobj(f'results/zim_sti_calib_{scenario}.obj')
            df = calib.df[:500]
            df['scenario'] = scenlabels[scenario]
            df['ng_p_treat'] = df['ng_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100
            df['ct_p_treat'] = df['ct_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100
            df['tv_p_treat'] = df['tv_p_symp']*df['p_symp_care']*int(scenario.strip('treat'))/100

            df = df.loc[:, df.columns != 'p_symp_care']
            cs_df = calib.df[:500].loc[:, calib.df.columns.isin(['index', 'p_symp_care'])]
            cs_df['scenario'] = scenario
            dfs += df
            cs_dfs += cs_df

            # Save results
            results[scenario] = calib.sim_results[:50]

        df = pd.concat(dfs)
        cs_df = pd.concat(cs_dfs)

        # Melt dataframe to long form
        vars = ['ng_p_symp', 'ct_p_symp', 'tv_p_symp', 'ng_p_treat', 'ct_p_treat', 'tv_p_treat']
        dfm = df.melt(id_vars=['index', 'mismatch', 'scenario'], value_vars=vars, var_name='variable', value_name='value')

        # Add column for disease
        dfm['disease'] = dfm['variable'].apply(lambda x: x.split('_')[0].upper())
        dfm['par'] = dfm['variable'].apply(lambda x: x[3:])

        # Save results for figure 3
        sc.saveobj('results/zim_sti_calib_df.obj', dfm)  # TOO BIG TO ADD TO REPO
        sc.saveobj('results/zim_sti_care_seeking.obj', cs_df)
        sc.saveobj('results/zim_sti_calib_res.obj', results)

    print('Done!')


