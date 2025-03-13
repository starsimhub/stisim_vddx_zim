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
n_trials = [1000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None
do_shrink = True  # Whether to shrink the calibration results


def build_sim(sim, calib_pars):

   # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        v = pars['value']
        if 'beta' in k :
            sim.diseases[k[:2]].pars[k[3:]] = v
        elif  'dur' in k:
            # sim.diseases[k[:2]].pars[k[3:]][0][0] = ss.dur(v, 'month')
            sim.diseases[k[:2]].pars['dur_symp2clear'][0][0] = ss.dur(v, 'month')
            sim.diseases[k[:2]].pars['dur_asymp2clear'][0][0] = ss.dur(v, 'month')
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
        ng_beta_m2f=dict(low=0.05, high=0.5, guess=0.1),
        ct_beta_m2f=dict(low=0.05, high=0.5, guess=0.1),
        tv_beta_m2f=dict(low=0.05, high=0.5, guess=0.1),
        ng_p_symp=dict(low=0.1, high=0.2, guess=0.15),
        ct_p_symp=dict(low=0.2, high=0.3, guess=0.25),
        tv_p_symp=dict(low=0.15, high=0.75, guess=0.45),
        ng_dur=dict(low=6, high=10, guess=8, step=0.5),
        ct_dur=dict(low=13, high=18, guess=15, step=0.5),
        # ng_dur_asymp2clear=dict(low=6, high=12, guess=8, step=0.5),
        # ct_dur_asymp2clear=dict(low=13, high=20, guess=15, step=0.5),
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
    sim = make_sim(scenario=scenario, **scenpars, start=1990, stop=2040, n_agents=5e3, verbose=-1, seed=1)
    data = pd.read_csv('data/zimbabwe_sti_data2.csv')

    weights = {'tv_new_infections': 0.2}

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

    for scenario in ['treat80']:  #, 'treat100']:  # 'treat50'
        sim, calib = run_calibration(scenario, n_trials=n_trials, n_workers=n_workers)
        if do_shrink:
            cal = calib.shrink(n_results=500)
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', cal)
        else:
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
        print(f'Best pars are {calib.best_pars}')
    print('Done!')


