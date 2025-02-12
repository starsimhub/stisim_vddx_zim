"""
Run calibration for the HIV model
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
import pandas as pd
from model import make_sim, make_scenpars


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [5000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None
scenario = 'treat80'


def build_sim(sim, calib_pars):

    ng = sim.diseases.ng
    ct = sim.diseases.ct
    tv = sim.diseases.tv

    # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        v = pars['value']
        if 'ng_' in k:
            k = k.replace('ng_', '')
            ng.pars[k] = v
        elif 'ct_' in k:
            k = k.replace('ct_', '')
            ct.pars[k] = v
        elif 'tv_' in k:
            k = k.replace('tv_', '')
            tv.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def run_calibration(scenario, n_trials=None, n_workers=None):

    # Define the calibration parameters
    calib_pars = dict(
        ng_beta_m2f=dict(low=0.05, high=0.3, guess=0.06),
        ct_beta_m2f=dict(low=0.02, high=0.1, guess=0.05),
        tv_beta_m2f=dict(low=0.08, high=0.3, guess=0.10),
        ng_eff_condom=dict(low=0.5, high=0.9, guess=0.8),
        ct_eff_condom=dict(low=0.5, high=0.9, guess=0.8),
        tv_eff_condom=dict(low=0.5, high=0.9, guess=0.8),
    )

    # Make the sim
    scenpars = make_scenpars(scenario)
    sim = make_sim(scenario=scenario, **scenpars, start=1990, stop=2040, n_agents=5e3, verbose=-1, seed=1)
    data = pd.read_csv('data/zimbabwe_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        build_fn=build_sim,
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

    to_run = [
        'run_calib',
        # 'load_calib'
    ]

    if 'run_calib' in to_run:
        sim, calib = run_calibration(scenario, n_trials=n_trials, n_workers=n_workers)

    if 'load_calib' in to_run:
        calib = sc.loadobj('results/zim_sti_calib.obj')
        df = calib.df
        sc.saveobj(f'results/zim_sti_calib_df.obj', df)
        res = calib.sim_results
        sc.saveobj(f'results/zim_sti_calib_res.obj', res)

