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
from model import make_sim


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [5000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None


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


def run_calibration(n_trials=None, n_workers=None):

    # Define the calibration parameters
    calib_pars = dict(
        ng_beta_m2f=dict(low=0.055, high=0.065, guess=0.06),
        ct_beta_m2f=dict(low=0.055, high=0.065, guess=0.06),
        tv_beta_m2f=dict(low=0.09, high=0.11, guess=0.10),
    )

    # Make the sim
    sim = make_sim(scenario='soc', start=1990, stop=2030, n_agents=5e3)
    data = pd.read_csv('data/zimbabwe_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        build_fn=build_sim,
        sim=sim,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=True, storage=storage, save_results=True,
    )

    calib.calibrate(load=True)
    sc.saveobj(f'results/zim_sti_calib.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return sim, calib


if __name__ == '__main__':

    to_run = [
        'run_calib',
        # 'load_calib'
    ]

    if 'run_calib' in to_run:
        sim, calib = run_calibration(n_trials=n_trials, n_workers=n_workers)

    if 'load_calib' in to_run:
        calib = sc.loadobj('results/zim_sti_calib.obj')
        df = calib.df
        sc.saveobj(f'results/zim_sti_calib_df.obj', df)
        res = calib.sim_results
        sc.saveobj(f'results/zim_sti_calib_res.obj', res)

