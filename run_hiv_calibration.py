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
n_trials = [1000, 2][debug]  # How many trials to run for calibration
n_workers = [40, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None


def build_sim(sim, calib_pars):

    hiv = sim.diseases.hiv
    nw = sim.networks.structuredsexual

    # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        v = pars['value']
        if 'hiv_' in k:  # HIV parameters
            k = k.replace('hiv_', '')  # Strip off indentifying part of parameter name
            hiv.pars[k] = v
        elif 'nw_' in k:  # Network parameters
            k = k.replace('nw_', '')  # As above
            if 'pair_form' in k:
                nw.pars[k].set(v)
            else:
                nw.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def run_calibration(n_trials=None, n_workers=None, do_save=True):

    # Define the calibration parameters
    calib_pars = dict(
        hiv_beta_m2f=dict(low=0.01, high=0.10, guess=0.05),
        nw_prop_f0 = dict(low=0.55, high=0.9, guess=0.85),
        nw_prop_m0 = dict(low=0.50, high=0.9, guess=0.81),
        nw_f1_conc = dict(low=0.01, high=0.2, guess=0.01),
        nw_m1_conc = dict(low=0.01, high=0.2, guess=0.01),
        nw_p_pair_form = dict(low=0.4, high=0.9, guess=0.5),
    )

    # Make the sim
    sim = make_sim(scenario='soc', start=1990, stop=2030, n_agents=5e3)
    data = pd.read_csv('data/zimbabwe_hiv_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        build_fn = build_sim,
        sim=sim,
        data=data,
        total_trials=n_trials, n_workers=n_workers, die=True, reseed=False, storage=storage,
    )

    calib.calibrate(load=True)
    sc.saveobj(f'results/zim_calib.obj', calib)
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
        calib = sc.loadobj('results/zim_calib.obj')
        df = calib.df
        sc.saveobj(f'results/zim_calib_df.obj', df)
        res = calib.sim_results
        sc.saveobj(f'results/zim_calib_res.obj', res)

