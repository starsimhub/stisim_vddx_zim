"""
Test calibration
"""

# Additions to handle numpy multithreading
import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

#%% Imports and settings
import sciris as sc
import stisim as sti
import pandas as pd
from model import make_sim


debug = False  # If True, this will do smaller runs that can be run locally for debugging
do_save = True

# Run settings for calibration (dependent on debug)
n_trials = [1000, 10][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]  # How many cores to use
storage = ["mysql://hpvsim_user@localhost/stisim_db", None][debug]  # Storage for calibrations


def run_calibration():

    sc.heading(f'Running calibration')

    # Define the calibration parameters
    calib_pars = dict(
        beta_m2f_ng = dict(low=0.02, high=0.08, guess=0.06, path=('diseases', 'ng', 'beta_m2f')),
        beta_m2f_ct = dict(low=0.02, high=0.08, guess=0.06, path=('diseases', 'ct', 'beta_m2f')),
        beta_m2f_tv = dict(low=0.04, high=0.12, guess=0.06, path=('diseases', 'tv', 'beta_m2f')),
    )

    # Make the sim
    sim = make_sim(verbose=-1)

    data = pd.read_csv('data/zimbabwe_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        calib_pars = calib_pars,
        sim = sim,
        data = data,
        total_trials = n_trials,
        n_workers = n_workers,
        die = True,
        debug = False,
    )

    # Perform the calibration
    sc.printcyan('\nPeforming calibration...')
    calib.calibrate(confirm_fit=False)

    return sim, calib


#%% Run as a script
if __name__ == '__main__':

    T = sc.tic()

    sim, calib = run_calibration()
    from utils import shrink_calib
    cal = shrink_calib(calib, n_results=100)
    sc.saveobj('results/calib.obj', cal)

    sc.toc(T)
    print('Done.')
