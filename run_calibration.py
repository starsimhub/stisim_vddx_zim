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
n_workers = [40, 1][debug]  # How many cores to use
storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations


def run_calibration():

    sc.heading(f'Running calibration with ')

    # Define the calibration parameters
    calib_pars = dict(
        diseases=dict(
            ng=dict(
                beta_m2f=[0.08, 0.06, 0.18],
            ),
            ct=dict(
                beta_m2f=[0.04, 0.02, 0.1],
            ),
            tv=dict(
                beta_m2f=[0.02, 0.005, 0.05],
            ),
        ),
    )

    # Make the sim
    sim = make_sim()

    data = pd.read_csv('data/zimbabwe_calib.csv')

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        sim=sim,
        data=data,
        total_trials=4, n_workers=2, die=True
    )

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
