"""
Test calibration
"""

#%% Imports and settings
import starsim as ss
import sciris as sc
import stisim as sti
import numpy as np
import pylab as pl
import pandas as pd
from model import make_sim

do_plot = 1
do_save = 0
n_agents = 2e3


#%% Define the tests


def test_calibration():

    sc.heading('Testing calibration')

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

    sim, calib = test_calibration()
    sc.saveobj('results/calib.obj', calib)


    sc.toc(T)
    print('Done.')
