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
import starsim as ss
import stisim as sti
import pandas as pd
import utils as ut
from model import make_sim, make_sim_pars


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [1000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None
do_shrink = True  # Whether to shrink the calibration results
make_stats = True  # Whether to make stats


def run_calibration(scenario, n_trials=None, n_workers=None, do_save=False):

    # Define the calibration parameters
    calib_par_dict = dict(
        treat50=dict(
            ng_p_symp=dict(low=0.15, high=0.2, guess=0.18),
            ct_p_symp=dict(low=0.25, high=0.3, guess=0.28),
            tv_p_symp=dict(low=0.5, high=0.75, guess=0.6),
            p_symp_care=dict(low=0.35, high=0.75, guess=0.55),
        ),
        treat80=dict(
            ng_p_symp=dict(low=0.12, high=0.8, guess=0.15),
            ct_p_symp=dict(low=0.22, high=0.27, guess=0.25),
            tv_p_symp=dict(low=0.3, high=0.6, guess=0.45),
            p_symp_care=dict(low=0.25, high=0.65, guess=0.45),
        ),
        treat100=dict(
            ng_p_symp=dict(low=0.1, high=0.5, guess=0.12),
            ct_p_symp=dict(low=0.2, high=0.25, guess=0.22),
            tv_p_symp=dict(low=0.15, high=0.5, guess=0.3),
            p_symp_care=dict(low=0.15, high=0.55, guess=0.35),
        ),
    )
    calib_pars = calib_par_dict[scenario]
    beta_pars = dict(
        ng_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        ct_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        tv_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
        # ng_dur=dict(low=6, high=10, guess=8, step=0.5),
        # ct_dur=dict(low=13, high=21, guess=15, step=0.5),
    )
    calib_pars = sc.mergedicts(calib_pars, beta_pars)
    # calib_pars = dict(
    #     # ng_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
    #     # ct_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
    #     # tv_beta_m2f=dict(low=0.02, high=0.2, guess=0.05),
    #     # ng_rel_beta_m2f=dict(low=1.5, high=3, guess=2),
    #     # ct_rel_beta_m2f=dict(low=1.5, high=3, guess=2),
    #     # tv_rel_beta_m2f=dict(low=1.5, high=3, guess=2),
    #     ng_p_symp=dict(low=0.1, high=0.2, guess=0.15),
    #     ct_p_symp=dict(low=0.2, high=0.3, guess=0.25),
    #     tv_p_symp=dict(low=0.15, high=0.75, guess=0.45),
    #     ng_dur=dict(low=6, high=10, guess=8, step=0.5),
    #     ct_dur=dict(low=13, high=21, guess=15, step=0.5),
    #     p_symp_care=dict(low=0.25, high=0.75, guess=0.5),
    # )

    # Extra results to save
    sres = sc.autolist()
    for dis in ['ng', 'ct', 'tv']:
        for res in ['prevalence', 'new_infections', 'n_infected']:
            for sk in ['', '_f', '_m']:
                sres += dis+'_'+res+sk

    # Make the sim
    sim = make_sim(scenario=scenario, use_calib=False, start=1990, stop=2040, verbose=-1, seed=1)
    data = pd.read_csv('data/zimbabwe_sti_data.csv')

    weights = dict(
        # ng_n_infected=0,
        # ct_n_infected=0,
        # tv_n_infected=0,
        # ng_new_infections=0,
        # ct_new_infections=0,
        # tv_new_infections=0,
        ng_prevalence=2,
        ct_prevalence=2,
        tv_prevalence=1,
    )

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        extra_results=sres,
        build_fn=make_sim_pars,
        weights=weights,
        sim=sim,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=False, storage=storage, save_results=True,
    )

    calib.calibrate(load=True)
    if do_save: sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return sim, calib


if __name__ == '__main__':

    # Settings
    sc.heading('Running STI calibration')

    # Run the calibration
    for scenario in ut.scenarios[1:]:
        sim, calib = run_calibration(scenario, n_trials=n_trials, n_workers=n_workers)
        print(f'Best pars are {calib.best_pars}')

        # Save the results
        print('Shrinking and saving...')
        if do_shrink:
            calib = calib.shrink(n_results=int(n_trials//20))  # Save 5% best results
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
        else:
            sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
        # Save the parameter dataframe
        sc.saveobj(f'results/zim_sti_pars_{scenario}.df', calib.df)

        if make_stats:
            print('Making stats...')
            from utils import percentiles
            df = calib.resdf
            df_stats = df.groupby(df.time).describe(percentiles=percentiles)
            sc.saveobj(f'results/zim_sti_calib_stats_{scenario}.df', df_stats)
            par_stats = calib.df.describe(percentiles=[0.05, 0.95])
            sc.saveobj(f'results/zim_sti_par_stats_{scenario}.df', par_stats)

    print('Done!')


