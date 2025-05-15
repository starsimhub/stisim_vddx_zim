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
import pandas as pd
import utils as ut
from model import make_sim, make_sim_pars


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [5000, 2][debug]  # How many trials to run for calibration
n_workers = [100, 1][debug]    # How many cores to use
storage = None
do_shrink = True  # Whether to shrink the calibration results
make_stats = True  # Whether to make stats
study_name = 'starsim_calibration'


def make_calibration(scenario, n_trials=None, n_workers=None, constrain=False):

    # Define the calibration parameters
    ckw = dict(suggest_type='suggest_float')
    calib_par_dict = dict(
        default=dict(
            ng_p_symp=dict(low=0.1, high=0.2, guess=0.15, **ckw),
            ct_p_symp=dict(low=0.2, high=0.3, guess=0.25, **ckw),
            tv_p_symp=dict(low=0.15, high=0.75, guess=0.45, **ckw),
            p_symp_care=dict(low=0.25, high=0.75, guess=0.5, **ckw),
            ng_beta_m2f=dict(low=0.02, high=0.25, guess=0.08, **ckw, log=True),
            ct_beta_m2f=dict(low=0.02, high=0.25, guess=0.06, **ckw, log=True),
            tv_beta_m2f=dict(low=0.02, high=0.25, guess=0.07, **ckw, log=True),
        ),
        treat50=dict(
            ng_p_symp=dict(low=0.15, high=0.25, guess=0.18, **ckw),
            ct_p_symp=dict(low=0.25, high=0.3, guess=0.28, **ckw),
            tv_p_symp=dict(low=0.5, high=0.75, guess=0.6, **ckw),
            p_symp_care=dict(low=0.45, high=0.75, guess=0.6, **ckw),
        ),
        treat80=dict(
            ng_p_symp=dict(low=0.11, high=0.19, guess=0.15, **ckw),
            ct_p_symp=dict(low=0.22, high=0.27, guess=0.25, **ckw),
            tv_p_symp=dict(low=0.3, high=0.6, guess=0.45, **ckw),
            p_symp_care=dict(low=0.35, high=0.65, guess=0.5, **ckw),
        ),
        treat100=dict(
            ng_p_symp=dict(low=0.09, high=0.15, guess=0.12, **ckw),
            ct_p_symp=dict(low=0.2, high=0.25, guess=0.22, **ckw),
            tv_p_symp=dict(low=0.15, high=0.5, guess=0.3, **ckw),
            p_symp_care=dict(low=0.25, high=0.55, guess=0.4, **ckw),
        ),
    )
    calib_pars = calib_par_dict['default']
    if constrain: calib_pars['p_symp_care'] = calib_par_dict[scenario]['p_symp_care']

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
        ng_n_infected=0,
        ct_n_infected=0,
        tv_n_infected=0,
        ng_new_infections=0,
        ct_new_infections=0,
        tv_new_infections=0,
        ng_prevalence=2,
        ct_prevalence_f_25_30=2,
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

    return sim, calib


def run_calibration(scenario, calib, n_trials=None, do_save=False, constrain=False):

    # Run the calibration
    printstr = f'Running calibration for {scenario}, {n_trials} trials'
    if constrain:
        printstr += ' with constraints'.upper()
    sc.heading(printstr)
    calib.calibrate(load=True)
    if do_save: sc.saveobj(f'results/zim_sti_calib_{scenario}.obj', calib)
    print(f'Best pars are {calib.best_pars}')

    return calib


if __name__ == '__main__':

    constrain = True  # Whether to constrain the p_symp_care parameter
    load_partial = False

    # Loop over scenarios and run calibrations for each
    for scenario in ['treat80']: #ut.scenarios:

        sc.heading(f'Running calibration: {scenario}')
        sim, calib = make_calibration(scenario, n_trials=n_trials, n_workers=n_workers, constrain=constrain)

        if load_partial:
            # Load a partially-run calibration study
            import optuna as op
            print(calib.run_args.study_name)
            study = op.load_study(storage=calib.run_args.storage, study_name=calib.run_args.study_name)
            output = study.optimize(calib.run_trial, n_trials=1)
            calib.best_pars = sc.objdict(study.best_params)
            calib.parse_study(study)
            print('Best pars:', calib.best_pars)

            # Tidy up
            calib.calibrated = True
            if not calib.run_args.keep_db:
                calib.remove_db()

        else:
            calib = run_calibration(scenario, calib, n_trials=n_trials, n_workers=n_workers, constrain=constrain)


        print(f'... finished calibration: {scenario}')
        print(f'Best pars are {calib.best_pars}')
        resfolder = 'results/'  #if not constrain else 'results/constrained'  # NB constrained not in repo

        # Save the results
        print('Shrinking and saving...')
        if do_shrink:
            sc.saveobj(f'{resfolder}/zim_sti_calib_{scenario}_BIG.obj', calib)
            calib = calib.shrink(n_results=int(n_trials//20))  # Save 5% best results
            sc.saveobj(f'{resfolder}/zim_sti_calib_{scenario}.obj', calib)
        else:
            sc.saveobj(f'{resfolder}/zim_sti_calib_{scenario}.obj', calib)
        # Save the parameter dataframe
        sc.saveobj(f'{resfolder}/zim_sti_pars_{scenario}.df', calib.df)

        if make_stats:
            print('Making stats...')
            from utils import percentiles
            df = calib.resdf
            df_stats = df.groupby(df.time).describe(percentiles=percentiles)
            sc.saveobj(f'{resfolder}/zim_sti_calib_stats_{scenario}.df', df_stats)
            par_stats = calib.df.describe(percentiles=[0.05, 0.95])
            sc.saveobj(f'{resfolder}/zim_sti_par_stats_{scenario}.df', par_stats)

    print('Done!')


