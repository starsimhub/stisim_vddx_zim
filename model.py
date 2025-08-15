"""
Create a model with HIV plus 4 co-circulating discharging STIs:
    - chlamydia, gonorrhea, trichomoniasis, and other (BV+)
Used for evaluation of etiological tests compared to syndromic management.
"""
import scipy.special

# %% Imports and settings
import starsim as ss
import stisim as sti

from analyzers import total_symptomatic as ts
from hiv_model import make_hiv, make_hiv_intvs
from interventions import make_testing
from plot_sims import *
import utils as ut


def make_stis():
    ng = sti.Gonorrhea(eff_condom=0.7)
    ct = sti.Chlamydia(eff_condom=0.8)
    tv = sti.Trichomoniasis(eff_condom=0.8)
    bv = sti.SimpleBV()
    return ng, ct, tv, bv


def make_sim_pars(sim, calib_pars):
    """
    Update the simulation parameters with the calibration parameters
    """
    def set_par(k, sim):
        if sim.initialized:
            return sim.diseases[k[:2]].pars
        else:
            idx = [d.name for d in sim.pars.diseases].index(k[:2])
            return sim.pars.diseases[idx].pars

    for k, pars in calib_pars.items():  # Loop over the calibration parameters
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        if isinstance(pars, dict):
            v = pars['value']
        elif sc.isnumber(pars):
            v = pars
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

        if 'beta_m2f' in k:
            set_par(k, sim)[k[3:]] = v
        elif 'dur' in k:
            set_par(k, sim)['dur_symp2clear'][0][0] = ss.dur(v, 'month')
            set_par(k, sim)['dur_asymp2clear'][0][0] = ss.dur(v, 'month')
        elif 'p_symp' in k and k != 'p_symp_care':
            set_par(k, sim)[k[3:]][0] = v
        elif 'p_symp_care' in k:
            for dis in ['ng', 'ct', 'tv']:
                if sim.initialized:
                    sim.diseases[dis].pars[k][0] = v
                else:
                    didx = [d.name for d in sim.pars.diseases].index(dis)
                    sim.pars.diseases[didx].pars[k][0] = v
        elif k in ['index', 'mismatch']:
            continue
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def make_sim(seed=1, n_agents=None, dt=1/12, start=1990, stop=2030, debug=False, verbose=1/12, add_stis=True,
             scenario='treat100', use_calib=False, calib_folder=None, par_idx=0, poc=False, analyzers=None,
             analyze_network=False):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1985: 8.691e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(10e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_data = pd.read_csv(f'data/asfr.csv')
    pregnancy = ss.Pregnancy(unit='month', fertility_rate=fertility_data)
    death_data = pd.read_csv(f'data/deaths.csv')
    death = ss.Deaths(death_rate=death_data, rate_units=1)

    ####################################################################################################################
    # People and networks
    ####################################################################################################################
    ppl = ss.People(n_agents, age_data=pd.read_csv(f'data/age_dist_{start}.csv', index_col='age')['value'])
    sexual = sti.StructuredSexual(
        prop_f0=0.79,
        prop_m0=0.83,
        f1_conc=0.16,
        m1_conc=0.11,
        p_pair_form=0.58,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet(unit='month')

    ####################################################################################################################
    # Diseases
    ####################################################################################################################
    hiv = make_hiv()
    diseases = [hiv]
    if add_stis:
        ng, ct, tv, bv = make_stis()
        stis = [ng, ct, tv, bv]
        diseases += stis  # Add the STIs to the list of diseases

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    intvs = make_hiv_intvs()
    if add_stis:
        intvs += make_testing(ng, ct, tv, bv, scenario=scenario, poc=poc, stop=stop)
        connectors = [sti.hiv_ng(hiv, ng), sti.hiv_ct(hiv, ct), sti.hiv_tv(hiv, tv)]
    else:
        connectors = []

    if analyzers is None and add_stis:
        analyzers = [sti.sw_stats(diseases=['ng', 'ct', 'tv', 'hiv']), ts()]

    # Add network analyzers
    if analyze_network:
        analyzers = sc.autolist(analyzers)
        analyzers += sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])
        analyzers += sti.RelationshipDurations()
        analyzers += sti.DebutAge()
        analyzers += sti.partner_age_diff()

    sim = ss.Sim(
        dt=dt,
        rand_seed=seed,
        total_pop=total_pop,
        start=start,
        stop=stop,
        people=ppl,
        diseases=diseases,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=intvs,
        analyzers=analyzers,
        connectors=connectors,
        verbose=verbose,
    )

    # Store scenario
    sim.scenario = scenario

    # If using calibration parameters, update the simulation
    if use_calib:
        calibname = scenario.strip('poc')
        if calib_folder is None:
            calib_folder = 'results'
        pars_df = sc.loadobj(f'{calib_folder}/zim_sti_pars_{calibname}.df')
        calib_pars = pars_df.iloc[par_idx].to_dict()
        sim.init()
        sim = make_sim_pars(sim, calib_pars)
        print(f'Using calibration parameters for scenario {scenario} and index {par_idx}')

    return sim


def run_msim(scenarios=None, use_calib=True, calib_folder=None, n_pars=1, seed=1, debug=False, do_save=True):

    # Mave individual sims
    sims = sc.autolist()

    if scenarios is None:
        scenarios = ut.scenarios

    for scenario in scenarios:
        for par_idx in range(n_pars):
            sim = make_sim(scenario=scenario, use_calib=use_calib, calib_folder=calib_folder, par_idx=par_idx, seed=seed, debug=debug, start=1990, stop=2026)
            if use_calib:
                print('Using calibration parameters:')
                print(f'ng_p_symp: {sim.diseases.ng.pars.p_symp}')
                print(f'p_symp_care: {sim.diseases.ct.pars.p_symp_care}')
            sim.par_idx = par_idx
            sims += sim

    sims = ss.parallel(sims).sims

    if do_save:
        dfs = sc.autolist()
        for sim in sims:
            scenario = sim.scenario
            par_idx = sim.par_idx
            df = sim.to_df(resample='year', use_years=True, sep='.')
            df['res_no'] = par_idx
            df['scenario'] = scenario
            dfs += df
        df = pd.concat(dfs)
        sc.saveobj(f'results/msim.df', df)

        # for sim in sims:
        #     scenario = sim.scenario
        #     df = sim.to_df(resample='year', use_years=True, sep='.')
        #     sc.saveobj(f'results/{scenario}_sim.df', df)

    return sims


def save_stats(sims, resfolder='results', scenario='treat80'):

    # Epi stats: save for all runs
    dfs = sc.autolist()
    for sim in sims:

        scenario = sim.scenario
        par_idx = sim.par_idx
        df = sim.to_df(resample='year', use_years=True, sep='.')

        # Save age/sex epi results
        age_bins = sim.diseases.ng.age_bins
        sex_labels = {'f': 'Female', 'm': 'Male'}
        for disease in ['ng', 'ct', 'tv', 'hiv']:
            for sex in ['f', 'm']:
                dd = dict()
                for ab1, ab2 in zip(age_bins[:-1], age_bins[1:]):
                    age = str(ab1) + '-' + str(ab2)
                    if ab1 == 65:
                        age = '65+'  # Combine the last two age groups
                    dd['age'] = [age]
                    dd['sex'] = sex_labels[sex]
                    dd['prevalence'] = sim.results[disease][f'prevalence_{sex}_{ab1}_{ab2}'][-1]
                    dd['new_infections'] = sim.results[disease][f'new_infections_{sex}_{ab1}_{ab2}'][-120:].mean()
                    dd['symp_prevalence'] = sim.results[disease][f'symp_prevalence_{sex}_{ab1}_{ab2}'][-1]
                    dd['disease'] = disease
                    dd['par_idx'] = par_idx
                    dd['scenario'] = scenario
                    dfs += pd.DataFrame(dd)
    epi_df = pd.concat(dfs)
    sc.saveobj(f'{resfolder}/epi_df_{scenario}.df', epi_df)

    # Save SW stats
    sim = [sim for sim in sims if sim.par_idx == 0][0]
    sw_res = sim.results['sw_stats']
    sw_df = sw_res.to_df(resample='year', use_years=True, sep='.')
    sc.saveobj(f'{resfolder}/sw_df_{scenario}.df', sw_df)

    # Save HIV results
    hiv_res = sim.results['total_symptomatic']
    hiv_df = hiv_res.to_df(resample='year', use_years=True, sep='.')
    hiv_df['timevec'] = df.timevec
    sc.saveobj(f'{resfolder}/hiv_df_{scenario}.df', hiv_df)

    return


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1  # 533833
    do_save = True
    do_run = True
    scenario = 'treat80'
    use_calib = True  # Whether to use the calibrated parameters

    # What to run
    to_run = [
        # 'hiv',
        'stis',
    ]

    if 'hiv' in to_run:
        sim = make_sim(add_stis=False, scenario=scenario, seed=seed, debug=debug, start=1990, stop=2041)
        sim.run()
        df = sim.to_df(resample='year', use_years=True, sep='_')
        if do_save: sc.saveobj(f'results/{scenario}_sim.df', df)

        # Process and plot
        df = sc.loadobj(f'results/{scenario}_sim.df')
        df.index = df.timevec
        plot_hiv_sims(df, start_year=1990, which='single')

    if 'stis' in to_run:

        if do_run:
            sim = make_sim(scenario=scenario, use_calib=use_calib, seed=seed, debug=debug, start=1990, stop=2041)
            if use_calib:
                print('Using calibration parameters:')
                print(f'ng_p_symp: {sim.diseases.ng.pars.p_symp}')
                print(f'p_symp_care: {sim.diseases.ct.pars.p_symp_care}')
            sim.run()

            df = sim.to_df(resample='year', use_years=True, sep='_')
            if do_save: sc.saveobj(f'results/{scenario}_sim.df', df)

        else:
            df = sc.loadobj(f'results/{scenario}_sim.df')

        plot_sti_sims(df, start_year=2000, end_year=2040, which='single', fext=scenario)
        plot_sti_tx(df, start_year=2000, fext=scenario, sex='f')



