"""
Create a model with HIV plus 4 co-circulating discharging STIs:
    - chlamydia, gonorrhea, trichomoniasis, and other (BV+)
Used for evaluation of etiological tests compared to syndromic management.
"""

# %% Imports and settings
import starsim as ss
import stisim as sti

from hiv_model import make_hiv, make_hiv_intvs
from interventions import make_testing
from plot_sims import *


def make_stis(p_symp=None, p_symp_care=None, ng=None, ct=None, tv=None):
    ng = sti.Gonorrhea(
        beta_m2f=ng['beta_m2f'],
        eff_condom=ng['eff_condom'],
        p_symp=[p_symp['ng'], 0.65],
        p_symp_care=[p_symp_care['ng'], 0.83]
    )
    ct = sti.Chlamydia(
        beta_m2f=ct['beta_m2f'],
        eff_condom=ct['eff_condom'],
        p_symp=[p_symp['ct'], 0.54],
        p_symp_care=[p_symp_care['ct'], 0.83]
    )
    tv = sti.Trichomoniasis(
        beta_m2f=tv['beta_m2f'],
        eff_condom=tv['eff_condom'],
        p_symp=[p_symp['tv'], 0.5],
        p_symp_care=[p_symp_care['tv'], 0.27]
    )
    bv = sti.BV()

    return ng, ct, tv, bv


def make_sim(seed=1, n_agents=None, dt=1/12, start=1990, stop=2030, debug=False, verbose=1/12, add_stis=True,
             scenario='treat100', p_symp=None, p_symp_care=None, poc=False, stipars=None, analyzers=None):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1985: 8.691e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(5e3), int(5e2)][debug]
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
    sexual = sti.FastStructuredSexual(
        prop_f0=0.8,
        prop_f2=0.05,
        prop_m0=0.65,
        f1_conc=0.05,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.3,
        p_pair_form=0.6,  # 0.6,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet(unit='month')

    ####################################################################################################################
    # Diseases
    ####################################################################################################################
    hiv = make_hiv()
    diseases = [hiv]
    if add_stis:
        ng, ct, tv, bv = make_stis(p_symp=p_symp, p_symp_care=p_symp_care, **stipars)
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
    # analyzers = [total_symptomatic()]  #, overtreatment_stats, coinfection_stats]

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

    # Store scenario and background BV rate for grouping
    sim.scenario = scenario

    return sim

# Define scenarios
def make_scens():
    scendict = sc.objdict(
        treat100=sc.objdict(
            p_symp=dict(ng=0.1, ct=0.2, tv=0.3),
            p_symp_care=dict(ng=0.75, ct=0.75, tv=0.6),
            stipars = dict(
                ng=dict(beta_m2f=0.2, eff_condom=0.8),
                ct=dict(beta_m2f=0.2, eff_condom=0.8),
                tv=dict(beta_m2f=0.1, eff_condom=0.8),
            ),
            poc=False,
        )
    )

    scendict['treat80'] = sc.dcp(scendict['treat100'])
    scendict['treat80'].p_symp = dict(ng=0.15, ct=0.3, tv=0.45)
    scendict['treat80'].p_symp_care = dict(ng=0.625, ct=0.625, tv=0.5)
    scendict['treat80'].stipars = dict(
        ng=dict(beta_m2f=0.2, eff_condom=0.8),
        ct=dict(beta_m2f=0.2, eff_condom=0.8),
        tv=dict(beta_m2f=0.1, eff_condom=0.8),
    )

    scendict['treat50'] = sc.dcp(scendict['treat100'])
    scendict['treat50'].p_symp = dict(ng=0.2, ct=0.4, tv=0.6)
    scendict['treat50'].stipars = dict(
        ng=dict(beta_m2f=0.2, eff_condom=0.8),
        ct=dict(beta_m2f=0.2, eff_condom=0.8),
        tv=dict(beta_m2f=0.1, eff_condom=0.8),
    )

    for scenario in scendict.keys():
        scendict[scenario+'poc'] = sc.dcp(scendict[scenario])
        scendict[scenario+'poc'].poc = True

    return scendict


def make_scenpars(scenario):
    scendict = make_scens()
    return scendict[scenario]


def load_calib_pars(scenario=None, calib=None, i=0):
    scenpars = make_scenpars(scenario)
    raw_calib_pars = calib.df.iloc[i].to_dict()

    # Overwrite
    diseases = ['ng', 'ct', 'tv']
    for disease in diseases:
        scenpars['p_symp'][disease] = raw_calib_pars[f'{disease}_p_symp']
        scenpars['p_symp_care'][disease] = raw_calib_pars['p_symp_care']
        scenpars['stipars'][disease]['beta_m2f'] = raw_calib_pars[f'{disease}_beta_m2f']
    return scenpars


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
        # 'plot_hiv'
    ]

    if 'hiv' in to_run:
        sim = make_sim(add_stis=False, scenario=scenario, seed=seed, debug=debug, start=1990, stop=2041)
        sim.run()
        df = sim.to_df(resample='year', use_years=True, sep='.')  # Use dots to separate columns
        if do_save: sc.saveobj(f'results/{scenario}_sim.df', df)

        # Process and plot
        df = sc.loadobj(f'results/{scenario}_sim.df')
        plot_hiv_sims(df, start_year=1990, which='single')

    if 'stis' in to_run:
        if use_calib:
            calib = sc.loadobj(f'results/zim_sti_calib_{scenario}.obj')
            scenpars = load_calib_pars(scenario=scenario, calib=calib, i=0)
        else:
            scenpars = make_scenpars(scenario)

        # Add analyzer
        analyzers = [sti.sw_stats(diseases=['ng', 'ct', 'tv'])]
        sim = make_sim(scenario=scenario, **scenpars, analyzers=analyzers, seed=seed, debug=debug, start=1990, stop=2041)
        sim.run()
        df = sim.to_df(resample='year', use_years=True, sep='.')
        if do_save: sc.saveobj(f'results/{scenario}_sim.df', df)

        # Process and plot
        df = sc.loadobj(f'results/{scenario}_sim.df')
        plot_hiv_sims(df, start_year=1990, which='single')
        plot_sti_sims(df, start_year=2000, end_year=2040, which='single', fext=scenario)
        plot_sti_tx(df, start_year=2000, fext=scenario, sex='f')

        # Save age/sex epi results
        dfs = sc.autolist()
        age_bins = sim.diseases.ng.age_bins
        for disease in ['ng', 'ct', 'tv']:
            for sex in ['f', 'm']:
                dd = dict()
                for ab1, ab2 in zip(age_bins[:-1], age_bins[1:]):
                    dd['age'] = [ab1]
                    dd['sex'] = sex
                    dd['prevalence'] = sim.results[disease][f'prevalence_{sex}_{ab1}_{ab2}'][-1]
                    dd['symp_prevalence'] = sim.results[disease][f'symp_prevalence_{sex}_{ab1}_{ab2}'][-1]
                    dd['disease'] = disease
                    dfs += pd.DataFrame(dd)
        epi_df = pd.concat(dfs)
        if do_save: sc.saveobj('results/epi_df.df', epi_df)

        # Save SW stats
        sw_res = sim.results['sw_stats']
        sw_df = sw_res.to_df(resample='year', use_years=True, sep='.')
        if do_save: sc.saveobj(f'results/{scenario}_sw.df', sw_df)


    if 'plot_hiv' in to_run:
        from utils import set_font
        import pylab as pl
        import seaborn as sns
        
        set_font(size=16)
        fig, ax = pl.subplots(1, 1, figsize=(10, 4))
        colors = ['#ee7989', '#ee7989', '#4682b4', '#4682b4']
        linestyles = ['--', '-', '--', '-']
        rdict = {'symp_prev_no_hiv_f': 'HIV- F', 'symp_prev_has_hiv_f': 'HIV+ F', 'symp_prev_no_hiv_m': 'HIV- M', 'symp_prev_has_hiv_m': 'HIV+ M'}

        cn = 0
        bi = 20*12
        for rname, rlabel in rdict.items():
            x = sim.timevec[bi:]
            y = pd.Series(sim.results.total_symptomatic[rname][bi:])
            y = y.rolling(10, min_periods=1).mean()
            ax.plot(x, y*100, color=colors[cn], ls=linestyles[cn], label=rlabel)
            ax.legend()
            ax.set_ylabel('')
            ax.set_xlabel('')
            ax.set_title('Prevalence of discharge (%)')
            cn += 1

        sc.figlayout()
        sc.savefig("figures/epi_hiv.png", dpi=100)



