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


def make_stis(p_symp=None, p_symp_care=None):
    ng = sti.Gonorrhea(
        beta_m2f=0.07,
        p_symp=[p_symp['ng'], 0.65],
        p_symp_care=[p_symp_care['ng'], 0.83]
    )
    ct = sti.Chlamydia(
        beta_m2f=0.06,
        p_symp=[p_symp['ct'], 0.54],
        p_symp_care=[p_symp_care['ct'], 0.83]
    )
    tv = sti.Trichomoniasis(
        beta_m2f=0.1,
        p_symp=[p_symp['tv'], 0.5],
        p_symp_care=[p_symp_care['tv'], 0.27]
    )
    bv = sti.BV()

    return ng, ct, tv, bv


def make_sim(seed=1, n_agents=None, dt=1/12, start=1980, stop=2030, debug=False, verbose=1/12, add_stis=True,
             scenario='treat100', p_symp=None, p_symp_care=None, prop_treat=None, poc=False):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1985: 8.691e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(5e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_data = pd.read_csv(f'data/asfr.csv')
    pregnancy = ss.Pregnancy(fertility_rate=fertility_data)
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
        ng, ct, tv, bv = make_stis(p_symp=p_symp, p_symp_care=p_symp_care)
        stis = [ng, ct, tv, bv]
        diseases += stis  # Add the STIs to the list of diseases

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    intvs = make_hiv_intvs()
    if add_stis:
        intvs += make_testing(ng, ct, tv, bv, prop_treat=prop_treat, poc=poc, stop=stop)
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
        # analyzers=analyzers,
        connectors=connectors,
        verbose=verbose,
    )

    # Store scenario and background BV rate for grouping
    sim.scenario = scenario

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_save = True
    do_run = True
    scenario = 'treat100'

    # Define scenarios
    scendict = sc.objdict(
        treat100=dict(prop_treat=1, p_symp=0.1, p_symp_care=0.75, poc=False),
        treat90=dict(prop_treat=0.9, p_symp=0.1, p_symp_care=5/6, poc=False),
        treat50=dict(prop_treat=0.5, p_symp=0.2, p_symp_care=0.75, poc=False),
        treat100poc=dict(prop_treat=1, p_symp=0.1, p_symp_care=0.75, poc=True),
        treat90poc=dict(prop_treat=0.9, p_symp=0.1, p_symp_care=5/6, poc=True),
        treat50poc=dict(prop_treat=0.5, p_symp=0.2, p_symp_care=0.75, poc=True),
    )
    scenpars = scendict[scenario]

    # What to run
    to_run = [
        # 'hiv',
        'stis',
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
        sim = make_sim(scenario=scenario, scenpars=scenpars, seed=seed, debug=debug, start=1990, stop=2041)
        sim.run()
        df = sim.to_df(resample='year', use_years=True, sep='.')
        if do_save: sc.saveobj(f'results/{scenario}_sim.df', df)

        # Process and plot
        df = sc.loadobj(f'results/{scenario}_sim.df')
        plot_hiv_sims(df, start_year=1990, which='single')
        plot_sti_sims(df, start_year=1990, end_year=2040, which='single')
        plot_sti_tx(df, start_year=1990)

    # from utils import set_font
    # import pylab as pl
    # import seaborn as sns
    #
    # epi_df = sc.loadobj('results/epi_df.df')
    # set_font(size=20)
    # fig, axes = pl.subplots(1, 3, figsize=(10, 4))
    # axes = axes.ravel()
    # colors = ['#ee7989', '#4682b4']
    #

    # # Save age/sex epi results
        # dfs = sc.autolist()
        # for disease in ['ng', 'ct', 'tv']:
        #     for sex in ['female', 'male']:
        #         dd = dict()
        #         dd['age'] = sim.diseases[disease].age_bins[:-1]
        #         dd['prevalence'] = sim.diseases[disease].age_sex_results['prevalence'][sex][:,-1]
        #         dd['symp_prevalence'] = sim.diseases[disease].age_sex_results['symp_prevalence'][sex][:,-1]
        #         dd['disease'] = disease
        #         dd['sex'] = sex
        #         dfs += pd.DataFrame(dd)
        # epi_df = pd.concat(dfs)
        # if do_save: sc.saveobj('results/epi_df.df', epi_df)

    # for pn, disease in enumerate(['ng', 'ct', 'tv']):
    #     ax = axes[pn]
    #     thisdf = epi_df.loc[(epi_df.disease == disease) & (epi_df.age > 1)]
    #     sns.barplot(data=thisdf, x="age", y="symp_prevalence", hue="sex", ax=ax, palette=colors)
    #     ax.set_title(disease.upper())
    #     ax.set_ylabel('')
    #     ax.set_xlabel('')
    #     # ax.legend_.remove()
    #
    # sc.figlayout()
    # sc.savefig("figures/epi.png", dpi=100)

    # set_font(size=16)
    # fig, ax = pl.subplots(1, 1, figsize=(10, 4))
    # colors = ['#ee7989', '#ee7989', '#4682b4', '#4682b4']
    # linestyles = ['--', '-', '--', '-']
    # rdict = {'symp_prev_no_hiv_f': 'HIV- F', 'symp_prev_has_hiv_f': 'HIV+ F', 'symp_prev_no_hiv_m': 'HIV- M', 'symp_prev_has_hiv_m': 'HIV+ M'}

    # cn = 0
    # bi = 20*12
    # for rname, rlabel in rdict.items():
    #     x = sim.timevec[bi:]
    #     y = pd.Series(sim.results.total_symptomatic[rname][bi:])
    #     y = y.rolling(10, min_periods=1).mean()
    #     ax.plot(x, y*100, color=colors[cn], ls=linestyles[cn], label=rlabel)
    #     ax.legend()
    #     ax.set_ylabel('')
    #     ax.set_xlabel('')
    #     ax.set_title('Prevalence of discharge (%)')
    #     cn += 1
    #
    # sc.figlayout()
    # sc.savefig("figures/epi_hiv.png", dpi=100)
    #

