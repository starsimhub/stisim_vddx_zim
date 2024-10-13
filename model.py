"""
Create a model with HIV plus 4 co-circulating discharging STIs:
    - chlamydia, gonorrhea, trichomoniasis, and other (BV+)
Used for evaluation of etiological tests compared to syndromic management.
"""

# %% Imports and settings
import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
import pandas as pd

from hiv_model import make_hiv, make_hiv_intvs
from interventions import make_testing
from utils import unneeded_results
from analyzers import total_symptomatic


def make_stis(bv_beta_m2f=0.2):
    ng = sti.Gonorrhea(
        beta_m2f=0.07,
        init_prev_data=pd.read_csv('data/init_prev_ng.csv'),
        rel_init_prev=0.2
    )
    ct = sti.Chlamydia(
        beta_m2f=0.07,
        init_prev_data=pd.read_csv('data/init_prev_ct.csv'),
        rel_init_prev=1.5
    )
    tv = sti.Trichomoniasis(
        beta_m2f=0.15,
        p_clear=[
            ss.bernoulli(p=0.1),
            ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
        ],
        init_prev_data=pd.read_csv('data/init_prev_tv.csv'),
        rel_init_prev=10
    )
    bv = sti.DischargingSTI(
        beta_m2f=bv_beta_m2f,
        init_prev_data=pd.read_csv('data/init_prev_bv.csv'),
    )

    return ng, ct, tv, bv


def make_sim(scenario='soc', seed=1, n_agents=None, bv_beta_m2f=0.15, dt=1/12, start=1980, stop=2030, debug=False, verbose=0.1):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(5e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_rates = {'fertility_rate': pd.read_csv(f'data/asfr.csv')}
    pregnancy = ss.Pregnancy(pars=fertility_rates)
    death_rates = {'death_rate': pd.read_csv(f'data/deaths.csv'), 'units': 1}
    death = ss.Deaths(death_rates)

    ####################################################################################################################
    # People and networks
    ####################################################################################################################
    ppl = ss.People(n_agents, age_data=pd.read_csv(f'data/age_dist_{start}.csv', index_col='age')['value'])
    sexual = sti.FastStructuredSexual(
        acts=ss.lognorm_ex(80, 30),
        prop_f1=0.2,
        prop_f2=0.05,
        prop_m1=0.2,
        f1_conc=0.05,
        f2_conc=0.25,
        m1_conc=0.15,
        m2_conc=0.3,
        p_pair_form=0.8,  # 0.6,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )
    maternal = ss.MaternalNet()

    ####################################################################################################################
    # Diseases
    ####################################################################################################################
    ng, ct, tv, bv = make_stis(bv_beta_m2f=bv_beta_m2f)
    stis = [ng, ct, tv, bv]
    hiv = make_hiv()
    diseases = stis + [hiv]

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    intvs = make_testing(ng, ct, tv, bv, scenario=scenario, stop=stop) + make_hiv_intvs()
    analyzers = [total_symptomatic()]  #, overtreatment_stats, coinfection_stats]

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
        connectors=[sti.hiv_ng(hiv, ng), sti.hiv_ct(hiv, ct), sti.hiv_tv(hiv, tv)],
        verbose=verbose,
    )

    # Store scenario and background BV rate for grouping
    sim.scenario = scenario
    sim.bv_beta_m2f = bv_beta_m2f

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_save = True

    if True:
        sim = make_sim(scenario='soc', seed=seed, debug=debug, start=1980, stop=2030)
        sim.run(verbose=0.1)
        df = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        if do_save: sc.saveobj('results/sim.df', df)
        # df['ng.rel_treat'].to_csv('results/Ciprofloxacin.csv')

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

    # Process and plot
    from plot_sims import *
    df = sc.loadobj('results/sim1.df')
    plot_sti_sims(df, start_year=2000, end_year=2040, which='single', fext='_alt')
    plot_sti_tx(df, start_year=2000, fext='_alt')
    # plot_hiv_sims(df, start_year=2000, which='single', fext='_alt')
    # plot_ng_sim(df, start_year=1990)

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


