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
import pylab as pl

from hiv_model import make_hiv, make_hiv_intvs
from interventions import make_testing


def make_stis():
    ng = sti.Gonorrhea(
        beta_m2f=0.055,
        beta_f2m=0.03,
        beta_m2c=0,
        init_prev_data=pd.read_csv('data/init_prev_ng.csv'),
        rel_init_prev=0.2
    )
    ct = sti.Chlamydia(
        beta_m2f=0.03,
        beta_f2m=0.015,
        beta_m2c=0,
        init_prev_data=pd.read_csv('data/init_prev_ct.csv'),
        rel_init_prev=0.8
    )
    tv = sti.Trichomoniasis(
        beta_m2f=0.012,
        beta_f2m=0.006,
        beta_m2c=0,
        p_clear=[
            ss.bernoulli(p=0.2),
            ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
        ],
        init_prev_data=pd.read_csv('data/init_prev_tv.csv'),
        rel_init_prev=2
    )
    bv = sti.DischargingSTI(
        beta_m2f=0.2,
        beta_f2m=0.1,
        beta_m2c=0,
        init_prev_data=pd.read_csv('data/init_prev_bv.csv'),
    )

    return ng, ct, tv, bv


def make_sim(scenario='soc', seed=1, n_agents=None, dt=1/12, start=1990, end=2030, debug=False, verbose=0.1):

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
    ng, ct, tv, bv = make_stis()
    stis = [ng, ct, tv, bv]
    hiv = make_hiv()
    diseases = stis + [hiv]

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    intvs = make_testing(ng, ct, tv, bv, scenario=scenario, end=end) + make_hiv_intvs()
    analyzers = []  #, overtreatment_stats, coinfection_stats]

    sim = ss.Sim(
        dt=dt,
        rand_seed=seed,
        total_pop=total_pop,
        start=start,
        end=end,
        people=ppl,
        diseases=diseases,
        networks=[sexual, maternal],
        demographics=[pregnancy, death],
        interventions=intvs,
        analyzers=analyzers,
        verbose=verbose,
    )

    # Store scenario for grouping
    sim.scenario = scenario

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1

    sim = make_sim(scenario='panel', seed=seed, debug=debug, end=2030)
    sim.run(verbose=0.1)
    sim.plot('ng')
    pl.show()

    # import sciris as sc
    # si = sc.findfirst(sim.results.yearvec, 2020)
    # ei = sc.findfirst(sim.results.yearvec, 2021)
    # (sim.results.syndromicmgmt.ng_only[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv_bv[si:ei].sum())
    #
    # (sim.results.syndromicmgmt.ct_only[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv_bv[si:ei].sum())
    #
    # (sim.results.syndromicmgmt.tv_only[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv_bv[si:ei].sum())
    #
    # (sim.results.syndromicmgmt.bv_only[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_tv_bv[si:ei].sum())
    #
    # # Number with two infections
    # (sim.results.syndromicmgmt.ng_ct[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.tv_bv[si:ei].sum())
    #
    # # Number with three infections
    # (sim.results.syndromicmgmt.ng_ct_tv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_ct_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ng_tv_bv[si:ei].sum()+
    # sim.results.syndromicmgmt.ct_tv_bv[si:ei].sum())
    #
    # # Number new symptomatic
    # (sim.results.ng.new_symptomatic[si:ei].sum()+
    # sim.results.ct.new_symptomatic[si:ei].sum()+
    # sim.results.tv.new_symptomatic[si:ei].sum()+
    # sim.results.bv.new_symptomatic[si:ei].sum())
