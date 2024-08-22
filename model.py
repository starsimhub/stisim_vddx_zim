"""
Create a model with HIV plus 4 co-circulating discharging STIs:
    - chlamydia, gonorrhea, trichomoniasis, and other (BV+)
Used for evaluation of etiological tests compared to syndromic management.
"""

# %% Imports and settings
# import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
import pandas as pd
import pylab as pl

from hiv_model import make_hiv, make_hiv_intvs
from analyzers import overtreatment_stats, coinfection_stats


def make_stis():
    gon = sti.Gonorrhea(
        beta_m2f=0.055,
        beta_f2m=0.03,
        init_prev_data=pd.read_csv('data/init_prev_ng.csv'),
        rel_init_prev=0.2
    )
    chlamydia = sti.Chlamydia(
        beta_m2f=0.03,
        beta_f2m=0.015,
        init_prev_data=pd.read_csv('data/init_prev_ct.csv'),
        rel_init_prev=0.8
    )
    trich = sti.Trichomoniasis(
        beta_m2f=0.012,
        beta_f2m=0.006,
        p_clear=[
            ss.bernoulli(p=0.2),
            ss.bernoulli(p=1),  # Men assumed to clear (https://sti.bmj.com/content/76/4/248)
        ],
        init_prev_data=pd.read_csv('data/init_prev_tv.csv'),
        rel_init_prev=2
    )
    vd = sti.DischargingSTI(
        beta_m2f=0.2,
        beta_f2m=0.1,
        init_prev_data=pd.read_csv('data/init_prev_vd.csv'),
    )
    stis = [gon, chlamydia, trich, vd]

    return stis


def make_testing(diseases, start=1980, end=2020):
    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.ti)
        tv_care = sim.diseases.tv.symptomatic & (sim.diseases.tv.ti_seeks_care == sim.ti)
        ct_care = sim.diseases.ct.symptomatic & (sim.diseases.ct.ti_seeks_care == sim.ti)
        vd_care = sim.diseases.vd.symptomatic & (sim.diseases.vd.ti_seeks_care == sim.ti)
        return (ng_care | tv_care | ct_care | vd_care).uids

    ng_tx = sti.GonorrheaTreatment()
    tv_tx = sti.STITreatment(disease='tv', name='tv_tx', label='tv_tx')
    ct_tx = sti.STITreatment(disease='ct', name='ct_tx', label='ct_tx')
    vd_tx = sti.STITreatment(disease='vd', name='vd_tx', label='vd_tx')

    treat_prob = pd.read_csv('data/treat_prob.csv')
    syndromic = sti.SyndromicMgmt(
        treat_prob_data=treat_prob,
        diseases=diseases,
        eligibility=seeking_care_discharge,
        treatments=[ng_tx, tv_tx, ct_tx, vd_tx],
    )
    intvs = [syndromic, ng_tx, tv_tx, ct_tx, vd_tx]
    return intvs


def make_sim(location='zimbabwe', seed=1, n_agents=None, dt=1/12, start=1990, end=2030, debug=False, verbose=0.1):

    total_pop = {1970: 5.203e6, 1980: 7.05e6, 1990: 9980999, 2000: 11.83e6}[start]
    if n_agents is None: n_agents = [int(10e3), int(5e2)][debug]
    if dt is None: dt = [1/12, 1][debug]

    ####################################################################################################################
    # Demographic modules
    ####################################################################################################################
    fertility_rates = {'fertility_rate': pd.read_csv(f'data/{location}_asfr.csv')}
    pregnancy = ss.Pregnancy(pars=fertility_rates)
    death_rates = {'death_rate': pd.read_csv(f'data/{location}_deaths.csv'), 'units': 1}
    death = ss.Deaths(death_rates)

    ####################################################################################################################
    # People and networks
    ####################################################################################################################
    ppl = ss.People(n_agents, age_data=pd.read_csv(f'data/{location}_age_{start}.csv', index_col='age')['value'])
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
        condom_data=pd.read_csv(f'data/{location}_condom_use.csv'),
    )
    maternal = ss.MaternalNet()

    ####################################################################################################################
    # Diseases
    ####################################################################################################################
    stis = make_stis()
    hiv = make_hiv()
    diseases = stis + hiv
    intvs = make_testing(stis, start=start, end=end) + make_hiv_intvs(end=end)
    analyzers = [overtreatment_stats, coinfection_stats]

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

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1

    sim = make_sim(seed=seed, debug=debug, dt=1/12)
    sim.run(verbose=0.1)
    sim.plot('ng')
    pl.show()

    # import sciris as sc
    # si = sc.findfirst(sim.results.yearvec, 2020)
    # ei = sc.findfirst(sim.results.yearvec, 2021)
    # (sim.results.coinfection_stats.ng_only[si:ei].mean()+
    # sim.results.coinfection_stats.ng_ct[si:ei].mean()+
    # sim.results.coinfection_stats.ng_tv[si:ei].mean()+
    # sim.results.coinfection_stats.ng_vd[si:ei].mean()+
    # sim.results.coinfection_stats.ng_ct_tv[si:ei].mean()+
    # sim.results.coinfection_stats.ng_tv_vd[si:ei].mean()+
    # sim.results.coinfection_stats.ng_ct_tv_vd[si:ei].mean())