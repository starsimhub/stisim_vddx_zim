"""
Gonorrhea model - mostly used for testing
"""

# %% Imports and settings
import sciris as sc
import numpy as np
import starsim as ss
import stisim as sti
import pandas as pd
import pylab as pl

from hiv_model import make_hiv, make_hiv_intvs
from utils import unneeded_results, set_font
from analyzers import total_symptomatic


def make_sim(scenario='soc', seed=1, n_agents=None, dt=1/12, start=1980, end=2030, debug=False, verbose=0.1):

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
    ng = sti.Gonorrhea(
        beta_m2f=0.5,
        beta_m2c=0,
        init_prev_data=pd.read_csv('data/init_prev_ng.csv'),
        rel_init_prev=0.2
    )
    hiv = make_hiv()
    diseases = [ng, hiv]

    ####################################################################################################################
    # Interventions and analyzers
    ####################################################################################################################
    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.ti)
        return (ng_care).uids

    ng_tx = sti.GonorrheaTreatment(
        name='ng_tx',
        rel_treat_unsucc=0.001,
        rel_treat_unneed=0.0001,
    )
    disease_treatment_map = {'ng': ng_tx}

    syndromic = sti.SymptomaticTesting(
        diseases=[ng],
        eligibility=seeking_care_discharge,
        treatments=[ng_tx],
        disease_treatment_map=disease_treatment_map,
    )

    intvs = [syndromic, ng_tx] + make_hiv_intvs()

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
        connectors=[sti.hiv_ng(hiv, ng)],
        verbose=verbose,
    )

    # Store scenario for grouping
    sim.scenario = scenario

    return sim


if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_save = True

    if True:
        sim = make_sim(scenario='soc', seed=seed, debug=debug, start=1990, end=2020)
        sim.run(verbose=0.1)
        df = sti.finalize_results(sim, modules_to_drop=unneeded_results)
        if do_save: sc.saveobj('results/sim_ng.df', df)

    # Process and plot
    df = sc.loadobj('results/sim_ng.df')
    from plot_sims import plot_ng_sim
    plot_ng_sim(df, start_year=1990, end_year=2020)


