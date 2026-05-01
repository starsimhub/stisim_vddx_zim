"""
Custom interventions for the discharge valuation
"""

import stisim as sti
import starsim as ss
import numpy as np
import pandas as pd
import sciris as sc


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


class SyndromicMgmt(sti.SyndromicMgmt):
    """
    Project-local extension of sti.SyndromicMgmt that adds co-infection
    burden results (new_sti1/2/3/4 per sex) needed for the VDDX-Zim analysis.
    All core logic lives in the base class.
    """

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        for sk in ['', 'f', 'm']:
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' - {sk.upper()}'
            results += [
                ss.Result('new_sti1' + skk, dtype=int, label='1 STI' + skl),
                ss.Result('new_sti2' + skk, dtype=int, label='2 STIs' + skl),
                ss.Result('new_sti3' + skk, dtype=int, label='3 STIs' + skl),
                ss.Result('new_sti4' + skk, dtype=int, label='4 STIs' + skl),
            ]
        self.define_results(*results)
        return

    def store_results(self):
        super().store_results()
        ti  = self.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti

        sexdict = {'': 'alive', 'f': 'female', 'm': 'male'}
        for sk, sl in sexdict.items():
            skk = '' if sk == '' else f'_{sk}'
            s = just_tested & ppl[sl]
            ng, ct, tv, bv = ppl.ng.infected, ppl.ct.infected, ppl.tv.infected, ppl.bv.infected
            self.results['new_sti1' + skk][ti] = count(s & ( ng & ~ct & ~tv & ~bv
                                                            | ~ng &  ct & ~tv & ~bv
                                                            | ~ng & ~ct &  tv & ~bv
                                                            | ~ng & ~ct & ~tv &  bv))
            self.results['new_sti2' + skk][ti] = count(s & ( ng &  ct & ~tv & ~bv
                                                            | ng & ~ct &  tv & ~bv
                                                            | ng & ~ct & ~tv &  bv
                                                            | ~ng &  ct &  tv & ~bv
                                                            | ~ng &  ct & ~tv &  bv
                                                            | ~ng & ~ct &  tv &  bv))
            self.results['new_sti3' + skk][ti] = count(s & ( ng &  ct &  tv & ~bv
                                                            | ng &  ct & ~tv &  bv
                                                            | ng & ~ct &  tv &  bv
                                                            | ~ng &  ct &  tv &  bv))
            self.results['new_sti4' + skk][ti] = count(s & ng & ct & tv & bv)
        return


class AMR(ss.Intervention):
    """ Sketch of AMR intervention, not in use and not tested """
    def __init__(self, amr_scen='baseline', **kwargs):
        super().__init__(**kwargs)

        self.years = [2027, 2030, 2035, 2040]
        self.ng_tx_eff = None
        if amr_scen == 'baseline':
            self.ng_tx_eff_pars = [0.96, 0.96, 0.96, 0.96]
        elif amr_scen == 'moderate':
            self.ng_tx_eff_pars = [0.96, 0.93, 0.90, 0.87]
        elif amr_scen == 'aggressive':
            self.ng_tx_eff_pars = [0.96, 0.9, 0.84, 0.78]
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.ng_tx_eff = sc.smoothinterp(self.t.yearvec, self.years, self.ng_tx_eff_pars, smoothness=0)
        return

    def step(self):
        self.sim.interventions.ng_tx.pars.base_treat_eff.set(self.ng_tx_eff[self.ti])
        return


# %%  Algorithms
def make_tx_mix(scenario):
    if 'treat100' in scenario:
        tx_mix_cerv = dict(
            all3=[1.00, 0.10],
            ngct=[0.00, 0.80],
            mtnz=[0.00, 0.00],
            none=[0.00, 0.10],
        )
        tx_mix_noncerv = sc.dcp(tx_mix_cerv)
    elif 'treat80' in scenario:
        tx_mix_cerv = dict(
            all3=[0.50, 0.10],
            ngct=[0.20, 0.80],
            mtnz=[0.15, 0.00],
            none=[0.15, 0.10],
        )
        tx_mix_noncerv = dict(
            all3=[0.40, 0.10],
            ngct=[0.10, 0.80],
            mtnz=[0.25, 0.00],
            none=[0.25, 0.10],
        )
    elif 'treat50' in scenario:
        tx_mix_cerv = dict(
            all3=[0.25, 0.10],
            ngct=[0.25, 0.80],
            mtnz=[0.25, 0.00],
            none=[0.25, 0.10],
        )
        tx_mix_noncerv = sc.dcp(tx_mix_cerv)
    elif 'treat30' in scenario:
        tx_mix_cerv = dict(
            all3=[0.10, 0.10],
            ngct=[0.20, 0.80],
            mtnz=[0.10, 0.00],
            none=[0.60, 0.10],
        )
        tx_mix_noncerv = sc.dcp(tx_mix_cerv)
    return tx_mix_cerv, tx_mix_noncerv


def neg_panel_mix(scenario):
    if 'treat100' in scenario:
        p_mtnz = 1.0
    elif 'treat80' in scenario:
        p_mtnz = 0.8
    elif 'treat50' in scenario:
        p_mtnz = 0.5
    elif 'treat30' in scenario:
        p_mtnz = 0.3
    return p_mtnz


def make_testing(ng, ct, tv, bv, scenario=None, poc=None, stop=2040):

    intv_year = 2027

    # Handle inputs
    synd_end = intv_year if poc else stop

    # Testing interventions
    def seeking_care_vds(sim):
        dis = sim.diseases
        female = sim.people.female
        ng_care = dis.ng.symptomatic & (dis.ng.ti_seeks_care == dis.ng.ti) & female
        tv_care = dis.tv.symptomatic & (dis.tv.ti_seeks_care == dis.tv.ti) & female
        ct_care = dis.ct.symptomatic & (dis.ct.ti_seeks_care == dis.ct.ti) & female
        bv_care = dis.bv.symptomatic & (dis.bv.ti_seeks_care == dis.bv.ti) & female
        return (ng_care | ct_care | tv_care | bv_care).uids

    def seeking_care_uds(sim):
        dis = sim.diseases
        male = sim.people.male
        ng_care = dis.ng.symptomatic & (dis.ng.ti_seeks_care == dis.ng.ti) & male
        tv_care = dis.tv.symptomatic & (dis.tv.ti_seeks_care == dis.tv.ti) & male
        ct_care = dis.ct.symptomatic & (dis.ct.ti_seeks_care == dis.ct.ti) & male
        return (ng_care | ct_care | tv_care).uids

    ng_tx = sti.GonorrheaTreatment(
        name='ng_tx',
        rel_treat_unsucc=0.005,
        rel_treat_unneed=0.0005,
    )
    ct_tx = sti.STITreatment(diseases='ct', name='ct_tx', label='ct_tx')
    metronidazole = sti.STITreatment(diseases=['tv', 'bv'], name='metronidazole', label='metronidazole')
    treatments = [ng_tx, ct_tx, metronidazole]
    outcome_treatment_map = dict(
        all3=treatments,
        ngct=[ng_tx, ct_tx],
        mtnz=[metronidazole],
        none=[],
    )
    tx_mix_cerv, tx_mix_noncerv = make_tx_mix(scenario)

    # Syndromic management of VDS
    syndromic_vds = SyndromicMgmt(
        name='syndromic_vds',
        label='syndromic_vds',
        tx_mix_cerv=tx_mix_cerv,
        tx_mix_noncerv=tx_mix_noncerv,
        stop=synd_end,
        diseases=[ng, ct, tv, bv],
        eligibility=seeking_care_vds,
        treatments=treatments,
        outcome_treatment_map=outcome_treatment_map,
    )

    syndromic_uds = SyndromicMgmt(
        name='syndromic_uds',
        label='syndromic_uds',
        tx_mix_cerv=tx_mix_cerv,
        tx_mix_noncerv=tx_mix_noncerv,
        stop=stop,
        diseases=[ng, ct, tv],
        eligibility=seeking_care_uds,
        treatments=treatments,
        outcome_treatment_map=outcome_treatment_map,
    )

    if not poc:
        intvs = [syndromic_vds, syndromic_uds, ng_tx, ct_tx, metronidazole]

    if poc:
        disease_treatment_map = {'ng': ng_tx, 'ct': ct_tx, 'tv': metronidazole}
        p_mtnz = neg_panel_mix(scenario)

        panel = sti.SymptomaticTesting(
            name='panel',
            label='panel',
            start=intv_year,
            diseases=[ng, ct, tv],
            eligibility=seeking_care_vds,
            treatments=treatments,
            disease_treatment_map=disease_treatment_map,
            p_mtnz=p_mtnz,
            negative_treatments=[metronidazole],
        )
        intvs = [syndromic_vds, syndromic_uds, panel, ng_tx, ct_tx, metronidazole]

    return intvs
