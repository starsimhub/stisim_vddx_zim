"""
Custom interventions for the discharge valuation
"""

import stisim as sti
import starsim as ss
import numpy as np
import pandas as pd


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


class SyndromicMgmt(sti.SymptomaticTesting):
    def init_results(self):
        super().init_results()
        self.define_results(
            ss.Result('new_ng_only', dtype=int, label='Only NG'),
            ss.Result('new_ct_only', dtype=int, label='Only CT'),
            ss.Result('new_tv_only', dtype=int, label='Only TV'),
            ss.Result('new_bv_only', dtype=int, label='Only BV'),
            ss.Result('new_ng_ct', dtype=int, label='NG & CT'),
            ss.Result('new_ng_tv', dtype=int, label='NG & TV'),
            ss.Result('new_ng_bv', dtype=int, label='NG & BV'),
            ss.Result('new_ct_tv', dtype=int, label='CT & TV'),
            ss.Result('new_ct_bv', dtype=int, label='CT & BV'),
            ss.Result('new_tv_bv', dtype=int, label='TV & BV'),
            ss.Result('new_ng_ct_tv', dtype=int, label='NG & CT & TV'),
            ss.Result('new_ng_ct_bv', dtype=int, label='NG & CT & BV'),
            ss.Result('new_ng_tv_bv', dtype=int, label='NG & TV & BV'),
            ss.Result('new_ct_tv_bv', dtype=int, label='CT & TV & BV'),
            ss.Result('new_ng_ct_tv_bv', dtype=int, label='NG & CT & TV & BV'),
            ss.Result('new_all_ng', dtype=int, label='All NG'),
            ss.Result('new_all_ct', dtype=int, label='All CT'),
            ss.Result('new_all_tv', dtype=int, label='All TV'),
            ss.Result('new_all_bv', dtype=int, label='All BV'),
            ss.Result('new_sti1', dtype=int, label='1 STI'),
            ss.Result('new_sti2', dtype=int, label='2 STIs'),
            ss.Result('new_sti3', dtype=int, label='3 STIs'),
            ss.Result('new_sti4', dtype=int, label='4 STIs'),

        )

        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti

        # Custom results
        self.results['new_ng_only'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_ct_only'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_tv_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_bv_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)

        self.results['new_ng_ct'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_ng_tv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_ng_bv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['new_ct_tv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_ct_bv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['new_tv_bv'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['new_ng_ct_tv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['new_ng_ct_bv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['new_ng_tv_bv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)
        self.results['new_ct_tv_bv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['new_ng_ct_tv_bv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['new_all_ng'][ti] = len((just_tested & ppl.ng.infected).uids)
        self.results['new_all_ct'][ti] = len((just_tested & ppl.ct.infected).uids)
        self.results['new_all_tv'][ti] = len((just_tested & ppl.tv.infected).uids)
        self.results['new_all_bv'][ti] = len((just_tested & ppl.bv.infected).uids)

        self.results['new_sti1'][ti] = self.results['new_ng_only'][ti] + self.results['new_ct_only'][ti] + self.results['new_tv_only'][ti] + self.results['new_bv_only'][ti]
        self.results['new_sti2'][ti] = (self.results['new_ng_ct'][ti]+
                                        self.results['new_ng_tv'][ti]+
                                        self.results['new_ng_bv'][ti]+
                                        self.results['new_ct_tv'][ti]+
                                        self.results['new_ct_bv'][ti]+
                                        self.results['new_tv_bv'][ti])
        self.results['new_sti3'][ti] = (self.results['new_ng_ct_tv'][ti]+
                                        self.results['new_ng_ct_bv'][ti]+
                                        self.results['new_ng_tv_bv'][ti]+
                                        self.results['new_ct_tv_bv'][ti])

        self.results['new_sti4'][ti] = self.results['new_ng_ct_tv_bv'][ti]

        return


# %%  Algorithms

def make_testing(ng, ct, tv, bv, prop_treat=None, poc=None, stop=2040):

    intv_year = 2028

    # Handle inputs
    synd_end = intv_year if poc else stop
    # Translate prop_treat into sens and spec
    sens = dict(  # Reflective of a treat-all approach
            ng=[prop_treat, prop_treat],
            ct=[prop_treat, prop_treat],
            tv=[prop_treat, prop_treat],
            bv=[prop_treat],
    )
    spec = dict(
            ng=[1-prop_treat, 1-prop_treat],
            ct=[1-prop_treat, 1-prop_treat],
            tv=[1-prop_treat, 1-prop_treat],
            bv=[1-prop_treat],
    )

    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.diseases.ng.ti)
        tv_care = sim.diseases.tv.symptomatic & (sim.diseases.tv.ti_seeks_care == sim.diseases.tv.ti)
        ct_care = sim.diseases.ct.symptomatic & (sim.diseases.ct.ti_seeks_care == sim.diseases.ct.ti)
        bv_care = sim.diseases.bv.symptomatic & (sim.diseases.bv.ti_seeks_care == sim.diseases.bv.ti)
        return (ng_care | ct_care | tv_care | bv_care).uids

    ng_tx = sti.GonorrheaTreatment(
        name='ng_tx',
        rel_treat_unsucc=0.005,
        rel_treat_unneed=0.0005,
    )
    ct_tx = sti.STITreatment(diseases='ct', name='ct_tx', label='ct_tx')
    metronidazole = sti.STITreatment(diseases=['tv', 'bv'], name='metronidazole', label='metronidazole')
    treatments = [ng_tx, ct_tx, metronidazole]
    disease_treatment_map = {
        'ng': ng_tx, 'ct': ct_tx, 'tv': metronidazole, 'bv': metronidazole
    }

    # Syndromic management
    syndromic = SyndromicMgmt(
        sens=sens,
        spec=spec,
        stop=synd_end,
        diseases=[ng, ct, tv, bv],
        eligibility=seeking_care_discharge,
        treatments=treatments,
        disease_treatment_map=disease_treatment_map,
    )

    if not poc:
        intvs = [syndromic, ng_tx, ct_tx, metronidazole]

    if poc:
        panel = SyndromicMgmt(
            sens = dict(
                ng=[0.95, 0.95],
                ct=[0.95, 0.95],
                tv=[0.95, 0.95],
                bv=[prop_treat],
            ),
            spec = dict(
                ng=[0.95, 0.95],
                ct=[0.95, 0.95],
                tv=[0.95, 0.95],
                bv=[1-prop_treat],
            ),
            start=intv_year,
            diseases=[ng, ct, tv, bv],
            eligibility=seeking_care_discharge,
            treatments=treatments,
            disease_treatment_map=disease_treatment_map,
        )
        intvs = [syndromic, panel, ng_tx, ct_tx, metronidazole]

    return intvs
