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
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_ng_only', npts, dtype=int, scale=True, label='Only NG'),
            ss.Result(self.name, 'new_ct_only', npts, dtype=int, scale=True, label='Only CT'),
            ss.Result(self.name, 'new_tv_only', npts, dtype=int, scale=True, label='Only TV'),
            ss.Result(self.name, 'new_bv_only', npts, dtype=int, scale=True, label='Only BV'),
            ss.Result(self.name, 'new_ng_ct', npts, dtype=int, scale=True, label='NG & CT'),
            ss.Result(self.name, 'new_ng_tv', npts, dtype=int, scale=True, label='NG & TV'),
            ss.Result(self.name, 'new_ng_bv', npts, dtype=int, scale=True, label='NG & BV'),
            ss.Result(self.name, 'new_ct_tv', npts, dtype=int, scale=True, label='CT & TV'),
            ss.Result(self.name, 'new_ct_bv', npts, dtype=int, scale=True, label='CT & BV'),
            ss.Result(self.name, 'new_tv_bv', npts, dtype=int, scale=True, label='TV & BV'),
            ss.Result(self.name, 'new_ng_ct_tv', npts, dtype=int, scale=True, label='NG & CT & TV'),
            ss.Result(self.name, 'new_ng_ct_bv', npts, dtype=int, scale=True, label='NG & CT & BV'),
            ss.Result(self.name, 'new_ng_tv_bv', npts, dtype=int, scale=True, label='NG & TV & BV'),
            ss.Result(self.name, 'new_ct_tv_bv', npts, dtype=int, scale=True, label='CT & TV & BV'),
            ss.Result(self.name, 'new_ng_ct_tv_bv', npts, dtype=int, scale=True, label='NG & CT & TV & BV'),
            ss.Result(self.name, 'new_all_ng', npts, dtype=int, scale=True, label='All NG'),
            ss.Result(self.name, 'new_all_ct', npts, dtype=int, scale=True, label='All CT'),
            ss.Result(self.name, 'new_all_tv', npts, dtype=int, scale=True, label='All TV'),
            ss.Result(self.name, 'new_all_bv', npts, dtype=int, scale=True, label='All BV'),
            ss.Result(self.name, 'new_sti1', npts, dtype=int, scale=True, label='1 STI'),
            ss.Result(self.name, 'new_sti2', npts, dtype=int, scale=True, label='2 STIs'),
            ss.Result(self.name, 'new_sti3', npts, dtype=int, scale=True, label='3 STIs'),
            ss.Result(self.name, 'new_sti4', npts, dtype=int, scale=True, label='4 STIs'),
        ]
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


class Panel(sti.SymptomaticTesting):
    def __init__(self, pars=None, treatments=None, diseases=None, disease_treatment_map=None, years=None, start=None, end=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(treatments=treatments, diseases=diseases, disease_treatment_map=disease_treatment_map, years=years, start=start, end=end, eligibility=eligibility, name=name, label=label, **kwargs)
        self.default_pars(
            sens=dict(
                ng=ss.bernoulli(0.95),
                ct=ss.bernoulli(0.95),
                tv=ss.bernoulli(0.95),
                bv=ss.bernoulli(0.56),
            ),
            spec=dict(
                ng=ss.bernoulli(0.9),
                ct=ss.bernoulli(0.9),
                tv=ss.bernoulli(0.9),
                bv=ss.bernoulli(0.71),
            ),
            dt_scale=False,
        )
        self.update_pars(pars, **kwargs)
        return
    #
    # def apply(self, sim, uids=None):
    #     """ Apply the testing intervention """
    #     # Apply if within the start years
    #     if (sim.year >= self.start) & (sim.year < self.end):
    #
    #         super().apply(sim, uids)
    #         dismissed_uids = self.ti_dismissed == sim.ti
    #         # Give BV treatment to anyone where nothing was detected
    #         # Question, should this be for women only? Retain symptomatic care for men. Reduce symptoms for BV for men
    #         self.sim.interventions.bv_tx.eligibility = dismissed_uids
    #
    #         # Update results
    #         self.update_results()
    #
    #     return


# %%  Algorithms

def make_testing(ng, ct, tv, bv, scenario='soc', end=2040):

    intv_year = 2000

    if scenario == 'soc':
        synd_end = end
    else:
        synd_end = intv_year

    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.ti)
        tv_care = sim.diseases.tv.symptomatic & (sim.diseases.tv.ti_seeks_care == sim.ti)
        ct_care = sim.diseases.ct.symptomatic & (sim.diseases.ct.ti_seeks_care == sim.ti)
        bv_care = sim.diseases.bv.symptomatic & (sim.diseases.bv.ti_seeks_care == sim.ti)
        return (ng_care | ct_care | tv_care | bv_care).uids

    ng_tx = sti.GonorrheaTreatment(
        name='ng_tx',
        rel_treat_unsucc=0.05,
        rel_treat_unneed=0.01,
    )
    ct_tx = sti.STITreatment(disease='ct', name='ct_tx', label='ct_tx')
    metronidazole = sti.STITreatment(disease=['tv', 'bv'], name='metronidazole', label='metronidazole')
    treatments = [ng_tx, ct_tx, metronidazole]
    disease_treatment_map = {
        'ng': ng_tx, 'ct': ct_tx, 'tv': metronidazole, 'bv': metronidazole
    }

    syndromic = SyndromicMgmt(
        end=synd_end,
        diseases=[ng, ct, tv, bv],
        eligibility=seeking_care_discharge,
        treatments=treatments,
        disease_treatment_map=disease_treatment_map,
    )

    intvs = [syndromic, ng_tx, ct_tx, metronidazole]

    if scenario == 'panel':
        panel = Panel(
            start=intv_year,
            eligibility=seeking_care_discharge,
            diseases=[ng, ct, tv, bv],
            treatments=treatments,
            disease_treatment_map=disease_treatment_map,
        )

        intvs += [panel]

    return intvs
