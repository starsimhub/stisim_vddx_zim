"""
Custom interventions for the discharge valuation
"""

import stisim as sti
import starsim as ss
import numpy as np
import pandas as pd


# %% Helper functions
def count(arr): return np.count_nonzero(arr)


class SyndromicMgmt(sti.SyndromicMgmt):
    def init_results(self):
        super().init_results()
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'ng_only', npts, dtype=int, scale=True, label='Only NG'),
            ss.Result(self.name, 'ct_only', npts, dtype=int, scale=True, label='Only CT'),
            ss.Result(self.name, 'tv_only', npts, dtype=int, scale=True, label='Only TV'),
            ss.Result(self.name, 'bv_only', npts, dtype=int, scale=True, label='Only BV'),
            ss.Result(self.name, 'ng_ct', npts, dtype=int, scale=True, label='NG & CT'),
            ss.Result(self.name, 'ng_tv', npts, dtype=int, scale=True, label='NG & TV'),
            ss.Result(self.name, 'ng_bv', npts, dtype=int, scale=True, label='NG & BV'),
            ss.Result(self.name, 'ct_tv', npts, dtype=int, scale=True, label='CT & TV'),
            ss.Result(self.name, 'ct_bv', npts, dtype=int, scale=True, label='CT & BV'),
            ss.Result(self.name, 'tv_bv', npts, dtype=int, scale=True, label='TV & BV'),
            ss.Result(self.name, 'ng_ct_tv', npts, dtype=int, scale=True, label='NG & CT & TV'),
            ss.Result(self.name, 'ng_ct_bv', npts, dtype=int, scale=True, label='NG & CT & BV'),
            ss.Result(self.name, 'ng_tv_bv', npts, dtype=int, scale=True, label='NG & TV & BV'),
            ss.Result(self.name, 'ct_tv_bv', npts, dtype=int, scale=True, label='CT & TV & BV'),
            ss.Result(self.name, 'ng_ct_tv_bv', npts, dtype=int, scale=True, label='NG & CT & TV & BV'),
            ss.Result(self.name, 'all_ng', npts, dtype=int, scale=True, label='All NG'),
            ss.Result(self.name, 'all_ct', npts, dtype=int, scale=True, label='All CT'),
            ss.Result(self.name, 'all_tv', npts, dtype=int, scale=True, label='All TV'),
            ss.Result(self.name, 'all_bv', npts, dtype=int, scale=True, label='All BV'),
        ]
        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti

        # Custom results
        self.results['ng_only'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_only'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['tv_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['bv_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_tv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_bv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_bv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['tv_bv'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_ct_bv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ng_tv_bv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv_bv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv_bv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['all_ng'][ti] = len((just_tested & ppl.ng.infected).uids)
        self.results['all_ct'][ti] = len((just_tested & ppl.ct.infected).uids)
        self.results['all_tv'][ti] = len((just_tested & ppl.tv.infected).uids)
        self.results['all_bv'][ti] = len((just_tested & ppl.bv.infected).uids)

        return


class Panel(sti.STITest):
    def __init__(self, treatments=None, diseases=None, pars=None, years=None, start=None, end=None, eligibility=None, product=None, name=None, label=None, **kwargs):
        super().__init__(years=years, start=start, end=end, eligibility=eligibility, product=product, name=name, label=label, **kwargs)
        self.default_pars(
            sens=dict(
                ng=ss.bernoulli(0.9),
                ct=ss.bernoulli(0.9),
                tv=ss.bernoulli(0.9),
            ),
            spec=dict(
                ng=ss.bernoulli(0.95),
                ct=ss.bernoulli(0.95),
                tv=ss.bernoulli(0.95),
            ),
        )
        self.update_pars(pars, **kwargs)

        # Store treatments and diseases
        self.treatments = treatments
        self.diseases = diseases
        self.disease_treatment_map = {t.disease:t for t in self.treatments}

        return

    def apply(self, sim, uids=None):
        """ Apply the testing intervention """
        # Apply if within the start years
        if (sim.year >= self.start) & (sim.year < self.end):

            for treatment in self.treatments:
                treatment.eligibility = ss.uids()  # Reset

            if uids is None:
                uids = self.get_testers(sim)
                self.ti_tested[uids] = sim.ti

            if len(uids):
                for disease in self.diseases:
                    inf = uids & disease.infected
                    sus = uids & disease.susceptible

                    sens = self.pars.sens[disease.name]
                    spec = self.pars.spec[disease.name]

                    true_pos, false_pos = sens.split(inf)
                    true_neg, false_neg = spec.split(sus)

                    tx = self.disease_treatment_map[disease.name]
                    tx.eligibility = true_pos | false_neg

            # Update results
            self.update_results()

        return


# %%  Algorithms

def make_testing(ng, ct, tv, bv, scenario='soc', end=2040):

    intv_year = 2000

    if scenario == 'soc':
        synd_end = intv_year
    else:
        synd_end = end

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
    tv_tx = sti.STITreatment(disease='tv', name='tv_tx', label='tv_tx')
    bv_tx = sti.STITreatment(disease='bv', name='bv_tx', label='bv_tx')
    treatments = [ng_tx, ct_tx, tv_tx, bv_tx]

    syndromic = SyndromicMgmt(
        end=synd_end,
        p_treat=0.8,
        diseases=[ng, ct, tv, bv],
        eligibility=seeking_care_discharge,
        treatments=treatments,
    )

    intvs = [syndromic, ng_tx, tv_tx, ct_tx, bv_tx]

    if scenario == 'panel':
        panel = Panel(
            start=intv_year,
            eligibility=seeking_care_discharge,
            diseases=[ng, ct, tv],
            treatments=treatments,
        )

        intvs += [panel]

    return intvs
