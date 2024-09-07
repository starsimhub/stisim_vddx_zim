"""
Custom interventions for the discharge valuation
"""

import stisim as sti
import starsim as ss
import numpy as np


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
            ss.Result(self.name, 'vd_only', npts, dtype=int, scale=True, label='Only VD'),
            ss.Result(self.name, 'ng_ct', npts, dtype=int, scale=True, label='NG & CT'),
            ss.Result(self.name, 'ng_tv', npts, dtype=int, scale=True, label='NG & TV'),
            ss.Result(self.name, 'ng_vd', npts, dtype=int, scale=True, label='NG & VD'),
            ss.Result(self.name, 'ct_tv', npts, dtype=int, scale=True, label='CT & TV'),
            ss.Result(self.name, 'ct_vd', npts, dtype=int, scale=True, label='CT & VD'),
            ss.Result(self.name, 'tv_vd', npts, dtype=int, scale=True, label='TV & VD'),
            ss.Result(self.name, 'ng_ct_tv', npts, dtype=int, scale=True, label='NG & CT & TV'),
            ss.Result(self.name, 'ng_ct_vd', npts, dtype=int, scale=True, label='NG & CT & VD'),
            ss.Result(self.name, 'ng_tv_vd', npts, dtype=int, scale=True, label='NG & TV & VD'),
            ss.Result(self.name, 'ct_tv_vd', npts, dtype=int, scale=True, label='CT & TV & VD'),
            ss.Result(self.name, 'ng_ct_tv_vd', npts, dtype=int, scale=True, label='NG & CT & TV & VD'),
            ss.Result(self.name, 'all_ng', npts, dtype=int, scale=True, label='All NG'),
            ss.Result(self.name, 'all_ct', npts, dtype=int, scale=True, label='All CT'),
            ss.Result(self.name, 'all_tv', npts, dtype=int, scale=True, label='All TV'),
            ss.Result(self.name, 'all_vd', npts, dtype=int, scale=True, label='All VD'),
        ]
        return

    def update_results(self):
        super().update_results()
        ti = self.sim.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti

        # Custom results
        self.results['ng_only'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ct_only'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['tv_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['vd_only'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_tv'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_vd'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ct_tv'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ct_vd'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['tv_vd'][ti] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct_tv'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_ct_vd'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ng_tv_vd'][ti] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ct_tv_vd'][ti] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct_tv_vd'][ti] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        self.results['all_ng'][ti] = len((just_tested & ppl.ng.infected).uids)
        self.results['all_ct'][ti] = len((just_tested & ppl.ct.infected).uids)
        self.results['all_tv'][ti] = len((just_tested & ppl.tv.infected).uids)
        self.results['all_vd'][ti] = len((just_tested & ppl.vd.infected).uids)

        return

# %%  Algorithms

def make_testing(diseases):

    # Testing interventions
    def seeking_care_discharge(sim):
        ng_care = sim.diseases.ng.symptomatic & (sim.diseases.ng.ti_seeks_care == sim.ti)
        tv_care = sim.diseases.tv.symptomatic & (sim.diseases.tv.ti_seeks_care == sim.ti)
        ct_care = sim.diseases.ct.symptomatic & (sim.diseases.ct.ti_seeks_care == sim.ti)
        vd_care = sim.diseases.vd.symptomatic & (sim.diseases.vd.ti_seeks_care == sim.ti)
        return (ng_care | tv_care | ct_care | vd_care).uids

    ng_tx = sti.GonorrheaTreatment(
        rel_treat_unsucc=0.05,
        rel_treat_unneed=0.01,
    )
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
