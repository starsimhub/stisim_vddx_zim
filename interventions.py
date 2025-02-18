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


class SyndromicMgmt(sti.STITest):
    def __init__(self, pars=None, treatments=None, diseases=None, outcome_treatment_map=None, treat_prob_data=None, years=None, start=None, stop=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(years=years, start=start, stop=stop, eligibility=eligibility, name=name, label=label)
        self.define_pars(
            tx_mix=dict(
                all3=[0.80, 0.05],
                ngct=[0.05, 0.80],
                mtnz=[0.05, 0.00],
                none=[0.10, 0.15],
            ),
            tx_dist_f=ss.choice(a=4),
            tx_dist_m=ss.choice(a=4),
            dt_scale=False,
        )
        self.update_pars(pars, **kwargs)
        self.fvals = [v[0] for v in self.pars.tx_mix.values()]
        self.mvals = [v[1] for v in self.pars.tx_mix.values()]

        # Store treatments and diseases
        self.treatments = sc.tolist(treatments)
        self.diseases = diseases
        if outcome_treatment_map is None:
            outcome_treatment_map = dict(
                all3=self.treatments,
                ngct=[self.treatments[0], self.treatments[1]],
                mtnz=[self.treatments[1]],
                none=[],
            )
        self.outcome_treatment_map = outcome_treatment_map

        self.define_states(
            ss.FloatArr('ti_referred'),
            ss.FloatArr('ti_dismissed'),
        )
        self.treat_prob_data = treat_prob_data
        self.treat_prob = None
        self.treated_by_uid = None

        # Interim results
        self.sti_results = sc.objdict(
            new_ng_only=0,
            new_ct_only=0,
            new_tv_only=0,
            new_bv_only=0,
            new_ng_ct=0,
            new_ng_tv=0,
            new_ng_bv=0,
            new_ct_tv=0,
            new_ct_bv=0,
            new_tv_bv=0,
            new_ng_ct_tv=0,
            new_ng_ct_bv=0,
            new_ng_tv_bv=0,
            new_ct_tv_bv=0,
            new_ng_ct_tv_bv=0,
            new_all_ng=0,
            new_all_ct=0,
            new_all_tv=0,
            new_all_bv=0,
        )
        for k in self.sti_results.keys():
            self.sti_results[k+'_f'] = 0
            self.sti_results[k+'_m'] = 0

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.pars.tx_dist_f.set(p=self.fvals)
        self.pars.tx_dist_m.set(p=self.mvals)
        return

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        sexkeys = ['', 'f', 'm']
        for sk in sexkeys:
            skk = '' if sk == '' else f'_{sk}'
            skl = '' if sk == '' else f' - {sk.upper()}'
            results += [
                ss.Result('new_care_seekers'+skk, dtype=int, label="Care seekers"+skl),
                ss.Result('new_tx0'+skk, dtype=int, label="No treatment"+skl),
                ss.Result('new_tx1'+skk, dtype=int, label="1 treatment"+skl),
                ss.Result('new_tx2'+skk, dtype=int, label="2 treatment"+skl),
                ss.Result('new_tx3'+skk, dtype=int, label="3 treatments"+skl),
                ss.Result('new_sti1'+skk, dtype=int, label='1 STI'+skl),
                ss.Result('new_sti2'+skk, dtype=int, label='2 STIs'+skl),
                ss.Result('new_sti3'+skk, dtype=int, label='3 STIs'+skl),
                ss.Result('new_sti4'+skk, dtype=int, label='4 STIs'+skl),
            ]
        self.define_results(*results)
        return

    def step(self, uids=None):
        sim = self.sim
        self.treated_by_uid = None

        # If this intervention has stopped, reset eligibility for all associated treatments
        if sim.now >= self.stop:
            for treatment in self.treatments:
                treatment.eligibility = ss.uids()  # Reset
            return

        if sim.now >= self.start:

            if uids is None:
                uids = self.check_eligibility()
                self.ti_tested[uids] = self.ti

            if len(uids):
                f_uids = uids[sim.people.female[uids]]
                m_uids = uids[sim.people.male[uids]]

                # Determine treatment outcome for each agent
                outcomes_f = self.pars.tx_dist_f.rvs(f_uids)
                outcomes_m = self.pars.tx_dist_m.rvs(m_uids)

                # Treatment outcomes
                outcomes = dict(
                    all3=f_uids[outcomes_f == 0] | m_uids[outcomes_m == 0],
                    ngct=f_uids[outcomes_f == 1] | m_uids[outcomes_m == 1],
                    mtnz=f_uids[outcomes_f == 2] | m_uids[outcomes_m == 2],
                    none=f_uids[outcomes_f == 3] | m_uids[outcomes_m == 3],
                )

                # Update treatment eligibility
                for outcome, txs in self.outcome_treatment_map.items():
                    for tx in txs:
                        tx.eligibility = tx.eligibility | outcomes[outcome]

                # Update states: time referred to treatment for anyone referred
                referred_uids = outcomes['all3'] | outcomes['ngct'] | outcomes['mtnz']
                dismissed_uids = outcomes['none']
                self.ti_referred[referred_uids] = self.ti
                self.ti_dismissed[dismissed_uids] = self.ti
                self.treated_by_uid = outcomes

            self.store_results()

            return

    def store_results(self):
        """
        This has a different name to the usual update_results method because we want to ensure
        that it is called BEFORE the treatments are applied, so that we record who was infected.
        """
        ti = self.ti
        ppl = self.sim.people
        just_tested = self.ti_tested == ti
        self.results['new_care_seekers'][ti] += count(just_tested)
        self.results['new_care_seekers_f'][ti] += count(just_tested & ppl.female)
        self.results['new_care_seekers_m'][ti] += count(just_tested & ppl.male)

        # Record the number of people who received 0-3 treatments
        sexdict = {'': 'alive', 'f': 'female', 'm': 'male'}
        if self.treated_by_uid is not None:
            for sk, sl in sexdict.items():
                skk = '' if sk == '' else f'_{sk}'
                self.results['new_tx0'+skk][ti] += count(ppl[sl][self.treated_by_uid['none']])
                self.results['new_tx1'+skk][ti] += count(ppl[sl][self.treated_by_uid['mtnz']])
                self.results['new_tx2'+skk][ti] += count(ppl[sl][self.treated_by_uid['ngct']])
                self.results['new_tx3'+skk][ti] += count(ppl[sl][self.treated_by_uid['all3']])

        # Record
        for sk, sl in sexdict.items():
            skk = '' if sk == '' else f'_{sk}'

            self.sti_results['new_ng_only'+skk] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ct_only'+skk] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_tv_only'+skk] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_bv_only'+skk] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)

            self.sti_results['new_ng_ct'+skk] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ng_tv'+skk] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ng_bv'+skk] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ct_tv'+skk] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ct_bv'+skk] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_tv_bv'+skk] = len((just_tested & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)

            self.sti_results['new_ng_ct_tv'+skk] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ng_ct_bv'+skk] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ng_tv_bv'+skk] = len((just_tested & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)
            self.sti_results['new_ct_tv_bv'+skk] = len((just_tested & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)

            self.sti_results['new_ng_ct_tv_bv'+skk] = len((just_tested & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected & ppl[sl]).uids)

            self.sti_results['new_all_ng'+skk] = len((just_tested & ppl.ng.infected & ppl[sl]).uids)
            self.sti_results['new_all_ct'+skk] = len((just_tested & ppl.ct.infected & ppl[sl]).uids)
            self.sti_results['new_all_tv'+skk] = len((just_tested & ppl.tv.infected & ppl[sl]).uids)
            self.sti_results['new_all_bv'+skk] = len((just_tested & ppl.bv.infected & ppl[sl]).uids)

            self.results['new_sti1'+skk][ti] = (self.sti_results['new_ng_only'+skk] +
                                                self.sti_results['new_ct_only'+skk] +
                                                self.sti_results['new_tv_only'+skk] +
                                                self.sti_results['new_bv_only'+skk])
            self.results['new_sti2'+skk][ti] = (self.sti_results['new_ng_ct'+skk] +
                                            self.sti_results['new_ng_tv'+skk] +
                                            self.sti_results['new_ng_bv'+skk] +
                                            self.sti_results['new_ct_tv'+skk] +
                                            self.sti_results['new_ct_bv'+skk] +
                                            self.sti_results['new_tv_bv'+skk])
            self.results['new_sti3'+skk][ti] = (self.sti_results['new_ng_ct_tv'+skk] +
                                            self.sti_results['new_ng_ct_bv'+skk] +
                                            self.sti_results['new_ng_tv_bv'+skk] +
                                            self.sti_results['new_ct_tv_bv'+skk])

            self.results['new_sti4'+skk][ti] = self.sti_results['new_ng_ct_tv_bv'+skk]

        return


# %%  Algorithms
def make_tx_mix(scenario):
    if scenario == 'treat100':
        tx_mix = dict(
            all3=[1.00, 0.00],
            ngct=[0.00, 1.00],
            mtnz=[0.00, 0.00],
            none=[0.00, 0.00],
        )
    elif scenario == 'treat80':
        tx_mix = dict(
            all3=[0.80, 0.00],
            ngct=[0.05, 0.80],
            mtnz=[0.05, 0.00],
            none=[0.10, 0.20],
        )
    elif scenario == 'treat50':
        tx_mix = dict(
            all3=[0.500, 0.00],
            ngct=[0.125, 0.50],
            mtnz=[0.125, 0.00],
            none=[0.250, 0.50],
        )
    return tx_mix


def make_testing(ng, ct, tv, bv, scenario=None, poc=None, stop=2040):

    intv_year = 2027

    # Handle inputs
    synd_end = intv_year if poc else stop

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
    outcome_treatment_map = dict(
        all3=treatments,
        ngct=[ng_tx, ct_tx],
        mtnz=[metronidazole],
        none=[],
    )
    tx_mix = make_tx_mix(scenario)

    # Syndromic management
    syndromic = SyndromicMgmt(
        tx_mix=tx_mix,
        stop=synd_end,
        diseases=[ng, ct, tv, bv],
        eligibility=seeking_care_discharge,
        treatments=treatments,
        outcome_treatment_map=outcome_treatment_map,
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
            name='panel',
            label='panel',
            start=intv_year,
            diseases=[ng, ct, tv, bv],
            eligibility=seeking_care_discharge,
            treatments=treatments,
            disease_treatment_map=disease_treatment_map,
        )
        intvs = [syndromic, panel, ng_tx, ct_tx, metronidazole]

    return intvs
