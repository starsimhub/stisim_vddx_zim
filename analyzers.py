"""
Analyzers for the discharging STI model
"""
import numpy as np

# %% Imports and settings
import starsim as ss
import stisim as sti
from utils import count


class overtreatment_stats(ss.Analyzer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'overtreatment_stats'
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.init_results()
        return

    def init_results(self):
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_ng_care_seekers', npts, dtype=int, scale=True, label='NG'),
            ss.Result(self.name, 'new_ct_care_seekers', npts, dtype=int, scale=True, label='CT'),
            ss.Result(self.name, 'new_tv_care_seekers', npts, dtype=int, scale=True, label='TV'),
            ss.Result(self.name, 'new_bv_care_seekers', npts, dtype=int, scale=True, label='BV'),
        ]

    def init_post(self):
        super().init_post()
        return

    def apply(self, sim):
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == sim.ti) | (ppl.ct.ti_seeks_care == sim.ti) | (ppl.tv.ti_seeks_care == sim.ti) | (ppl.bv.ti_seeks_care == sim.ti)
        for k in ['ng', 'ct', 'tv', 'bv']:
            self.results[f'new_{k}_care_seekers'][sim.ti] = len((care & ~ppl[k].infected).uids)
        return


class coinfection_stats(ss.Analyzer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'coinfection_stats'
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.init_results()
        return

    def init_results(self):
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
            ss.Result(self.name, 'ng_ct_bv', npts, dtype=int, scale=True, label='NG & CT & BD'),
            ss.Result(self.name, 'ng_tv_bv', npts, dtype=int, scale=True, label='NG & TV & BD'),
            ss.Result(self.name, 'ct_tv_bv', npts, dtype=int, scale=True, label='CT & TV & BD'),
            ss.Result(self.name, 'ng_ct_tv_bv', npts, dtype=int, scale=True, label='NG & CT & TV & BD'),
        ]

    def init_post(self):
        super().init_post()
        return

    def apply(self, sim):
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == sim.ti) | (ppl.ct.ti_seeks_care == sim.ti) | (ppl.tv.ti_seeks_care == sim.ti) | (ppl.bv.ti_seeks_care == sim.ti)

        self.results['ng_only'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_only'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['tv_only'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['bv_only'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_tv'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_bv'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_bv'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['tv_bv'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_ct_bv'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ng_tv_bv'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv_bv'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv_bv'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        not_infected = care & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected
        if not_infected.any():
            errormsg = 'Should not be seeking care if not infected.'
            raise ValueError(errormsg)

        return


class total_symptomatic(ss.Analyzer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'total_symptomatic'
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.init_results()
        return

    def init_results(self):
        npts = self.sim.npts
        self.results += [
            ss.Result(self.name, 'new_symptoms', npts, dtype=int, scale=True, label='Symptomatic incidence'),
            ss.Result(self.name, 'n_symptomatic', npts, dtype=int, scale=True, label='Number with symptoms'),
            ss.Result(self.name, 'symp_prev', npts, dtype=float, scale=False, label='Adult symptomatic prevalence'),
            ss.Result(self.name, 'symp_prev_f', npts, dtype=float, scale=False, label='Vaginal discharge prevalence'),
            ss.Result(self.name, 'symp_prev_m', npts, dtype=float, scale=False, label='Urethral discharge prevalence'),
        ]
        return

    def apply(self, sim):
        adults = (sim.people.age >= 15) & (sim.people.age <= 65)
        women = adults & sim.people.female
        men = adults & sim.people.male

        new_symptoms = (sim.people.ng.ti_symptomatic == sim.ti) | (sim.people.ct.ti_symptomatic == sim.ti) | (sim.people.tv.ti_symptomatic == sim.ti) | (sim.people.bv.ti_symptomatic == sim.ti)
        any_symptoms = sim.people.ng.symptomatic | sim.people.ct.symptomatic | sim.people.tv.symptomatic | sim.people.bv.symptomatic

        n_symp = any_symptoms & adults
        n_symp_f = any_symptoms & women
        n_symp_m = any_symptoms & men

        self.results['new_symptoms'][sim.ti] = np.count_nonzero(new_symptoms)
        self.results['n_symptomatic'][sim.ti] = np.count_nonzero(any_symptoms)
        self.results['symp_prev'][sim.ti] = np.count_nonzero(n_symp) / np.count_nonzero(adults)
        self.results['symp_prev_f'][sim.ti] = np.count_nonzero(n_symp_f) / np.count_nonzero(women)
        self.results['symp_prev_m'][sim.ti] = np.count_nonzero(n_symp_m) / np.count_nonzero(men)

        for disease in ['ng', 'ct', 'tv', 'bv']:
            if self.results['symp_prev_f'][sim.ti] < sim.results[disease]['female_symp_adult_prevalence'][sim.ti]:
                errormsg = f'Overall symptomatic prevalence should not be lower than for disease {disease}'
                raise ValueError(errormsg)

