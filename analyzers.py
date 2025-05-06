"""
Analyzers for the discharging STI model
"""
import numpy as np

# %% Imports and settings
import starsim as ss
import sciris as sc


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
        self.define_results(
            ss.Result('new_ng_care_seekers', dtype=int, label='NG'),
            ss.Result('new_ct_care_seekers', dtype=int, label='CT'),
            ss.Result('new_tv_care_seekers', dtype=int, label='TV'),
            ss.Result('new_bv_care_seekers', dtype=int, label='BV'),
        )

    def step(self):
        sim = self.sim
        ti = self.ti
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == ti) | (ppl.ct.ti_seeks_care == ti) | (ppl.tv.ti_seeks_care == ti) | (ppl.bv.ti_seeks_care == ti)
        for k in ['ng', 'ct', 'tv', 'bv']:
            self.results[f'new_{k}_care_seekers'][ti] = len((care & ~ppl[k].infected).uids)
        return


class coinfection_stats(ss.Analyzer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'coinfection_stats'
        return

    def init_results(self):
        self.define_results(
            ss.Result('ng_only', dtype=int, label='Only NG'),
            ss.Result('ct_only', dtype=int, label='Only CT'),
            ss.Result('tv_only', dtype=int, label='Only TV'),
            ss.Result('bv_only', dtype=int, label='Only BV'),
            ss.Result('ng_ct', dtype=int, label='NG & CT'),
            ss.Result('ng_tv', dtype=int, label='NG & TV'),
            ss.Result('ng_bv', dtype=int, label='NG & BV'),
            ss.Result('ct_tv', dtype=int, label='CT & TV'),
            ss.Result('ct_bv', dtype=int, label='CT & BV'),
            ss.Result('tv_bv', dtype=int, label='TV & BV'),
            ss.Result('ng_ct_tv', dtype=int, label='NG & CT & TV'),
            ss.Result('ng_ct_bv', dtype=int, label='NG & CT & BD'),
            ss.Result('ng_tv_bv', dtype=int, label='NG & TV & BD'),
            ss.Result('ct_tv_bv', dtype=int, label='CT & TV & BD'),
            ss.Result('ng_ct_tv_bv', dtype=int, label='NG & CT & TV & BD'),
        )

    def step(self):
        sim = self.sim
        ti = self.ti
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == ti) | (ppl.ct.ti_seeks_care == ti) | (ppl.tv.ti_seeks_care == ti) | (ppl.bv.ti_seeks_care == ti)

        self.results['ng_only'][ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_only'][ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['tv_only'][ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['bv_only'][ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct'][ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_tv'][ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_bv'][ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv'][ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ct_bv'][ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['tv_bv'][ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv'][ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.bv.infected).uids)
        self.results['ng_ct_bv'][ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ng_tv_bv'][ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)
        self.results['ct_tv_bv'][ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

        self.results['ng_ct_tv_bv'][ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.bv.infected).uids)

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

    def init_results(self):
        self.define_results(
            ss.Result('new_symptoms', dtype=int, label='Symptomatic incidence'),
            ss.Result('n_symptomatic', dtype=int, label='Number with symptoms'),
            ss.Result('symp_prev', scale=False, label='Adult symptomatic prevalence'),
            ss.Result('symp_prev_f', scale=False, label='Vaginal discharge prevalence'),
            ss.Result('symp_prev_m', scale=False, label='Urethral discharge prevalence'),
            ss.Result('symp_prev_no_hiv', scale=False, label="Symptomatic prevalence HIV-"),
            ss.Result('symp_prev_has_hiv', scale=False, label="Symptomatic prevalence HIV+"),
            ss.Result('symp_prev_no_hiv_f', scale=False, label="Symptomatic prevalence HIV- F"),
            ss.Result('symp_prev_has_hiv_f', scale=False, label="Symptomatic prevalence HIV+ F"),
            ss.Result('symp_prev_no_hiv_m', scale=False, label="Symptomatic prevalence HIV- M"),
            ss.Result('symp_prev_has_hiv_m', scale=False, label="Symptomatic prevalence HIV+ M"),
        )
        return

    @staticmethod
    def cond_prob(numerator, denominator):
        numer = len((denominator & numerator).uids)
        denom = len(denominator.uids)
        out = sc.safedivide(numer, denom)
        return out

    def step(self):
        sim = self.sim
        ti = self.ti
        adults = (sim.people.age >= 15) & (sim.people.age <= 65)
        women = adults & sim.people.female
        men = adults & sim.people.male
        hiv = sim.diseases.hiv

        new_symptoms = (sim.people.ng.ti_symptomatic == self.ti) | (sim.people.ct.ti_symptomatic == self.ti) | (sim.people.tv.ti_symptomatic == self.ti) | (sim.people.bv.ti_symptomatic == self.ti)
        any_symptoms = sim.people.ng.symptomatic | sim.people.ct.symptomatic | sim.people.tv.symptomatic | sim.people.bv.symptomatic
        has_hiv = adults & hiv.infected  # Adults with HIV

        n_symp = any_symptoms & adults
        n_symp_f = any_symptoms & women
        n_symp_m = any_symptoms & men
        has_hiv_f = has_hiv & women  # Women with HIV
        has_hiv_m = has_hiv & men  # Men with HIV
        no_hiv = adults & hiv.susceptible  # Adults without HIV
        no_hiv_f = no_hiv & women  # Women without HIV
        no_hiv_m = no_hiv & men  # Men without HIV

        self.results['new_symptoms'][ti] = np.count_nonzero(new_symptoms)
        self.results['n_symptomatic'][ti] = np.count_nonzero(any_symptoms)
        self.results['symp_prev'][ti] = np.count_nonzero(n_symp) / np.count_nonzero(adults)
        self.results['symp_prev_f'][ti] = np.count_nonzero(n_symp_f) / np.count_nonzero(women)
        self.results['symp_prev_m'][ti] = np.count_nonzero(n_symp_m) / np.count_nonzero(men)

        self.results['symp_prev_no_hiv'][ti] = self.cond_prob(n_symp, no_hiv)
        self.results['symp_prev_has_hiv'][ti] = self.cond_prob(n_symp, has_hiv)
        self.results['symp_prev_no_hiv_f'][ti] = self.cond_prob(n_symp_f, no_hiv_f)
        self.results['symp_prev_has_hiv_f'][ti] = self.cond_prob(n_symp_f, has_hiv_f)
        self.results['symp_prev_no_hiv_m'][ti] = self.cond_prob(n_symp_m, no_hiv_m)
        self.results['symp_prev_has_hiv_m'][ti] = self.cond_prob(n_symp_m, has_hiv_m)

        for disease in ['ng', 'ct', 'tv']:
            if self.results['symp_prev_f'][ti] < sim.results[disease]['symp_prevalence_f'][ti]:
                errormsg = f'Overall symptomatic prevalence should not be lower than for disease {disease}'
                raise ValueError(errormsg)

