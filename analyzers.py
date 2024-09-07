"""
Analyzers for the discharging STI model
"""

# %% Imports and settings
import starsim as ss


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
            ss.Result(self.name, 'n_treated', npts, dtype=int, scale=False, label='Total treated'),
            ss.Result(self.name, 'ng', npts, dtype=int, scale=True, label='NG'),
            ss.Result(self.name, 'ct', npts, dtype=int, scale=True, label='CT'),
            ss.Result(self.name, 'tv', npts, dtype=int, scale=True, label='TV'),
            ss.Result(self.name, 'bv', npts, dtype=int, scale=True, label='BV'),
        ]

    def init_post(self):
        super().init_post()
        return

    def apply(self, sim):
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == sim.ti) | (ppl.ct.ti_seeks_care == sim.ti) | (ppl.tv.ti_seeks_care == sim.ti) | (ppl.bv.ti_seeks_care == sim.ti)
        for k in ['ng', 'ct', 'tv', 'bv']:
            self.results[k][sim.ti] = len((care & ~ppl[k].infected).uids)
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


