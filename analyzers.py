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
            ss.Result(self.name, 'vd', npts, dtype=int, scale=True, label='VD'),
        ]

    def init_post(self):
        super().init_post()
        return

    def apply(self, sim):
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == sim.ti) | (ppl.ct.ti_seeks_care == sim.ti) | (ppl.tv.ti_seeks_care == sim.ti) | (ppl.vd.ti_seeks_care == sim.ti)
        for k in ['ng', 'ct', 'tv', 'vd']:
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
        ]

    def init_post(self):
        super().init_post()
        return

    def apply(self, sim):
        ppl = sim.people
        care = (ppl.ng.ti_seeks_care == sim.ti) | (ppl.ct.ti_seeks_care == sim.ti) | (ppl.tv.ti_seeks_care == sim.ti) | (ppl.vd.ti_seeks_care == sim.ti)

        self.results['ng_only'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ct_only'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['tv_only'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['vd_only'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_tv'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_vd'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ct_tv'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ct_vd'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['tv_vd'][sim.ti] = len((care & ~ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct_tv'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ~ppl.vd.infected).uids)
        self.results['ng_ct_vd'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ~ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ng_tv_vd'][sim.ti] = len((care & ppl.ng.infected & ~ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)
        self.results['ct_tv_vd'][sim.ti] = len((care & ~ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        self.results['ng_ct_tv_vd'][sim.ti] = len((care & ppl.ng.infected & ppl.ct.infected & ppl.tv.infected & ppl.vd.infected).uids)

        not_infected = care & ~ppl.ng.infected & ~ppl.ct.infected & ~ppl.tv.infected & ~ppl.vd.infected
        if not_infected.any():
            errormsg = 'Should not be seeking care if not infected.'
            raise ValueError(errormsg)

        return


