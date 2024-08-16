"""
Analyzers for the discharging STI model
"""

# %% Imports and settings
import sciris as sc
import starsim as ss


class coinfection_stats(ss.Analyzer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'coinfection_stats'
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.init_results()
        return

    def init_post(self):
        super().init_post()
        return

    def finalize(self, sim):
        super().finalize()
        ppl = sim.people

        return


