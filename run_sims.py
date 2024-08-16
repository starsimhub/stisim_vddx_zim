# %% Imports and settings
import starsim as ss
from model import make_sim


def run_sims(seed=None, n_runs=None):
    sims = [make_sim(seed=seed + i, verbose=0.01) for i in range(n_runs)]
    sims = ss.parallel(sims).sims
    return sims


if __name__ == '__main__':

    # SETTINGS
    debug = False
    n_runs = [100, 2][debug]  # Number of runs when using multisim
    seed = 1
    sims = run_sims(seed=seed, n_runs=n_runs)

