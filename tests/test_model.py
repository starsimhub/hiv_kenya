"""
Smoke test for the Kenya HIV model.
"""

import sys; sys.path.insert(0, '..')
import numpy as np
import sciris as sc
from hiv_model import make_sim


def test_hiv_model():
    sim = make_sim(n_agents=1000, use_calib=False, verbose=-1)
    sim.run()

    prev = sim.results.hiv['prevalence_15_49']
    art = sim.results.hiv.n_on_art
    assert np.all(prev >= 0), 'Prevalence should be non-negative'
    assert prev[-1] > prev[0], 'Prevalence should increase during the sim'
    assert art[0] == 0, 'No one on ART at simulation start'
    assert art[-1] > 0, 'People on ART at simulation end'
    return sim


if __name__ == '__main__':
    T = sc.timer()
    sim = test_hiv_model()
    T.toc()
