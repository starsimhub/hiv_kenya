"""
Smoke test for the Kenya HIV model.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import sciris as sc
from hiv_model import make_sim


def test_hiv_model():
    # 5000 agents keeps the smoke test stable given fp.Sim scales up to
    # Kenya-national population (~13k per agent); smaller n_agents leaves
    # too few raw-agent events to survive stochastic rounding on ART.
    sim = make_sim(n_agents=5000, use_calib=False, verbose=-1)
    sim.run()

    prev = sim.results.hiv['prevalence_15_49']
    art = sim.results.hiv.n_on_art
    assert np.all(prev >= 0), 'Prevalence should be non-negative'
    assert np.max(prev) > prev[0], 'Prevalence should peak above starting level during the sim'
    assert art[0] == 0, 'No one on ART at simulation start'
    assert np.max(art) > 0, 'Someone should be on ART at some point in the sim'
    return sim


if __name__ == '__main__':
    T = sc.timer()
    sim = test_hiv_model()
    T.toc()
