"""
Test that the HIV Kenya model runs successfully.
"""

import sys; sys.path.insert(0, '..') # Add parent folder to get hiv_model
import numpy as np
import sciris as sc
from hiv_model import make_sim


def test_hiv_model(do_plot=False):
    """Test that the model can be created and run without errors."""
    
    # Create and run the sim
    sim = make_sim(n_agents=1000)
    sim.run()
    
    # Do simple checks
    res = sim.results.hiv
    prev = res.prevalence
    art = res.n_on_art
    assert np.all(prev > 0), 'Expect nonzero prevalence at all timepoints'
    assert prev[-1] > prev[0], 'Expect prevalence to increase during the sim'
    assert art[0] == 0, 'Expect no one on ART at simulation start'
    assert art[-1] > 0, 'Expect people on ART at simulation end'
    
    if do_plot:
        sim.plot('hiv_prevalence_15_49')
    
    return sim


if __name__ == '__main__':
    T = sc.timer()
    do_plot = True
    
    sim = test_hiv_model(do_plot=do_plot)
    
    T.toc()
