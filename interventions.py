"""
Create interventions
"""

# %% Imports and settings
import numpy as np
import starsim as ss
import pandas as pd
import stisim as sti
import sciris as sc


def get_testing_products():
    """
    Define HIV products and testing interventions
    """
    scaleup_years = np.arange(1990, 2021)  # Years for testing
    years = np.arange(1990, 2051)  # Years for simulation
    n_years = len(scaleup_years)
    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_years), np.linspace(0.75, 0.85, len(years) - n_years)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_years), np.linspace(0.85, 0.95, len(years) - n_years)])
    gp_prob = np.concatenate([np.linspace(0, 0.1, n_years), np.linspace(0.1, 0.1, len(years) - n_years)])

    # FSW agents who haven't been diagnosed or treated yet
    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    fsw_testing = sti.HIVTest(
        years=years,
        test_prob_data=fsw_prob,
        name='fsw_testing',
        eligibility=fsw_eligibility,
        label='fsw_testing',
    )

    # Non-FSW agents who haven't been diagnosed or treated yet
    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    other_testing = sti.HIVTest(
        years=years,
        test_prob_data=gp_prob,
        name='other_testing',
        eligibility=other_eligibility,
        label='other_testing',
    )

    # Agents whose CD4 count is below 200.
    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    low_cd4_testing = sti.HIVTest(
        years=years,
        test_prob_data=low_cd4_prob,
        name='low_cd4_testing',
        eligibility=low_cd4_eligibility,
        label='low_cd4_testing',
    )

    return fsw_testing, other_testing, low_cd4_testing


def make_hiv_intvs():

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    # n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    fsw_testing, other_testing, low_cd4_testing = get_testing_products()
    art = sti.ART(coverage_data=n_art, future_coverage={'year': 2024, 'prop': 0.97})
    # vmmc = sti.VMMC(coverage_data=n_vmmc)
    prep = sti.Prep(
        coverage=[0, 0.01, 0.5, 0.8],
        years=[2004, 2005, 2015, 2025],
        eff_prep=0.8,
    )
    interventions = [
        fsw_testing,
        other_testing,
        low_cd4_testing
    ]
    interventions += [art, prep]

    return interventions

