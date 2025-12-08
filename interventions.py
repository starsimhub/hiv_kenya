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

    partner_testing = sti.HIVTest(
        years=years,
        test_prob_data=0.9,
        name='partner_testing',
        eligibility=ss.uids(),  # Set by partner notification intervention
        label='partner_testing',
    )

    return fsw_testing, other_testing, low_cd4_testing, partner_testing


class PartnerNotification(ss.Intervention):
    """ Notify partners of people testing positive to HIV """

    def __init__(self, pars=None, eligibility=None, name=None, label=None, start=2026, **kwargs):
        super().__init__(eligibility=eligibility, name=name, label=label)
        self.define_pars(
            p_notify=dict(
                current=ss.bernoulli(p=0.5),  # Probability of notifying current partners
                previous=ss.bernoulli(p=0.05),  # Probability of notifying previous partners
            ),
            p_attends=dict(
                current=ss.bernoulli(p=0.5),  # Probability that current partners will attend
                previous=ss.bernoulli(p=0.01),  # Probability that previous partners will attend
            ),
        )
        self.update_pars(pars, **kwargs)

        # Store the current and prior network
        self.nws = None  # Initialized in init_pre
        self.start = start

        self.define_states(
            ss.FloatArr('ti_notified')
        )

        self.contacts = sc.objdict(
            current=sc.objdict(mf=None, fm=None),  # Current partners notified
            previous=sc.objdict(mf=None, fm=None),  # Current partners notified
        )

        return

    def init_pre(self, sim):
        super().init_pre(sim)
        self.nws = dict(
            current=sim.networks.structuredsexual,  # Current sexual network
            previous=sim.networks.priorpartners,  # Prior sexual network
        )

    def identify_contacts(self, uids):
        # Return UIDs of people that have been identified as contacts and should be notified

        # Find contacts
        for nwtype, nw in self.nws.items():
            m_edge_inds = np.isin(nw.p1, uids).nonzero()[-1]
            f_edge_inds = np.isin(nw.p2, uids).nonzero()[-1]
            m_idx = nw.p1[m_edge_inds]
            f_partners = nw.p2[m_edge_inds]
            f_idx = nw.p2[f_edge_inds]
            m_partners = nw.p1[f_edge_inds]

            # Females notified and attending
            notified_f = self.pars.p_notify[nwtype].filter(f_partners)
            attending_f = self.pars.p_attends[nwtype].filter(notified_f)
            fp_edge_inds = np.isin(nw.p2, attending_f).nonzero()[-1]
            successful_m_idx = nw.p1[fp_edge_inds] & m_idx
            mf_pairs = list(zip(successful_m_idx, attending_f))

            # Males notified and attending
            notified_m = self.pars.p_notify[nwtype].filter(m_partners)
            attending_m = self.pars.p_attends[nwtype].filter(notified_m)
            mp_edge_inds = np.isin(nw.p1, attending_m).nonzero()[-1]
            successful_f_idx = nw.p2[mp_edge_inds] & f_idx
            fm_pairs = list(zip(successful_f_idx, attending_m))

            # Store contacts
            self.ti_notified[attending_m | attending_f] = self.ti
            self.contacts[nwtype].mf = mf_pairs
            self.contacts[nwtype].fm = fm_pairs

        return

    def step(self):
        sim = self.sim

        if self.t.now('year') >= self.start:
            self.sim.interventions['partner_testing'].eligibility = ss.uids()  # Reset

            index_cases = self.eligibility(sim)

            if len(index_cases) > 0:
                self.identify_contacts(index_cases)

                # In this scenario, we test partners
                eligible_partners = []
                for pstatus in ['current', 'previous']:
                    for pairkey, clist in self.contacts[pstatus].items():
                        if len(clist) > 0:
                            eligible_partners += [partners[1] for partners in clist if not sim.diseases.hiv.diagnosed[partners[1]]]
                eligible_partners = ss.uids(set(eligible_partners))
                self.sim.interventions['partner_testing'].eligibility = eligible_partners

        return


def make_hiv_intvs(pn_pars=None):

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    # n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    fsw_testing, other_testing, low_cd4_testing, partner_testing = get_testing_products()
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

    if pn_pars is not None:
        # Optionally add partner notification
        # Partner treatment eligibility
        def just_diagnosed(sim):
            """ Return UIDs of people who have just been diagnosed - could add recency test here too  """
            new_diagnoses = (sim.diseases.hiv.ti_diagnosed == sim.diseases.hiv.ti).uids
            # previous_index = ... # Exclude people who were the original index case
            return new_diagnoses

        pn = PartnerNotification(
            **pn_pars,
            eligibility=just_diagnosed,
            name='notify_partners',
            label='notify_partners',
        )

        interventions += [
            pn,
            partner_testing,
        ]

    interventions += [art, prep]

    return interventions

