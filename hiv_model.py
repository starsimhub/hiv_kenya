"""
Create HIV model and interventions for Kenya
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import pandas as pd
import starsim as ss
import stisim as sti


# Kenya-specific parameters
sim_pars = dict(
    start=1985,
    stop=2030,
    n_agents=10e3,
    use_migration=True,
)

sti_pars = dict(
    hiv=dict(
        beta_m2f=0.012,
        eff_condom=0.95,
        rel_init_prev=0.5,
    )
)

nw_pars = dict(
    prop_f0=0.8,
    prop_m0=0.75,
    f1_conc=0.15,
    m1_conc=0.15,
    p_pair_form=0.5,
)


def make_custom_interventions(test_years=None):
    """
    Create custom interventions for Kenya: HIV testing, ART, and PrEP.

    Note: ART is created here (rather than auto-loaded from art_coverage.csv)
    to allow setting future_coverage.

    Args:
        test_years (array): Years for testing coverage. Default: 1990-2050.

    Returns:
        list: List of intervention instances
    """
    if test_years is None:
        test_years = np.arange(1990, 2051)

    scaleup_end = min(2020, test_years[-1])
    scaleup_years = np.arange(test_years[0], scaleup_end + 1)
    n_scaleup = len(scaleup_years)
    n_future = len(test_years) - n_scaleup

    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_scaleup), np.linspace(0.75, 0.85, n_future)])
    gp_prob = np.concatenate([np.linspace(0, 0.1, n_scaleup), np.linspace(0.1, 0.1, n_future)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_scaleup), np.linspace(0.85, 0.95, n_future)])

    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    # Testing
    testing = [
        sti.HIVTest(years=test_years, test_prob_data=fsw_prob, name='fsw_testing', eligibility=fsw_eligibility, label='fsw_testing'),
        sti.HIVTest(years=test_years, test_prob_data=gp_prob, name='other_testing', eligibility=other_eligibility, label='other_testing'),
        sti.HIVTest(years=test_years, test_prob_data=low_cd4_prob, name='low_cd4_testing', eligibility=low_cd4_eligibility, label='low_cd4_testing'),
    ]

    # ART
    data_path = sc.thispath() / 'data'
    n_art = pd.read_csv(data_path / 'n_art.csv').set_index('year')
    art = sti.ART(coverage_data=n_art, future_coverage={'year': 2024, 'prop': 0.97})

    # PrEP
    prep = sti.Prep(
        coverage=[0, 0.01, 0.5, 0.8],
        years=[2004, 2005, 2015, 2025],
        eff_prep=0.8,
    )

    return testing + [art, prep]


def make_sim_pars(sim, calib_pars):
    """Apply calibration parameters to a simulation."""
    if not sim.initialized: sim.init()
    hiv = sim.diseases.hiv
    nw = sim.networks.structuredsexual

    for k, pars in calib_pars.items():
        if k == 'rand_seed':
            sim.pars.rand_seed = v
            continue

        elif k in ['index', 'mismatch']:
            continue

        if isinstance(pars, dict):
            v = pars['value']
        elif sc.isnumber(pars):
            v = pars
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

        if 'hiv_' in k:
            k = k.replace('hiv_', '')
            hiv.pars[k] = v
        elif 'nw_' in k:
            k = k.replace('nw_', '')
            if 'pair_form' in k:
                nw.pars[k].set(v)
            else:
                nw.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def make_sim(**kwargs):
    """
    Create a Kenya HIV simulation.

    Uses data_path to auto-load init_prev and condom_use data via DataLoader.
    Custom interventions (testing, ART with future_coverage, PrEP) are created
    separately and merged with any user-provided interventions.

    Args:
        **kwargs: Override any Sim parameters (e.g. verbose, rand_seed, stop, analyzers, interventions)

    Returns:
        sti.Sim: Configured simulation instance
    """
    intvs = make_custom_interventions()
    user_intvs = sc.tolist(kwargs.pop('interventions', []))

    # Default analyzers
    analyzers = sc.tolist(kwargs.pop('analyzers', []))
    analyzers += sti.sw_stats(diseases=['hiv'])

    sim = sti.Sim(
        location='kenya',
        diseases='hiv',
        data_path=sc.thispath() / 'data',
        sim_pars=sim_pars,
        nw_pars=nw_pars,
        sti_pars=sti_pars,
        interventions=intvs + user_intvs,
        analyzers=analyzers,
        **kwargs,
    )
    return sim


def run_msim(use_calib=True, n_pars=1, do_save=True):
    """Run multiple simulations, optionally applying calibration parameters."""
    calib = sc.loadobj('results/kenya_hiv_calib.obj') if use_calib else None

    sims = sc.autolist()
    for par_idx in range(n_pars):
        sim = make_sim(verbose=-1)
        if use_calib:
            calib_pars = calib.df.iloc[par_idx].to_dict()
            sim.init()
            sim = make_sim_pars(sim, calib_pars)
            print(f'Using calibration parameters for index {par_idx}')
        sim.par_idx = par_idx
        sims += sim
    sims = ss.parallel(sims).sims

    if do_save:
        dfs = sc.autolist()
        for sim in sims:
            par_idx = sim.par_idx
            df = sim.to_df(resample='year', use_years=True, sep='.')
            df['res_no'] = par_idx
            dfs += df
        df = pd.concat(dfs)
        sc.saveobj(f'results/msim.df', df)

    return sims


def save_stats(sims, resfolder='results'):
    """Save age/sex stratified epi stats and SW stats."""

    dfs = sc.autolist()
    for sim in sims:

        par_idx = sim.par_idx

        # Save age/sex epi results
        age_bins = sim.diseases.hiv.age_bins
        sex_labels = {'f': 'Female', 'm': 'Male'}
        disease = 'hiv'
        for sex in ['f', 'm']:
            dd = dict()
            for ab1, ab2 in zip(age_bins[:-1], age_bins[1:]):
                age = str(ab1) + '-' + str(ab2)
                if ab1 == 65:
                    age = '65+'  # Combine the last two age groups
                dd['age'] = [age]
                dd['sex'] = sex_labels[sex]
                dd['prevalence'] = sim.results[disease][f'prevalence_{sex}_{ab1}_{ab2}'][-1]
                dd['new_infections'] = sim.results[disease][f'new_infections_{sex}_{ab1}_{ab2}'][-120:].mean()
                dd['par_idx'] = par_idx
                dfs += pd.DataFrame(dd)
    epi_df = pd.concat(dfs)
    sc.saveobj(f'{resfolder}/epi_df.df', epi_df)

    # Save SW stats
    sim = [sim for sim in sims if sim.par_idx == 0][0]
    sw_res = sim.results['sw_stats']
    sw_df = sw_res.to_df(resample='year', use_years=True, sep='.')
    sc.saveobj(f'{resfolder}/sw_df.df', sw_df)

    return


# %% Run as a script
if __name__ == '__main__':

    # SETTINGS
    debug = False
    seed = 1
    do_save = True
    do_run = True
    do_plot = True
    use_calib = True

    to_run = [
        # 'run_sim',
        'run_msim',
    ]

    if 'run_sim' in to_run:

        if do_run:
            sim = make_sim(rand_seed=seed, verbose=1/12)
            if use_calib:
                calib = sc.loadobj('results/kenya_hiv_calib.obj')
                calib_pars = calib.df.iloc[0].to_dict()
                sim.init()
                sim = make_sim_pars(sim, calib_pars)
                print('Using calibration parameters')
            sim.run()
            df = sim.to_df(resample='year', use_years=True, sep='.')
            df.index = df['timevec']
            if do_save:
                sc.saveobj(f'results/kenya_sim.df', df)
                sc.saveobj(f'results/kenya.sim', sim)
        else:
            df = sc.loadobj(f'results/kenya_sim.df')

        if do_plot:
            from plot_sims import plot_hiv_sims
            plot_hiv_sims(df, start_year=1985, title='hiv_plots')

    if 'run_msim' in to_run:
        n_pars = 50 if not debug else 2
        if do_run:
            sims = run_msim(use_calib=use_calib, n_pars=n_pars, do_save=do_save)
        else:
            sims = None

        if do_save and sims is not None:
            save_stats(sims, resfolder='results')
