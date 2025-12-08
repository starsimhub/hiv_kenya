"""
Create HIV model and interventions
"""

# %% Imports and settings
import sciris as sc
import pandas as pd
import stisim as sti
import starsim as ss
from interventions import make_hiv_intvs


def make_sim_pars(sim, calib_pars):
    if not sim.initialized: sim.init()
    hiv = sim.diseases.hiv
    nw = sim.networks.structuredsexual

    # Apply the calibration parameters
    for k, pars in calib_pars.items():  # Loop over the calibration parameters
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

        if 'hiv_' in k:  # HIV parameters
            k = k.replace('hiv_', '')  # Strip off indentifying part of parameter name
            hiv.pars[k] = v
        elif 'nw_' in k:  # Network parameters
            k = k.replace('nw_', '')  # As above
            if 'pair_form' in k:
                nw.pars[k].set(v)
            else:
                nw.pars[k] = v
        else:
            raise NotImplementedError(f'Parameter {k} not recognized')

    return sim


def make_sim(seed=1, stop=2030, verbose=1/12, analyzers=None, use_calib=True, analyze_network=False, par_idx=0):

    nw = sti.StructuredSexual(
        prop_f0=0.79,
        prop_m0=0.75,
        f1_conc=0.15,
        m1_conc=0.15,
        p_pair_form=0.5,
        condom_data=pd.read_csv(f'data/condom_use.csv'),
    )

    hiv = sti.HIV(
        beta_m2f=0.012,
        eff_condom=0.5,
        init_prev_data=pd.read_csv('data/init_prev_hiv.csv'),
        rel_init_prev=.5,
    )

    intvs = make_hiv_intvs()

    # Add network analyzers
    analyzers = sc.autolist(analyzers)
    analyzers += sti.sw_stats(diseases=['hiv'])
    if analyze_network:
        analyzers += sti.NetworkDegree(relationship_types=['partners', 'stable', 'casual'])
        analyzers += sti.RelationshipDurations()
        analyzers += sti.DebutAge()
        analyzers += sti.partner_age_diff()

    sim = sti.Sim(
        n_agents=10e3,
        start=1985,
        stop=stop,
        use_migration=True,
        datafolder='data/',
        demographics='kenya',
        diseases=hiv,
        rand_seed=seed,
        networks=[nw, ss.MaternalNet()],
        interventions=intvs,
        analyzers=analyzers,
        verbose=verbose,
    )

    # If using calibration parameters, update the simulation
    if use_calib:
        calib = sc.loadobj('results/kenya_hiv_calib.obj')
        calib_pars = calib.df.iloc[par_idx].to_dict()
        sim.init()
        sim = make_sim_pars(sim, calib_pars)
        print(f'Using calibration parameters for index {par_idx}')

    return sim


def run_msim(use_calib=True, n_pars=1, do_save=True):

    # Mave individual sims
    sims = sc.autolist()

    for par_idx in range(n_pars):
        sim = make_sim(use_calib=use_calib, par_idx=par_idx, verbose=-1)
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

    # Epi stats: save for all runs
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
    use_calib = False

    to_run = [
        'run_sim',
        # 'run_msim',
    ]

    if 'run_sim' in to_run:

        if do_run:
            sim = make_sim(use_calib=use_calib, analyze_network=False, verbose=1/12)
            sim.run()
            df = sim.to_df(resample='year', use_years=True, sep='.')  # Use dots to separate columns
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
