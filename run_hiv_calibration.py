"""
Run calibration for the HIV model
""" 
 
# Additions to handle numpy multithreading
import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# %% Imports and settings
import sciris as sc
import stisim as sti
import pandas as pd
from hiv_model import make_sim, make_sim_pars


# Run settings
debug = False  # If True, this will do smaller runs that can be run locally for debugging
n_trials = [1000, 2][debug]  # How many trials to run for calibration
n_workers = [50, 1][debug]    # How many cores to use
# storage = ["mysql://hpvsim_user@localhost/hpvsim_db", None][debug]  # Storage for calibrations
storage = None
do_shrink = True  # Whether to shrink the calibration results
make_stats = True  # Whether to make stats


def run_calibration(n_trials=None, n_workers=None, do_save=True):

    # Define the calibration parameters
    calib_pars = dict(
        hiv_beta_m2f=dict(low=0.008, high=0.02, guess=0.012),
        nw_prop_f0 = dict(low=0.55, high=0.9, guess=0.85),
        nw_prop_m0 = dict(low=0.50, high=0.9, guess=0.81),
        nw_f1_conc = dict(low=0.01, high=0.2, guess=0.01),
        nw_m1_conc = dict(low=0.01, high=0.2, guess=0.01),
        nw_p_pair_form = dict(low=0.4, high=0.9, guess=0.5),
    )

    # Make the sim
    sim = make_sim(verbose=-1, use_calib=False)
    data = pd.read_csv('data/kenya_hiv_calib.csv')
    extra_results = ['hiv.n_diagnosed', 'hiv.n_on_art', 'n_alive']

    # Make the calibration
    calib = sti.Calibration(
        calib_pars=calib_pars,
        build_fn=make_sim_pars,
        sim=sim,
        extra_results=extra_results,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=False, storage=storage, save_results=True,
    )

    calib.calibrate(load=True)

    return sim, calib


if __name__ == '__main__':

    sim, calib = run_calibration(n_trials=n_trials, n_workers=n_workers)
    print(f'Best pars are {calib.best_pars}')

    # Save the results
    print('Shrinking and saving...')
    if do_shrink:
        calib = calib.shrink(n_results=500)
        sc.saveobj(f'results/kenya_hiv_calib.obj', calib)
    else:
        sc.saveobj(f'results/kenya_hiv_calib.obj', calib)

    # Make stats
    if make_stats:
        print('Making stats...')
        from utils import percentiles
        df = calib.resdf
        df_stats = df.groupby(df.time).describe(percentiles=percentiles)
        sc.saveobj(f'results/kenya_hiv_calib_stats.df', df_stats)
        par_stats = calib.df.describe(percentiles=[0.05, 0.95])
        sc.saveobj(f'results/kenya_hiv_par_stats.df', par_stats)

    print('Done!')