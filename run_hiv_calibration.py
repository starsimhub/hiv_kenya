"""
Optuna calibration entry point for the Kenya HIV model.
"""

import os
os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

import pandas as pd
import sciris as sc
import stisim as sti

from hiv_model import make_sim
from utils import percentiles


debug = False
n_trials = [1000, 2][debug]
n_workers = [50, 1][debug]
do_shrink = True
shrink_to = 500


def run_calibration(n_trials=n_trials, n_workers=n_workers):
    calib_pars = dict(
        hiv=dict(
            beta_m2f=dict(low=0.008, high=0.02, guess=0.012),
            eff_condom=dict(low=0.5, high=0.95, guess=0.75),
        ),
        structuredsexual=dict(
            prop_f0=dict(low=0.55, high=0.9, guess=0.85),
            prop_m0=dict(low=0.50, high=0.9, guess=0.81),
            f1_conc=dict(low=0.01, high=0.2, guess=0.01),
            m1_conc=dict(low=0.01, high=0.2, guess=0.01),
            p_pair_form=dict(low=0.4, high=0.9, guess=0.5),
        ),
    )

    sim = make_sim(verbose=-1, use_calib=False)
    data = pd.read_csv('data/kenya_hiv_calib.csv')
    extra_results = ['hiv.n_diagnosed', 'hiv.n_on_art', 'n_alive']

    calib = sti.Calibration(
        calib_pars=calib_pars,
        sim=sim,
        extra_results=extra_results,
        data=data,
        total_trials=n_trials, n_workers=n_workers,
        die=True, reseed=False, storage=None, save_results=True,
    )
    calib.calibrate(load=True)
    print(f'Best pars: {calib.best_pars}')

    return sim, calib


if __name__ == '__main__':
    sim, calib = run_calibration(n_trials=n_trials, n_workers=n_workers)

    print('Shrinking and saving...')
    if do_shrink:
        calib = calib.shrink(n_results=shrink_to)
    sc.saveobj('results/kenya_hiv_calib.obj', calib)

    print('Making stats...')
    df_stats = calib.resdf.groupby(calib.resdf.time).describe(percentiles=percentiles)
    sc.saveobj('results/kenya_hiv_calib_stats.df', df_stats)
    par_stats = calib.df.describe(percentiles=[0.05, 0.95])
    sc.saveobj('results/kenya_hiv_par_stats.df', par_stats)

    print('Done!')
