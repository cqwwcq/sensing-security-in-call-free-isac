"""Shared helpers of the sweep drivers: parallel runs, checkpoints, saving.

The MATLAB drivers run the channel realizations in a parfor loop and save
only when the sweep finishes.  Here every finished realization is written to
a checkpoint folder next to the output file, so an interrupted sweep resumes
where it stopped when started again with the same output file and settings.
"""
import argparse
import hashlib
import json
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from .config import apply_overrides, cfg_to_dict, parse_overrides, sim_config


def limit_threads():
    """One BLAS thread per process; the sweeps parallelize over realizations."""
    for var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ.setdefault(var, '1')


def driver_args(description, default_out, extra=None):
    """Command line common to all drivers: --out, --set key=value, --workers."""
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--out', default=default_out, help='output .npz file')
    ap.add_argument('--set', nargs='*', default=[], metavar='KEY=VALUE',
                    help='override settings of sim_config / the driver, '
                         'e.g. --set num_runs=1 "rho_list=[0.75, 1.05]"')
    ap.add_argument('--workers', type=int, default=None,
                    help='parallel processes (default: cfg.workers, 0 = all cores)')
    for args, kw in (extra or []):
        ap.add_argument(*args, **kw)
    return ap.parse_args()


def make_cfg(driver_defaults, overrides):
    """sim_config() plus the driver's own settings, then the overrides."""
    cfg = sim_config()
    for key, val in driver_defaults.items():
        setattr(cfg, key, val)
    return apply_overrides(cfg, parse_overrides(overrides))


def _json_default(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    return str(x)


def cfg_json(cfg):
    d = cfg_to_dict(cfg)
    return json.dumps(d, default=_json_default, sort_keys=True)


def run_all(one_run, cfg, geo, runs, out_file, workers=None, label='run'):
    """Run one_run(cfg, geo, r) for r in runs, in parallel, with checkpoints.

    Returns the list of per-run results in the order of runs.
    """
    d = json.loads(cfg_json(cfg))
    d.pop('workers', None)
    tag = hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()[:10]
    ckpt = f'{os.path.splitext(out_file)[0]}.parts_{tag}'
    os.makedirs(ckpt, exist_ok=True)
    path = {r: os.path.join(ckpt, f'{label}_{r:05d}.pkl') for r in runs}
    todo = [r for r in runs if not os.path.exists(path[r])]
    done = len(runs) - len(todo)
    if done:
        print(f'Resuming: {done} of {len(runs)} realizations already in {ckpt}', flush=True)

    nproc = cfg.workers if workers is None else workers
    if nproc == 0:
        nproc = os.cpu_count() or 1
    t0 = time.time()

    def store(r, res):
        tmp = path[r] + '.tmp'
        with open(tmp, 'wb') as f:
            pickle.dump(res, f)
        os.replace(tmp, path[r])

    def report(n):
        el = time.time() - t0
        rate = el / max(n, 1)
        left = rate * (len(todo) - n)
        print(f'  {label} {done + n}/{len(runs)} done, {el / 60:.1f} min elapsed, '
              f'about {left / 60:.1f} min left', flush=True)

    if nproc <= 1 or len(todo) <= 1:
        for n, r in enumerate(todo, 1):
            store(r, one_run(cfg, geo, r))
            report(n)
    else:
        with ProcessPoolExecutor(max_workers=nproc, initializer=limit_threads) as ex:
            futs = {ex.submit(one_run, cfg, geo, r): r for r in todo}
            for n, fut in enumerate(as_completed(futs), 1):
                store(futs[fut], fut.result())
                report(n)

    results = []
    for r in runs:
        with open(path[r], 'rb') as f:
            results.append(pickle.load(f))
    return results, time.time() - t0, ckpt


def save_npz(out_file, data, cfg):
    os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
    data = dict(data)
    data['cfg_json'] = np.array(cfg_json(cfg))
    np.savez_compressed(out_file, **data)


def load_npz(path):
    """Results as a dict; data['cfg'] is the saved configuration dict."""
    with np.load(path, allow_pickle=False) as f:
        data = {k: f[k] for k in f.files}
    data['cfg'] = json.loads(str(data.pop('cfg_json')))
    return data


def here():
    return os.path.dirname(os.path.abspath(sys.argv[0]))
