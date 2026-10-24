"""Synchronization sweep, Reviewer 2, Comment 4 (port of run_sync_sweep.m).

Example:
    python run_sync_sweep.py --set num_runs=2
"""
import os

from isac_sim.runner import limit_threads

limit_threads()

from isac_sim import sync  # noqa: E402
from isac_sim.geometry import setup_geometry  # noqa: E402
from isac_sim.runner import driver_args, make_cfg, run_all, save_npz  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    args = driver_args(__doc__, os.path.join(HERE, 'results', 'sync_results.npz'))
    cfg = make_cfg(sync.DRIVER_DEFAULTS, args.set)
    geo = setup_geometry(cfg)
    results, elapsed, ckpt = run_all(sync.one_run, cfg, geo, list(range(1, cfg.num_runs + 1)),
                                     args.out, args.workers)
    S = sync.collect(results, cfg)
    S['elapsed'] = elapsed
    save_npz(args.out, S, cfg)
    print(f"Saved {args.out} ({S['n_kept']} of {cfg.num_runs} runs kept, {elapsed:.0f} s).")
    sync.print_results(S)


if __name__ == '__main__':
    main()
