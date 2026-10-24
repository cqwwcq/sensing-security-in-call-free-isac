"""Outage sweep behind Fig. 4 and Reviewer 1, Comment 5 (port of run_outage_sweep.m).

Example:
    python run_outage_sweep.py --set num_runs=2 num_err=500
"""
import os

from isac_sim.runner import limit_threads

limit_threads()

from isac_sim import outage  # noqa: E402
from isac_sim.geometry import setup_geometry  # noqa: E402
from isac_sim.runner import driver_args, make_cfg, run_all, save_npz  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    args = driver_args(__doc__, os.path.join(HERE, 'results', 'outage_results.npz'))
    cfg = make_cfg(outage.DRIVER_DEFAULTS, args.set)
    geo = setup_geometry(cfg)
    results, elapsed, ckpt = run_all(outage.one_run, cfg, geo, list(range(1, cfg.num_runs + 1)),
                                     args.out, args.workers)
    S = outage.collect(results, cfg)
    S['elapsed'] = elapsed
    save_npz(args.out, S, cfg)
    print(f'Saved {args.out} ({elapsed:.0f} s; checkpoints in {ckpt}).')
    outage.print_results(S)


if __name__ == '__main__':
    main()
