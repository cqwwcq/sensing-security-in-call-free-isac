"""Communication tradeoff, Reviewer 3, Comment 3 (port of run_comm_tradeoff.m).

Examples:
    python run_comm_tradeoff.py
    python run_comm_tradeoff.py --out results/comm_tradeoff_rho030.npz --set rho_ref=0.3
"""
import os

from isac_sim.runner import limit_threads

limit_threads()

from isac_sim import tradeoff  # noqa: E402
from isac_sim.geometry import setup_geometry  # noqa: E402
from isac_sim.runner import driver_args, make_cfg, run_all, save_npz  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    args = driver_args(__doc__, os.path.join(HERE, 'results', 'comm_tradeoff.npz'))
    cfg = make_cfg(tradeoff.DRIVER_DEFAULTS, args.set)
    geo = setup_geometry(cfg)
    results, elapsed, ckpt = run_all(tradeoff.one_run, cfg, geo, list(range(1, cfg.num_runs + 1)),
                                     args.out, args.workers)
    S = tradeoff.collect(results, cfg)
    S['elapsed'] = elapsed
    save_npz(args.out, S, cfg)
    print(f"Saved {args.out} ({S['n_kept']} of {cfg.num_runs} runs kept, {elapsed:.0f} s).")
    tradeoff.print_results(S, cfg)


if __name__ == '__main__':
    main()
