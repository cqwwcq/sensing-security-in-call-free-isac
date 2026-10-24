"""Detection sweep behind Figs. 2-3 (port of run_detection_sweep.m).

Examples:
    python run_detection_sweep.py
    python run_detection_sweep.py --out results/smoke.npz --set num_runs=1 "rho_list=[0.75, 1.05]"
"""
import os

from isac_sim.runner import limit_threads

limit_threads()

from isac_sim import detection  # noqa: E402
from isac_sim.geometry import setup_geometry  # noqa: E402
from isac_sim.runner import driver_args, make_cfg, run_all, save_npz  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    args = driver_args(__doc__, os.path.join(HERE, 'results', 'detection_sweep.npz'))
    cfg = make_cfg(detection.DRIVER_DEFAULTS, args.set)
    assert cfg.adversary == 1, 'The attack treats user 1 as the adversary.'
    if not detection.HAS_ZOU:
        print('isac_sim/zou.py not found: skipping the MI-minimizing baseline.')
    geo = setup_geometry(cfg)
    results, elapsed, ckpt = run_all(detection.one_run, cfg, geo, list(range(1, cfg.num_runs + 1)),
                                     args.out, args.workers)
    S = detection.collect(results, cfg)
    S['elapsed'] = elapsed
    save_npz(args.out, S, cfg)
    print(f"Saved {args.out} ({S['num_ok']} of {cfg.num_runs} runs kept, {elapsed:.0f} s; "
          f"checkpoints in {ckpt}).")


if __name__ == '__main__':
    main()
