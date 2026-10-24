"""Aggregate a detection sweep into the plotted values and print the numbers
quoted in the paper and the letter (port of make_detection_results.m).

Example:
    python make_detection_results.py results/detection_sweep.npz results/detection_results.npz
"""
import argparse
import json
import os

import numpy as np

from isac_sim import detection
from isac_sim.runner import load_npz

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('sweep', nargs='?', default=os.path.join(HERE, 'results', 'detection_sweep.npz'))
    ap.add_argument('out', nargs='?', default=os.path.join(HERE, 'results', 'detection_results.npz'))
    args = ap.parse_args()
    S = load_npz(args.sweep)
    cfg = S.pop('cfg')
    O = detection.make_results(S, cfg)
    np.savez_compressed(args.out, cfg_json=np.array(json.dumps(cfg)), **O)
    print(f'Saved {args.out}')
    detection.print_results(O, cfg)


if __name__ == '__main__':
    main()
