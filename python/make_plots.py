"""Figures of the paper and of the response letter (port of the plot_*.m scripts).

By default the figures go to results/figures/.  With --paper, Figs. 2-4 are
written to the paper's figure folder instead, replacing the published PDFs
(as the MATLAB plot scripts do).

Examples:
    python make_plots.py                      # all figures whose inputs exist
    python make_plots.py fig2 fig3 --paper
"""
import argparse
import os

from isac_sim import plots
from isac_sim.runner import load_npz

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER_FIG = os.path.join(HERE, '..', 'Fake-Peak Artificial Noise for Sensing-Location Privacy '
                         'in Cell-Free Massive MIMO ISAC', 'figure')
ALL = ['fig2', 'fig3', 'fig4', 'pd', 'mh', 'scatter']


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('which', nargs='*', choices=ALL + [[]], default=[],
                    help='figures to make (default: all)')
    ap.add_argument('--results', default=os.path.join(HERE, 'results'))
    ap.add_argument('--figdir', default=os.path.join(HERE, 'results', 'figures'))
    ap.add_argument('--paper', action='store_true', help='write Figs. 2-4 into the paper folder')
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)
    res = lambda name: os.path.join(args.results, name)  # noqa: E731
    paper = {'fig2': 'localization_success_tvt.pdf', 'fig3': 'sensing_snr_final.pdf',
             'fig4': 'outage_vs_nmse.pdf'}

    def target(key, default):
        if args.paper and key in paper:
            return os.path.join(PAPER_FIG, paper[key])
        return os.path.join(args.figdir, paper.get(key, default))

    for key in args.which or ALL:
        need = {'fig4': 'outage_results.npz', 'scatter': 'detection_sweep.npz'}.get(
            key, 'detection_results.npz')
        if not os.path.exists(res(need)):
            print(f'{key}: {res(need)} not found, skipped')
            continue
        D = load_npz(res(need))
        cfg = D.pop('cfg')
        if key == 'fig2':
            plots.fig2(D, target(key, ''))
        elif key == 'fig3':
            plots.fig3(D, target(key, ''))
        elif key == 'fig4':
            plots.fig4(D, target(key, ''))
        elif key == 'pd':
            plots.pd_figure(D, cfg, target(key, 'detection_pd_response.pdf'))
        elif key == 'mh':
            plots.mh_figure(D, target(key, 'multi_hypothesis_response.pdf'))
        elif key == 'scatter':
            plots.scatter_figure(D, cfg, target(key, 'location_scatter_response.pdf'))


if __name__ == '__main__':
    main()
