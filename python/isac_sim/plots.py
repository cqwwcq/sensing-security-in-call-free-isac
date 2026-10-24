"""Figures of the paper and of the response letter.

Python port of common/paper_fig_style.m and the plot_*.m scripts:
  fig2     Fig. 2, localization success vs eta     (plot_localization_success.m)
  fig3     Fig. 3, sensing SNR vs eta               (plot_sensing_snr.m)
  fig4     Fig. 4, outage vs CSI NMSE               (plot_outage.m)
  pd       letter, Reviewer 3 C2, Pd vs RCS change  (plot_detection_pd.m)
  mh       letter, Reviewer 2 C5, multi-hypothesis  (plot_multi_hypothesis.m)
  scatter  letter, Reviewer 1 C3, estimated locations (plot_location_scatter.m)
"""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import chi2, ncx2  # noqa: E402

STY = {'blue': (0.00, 0.45, 0.74), 'green': (0.47, 0.67, 0.19),
       'purple': (0.49, 0.18, 0.56), 'orange': (0.85, 0.33, 0.10),
       'gray': (0.40, 0.40, 0.40), 'lw': 1.2, 'ms': 5,
       'size_cm': (12.8, 10.2)}       # about the size of the published PDFs


def _style():
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'Liberation Serif', 'Nimbus Roman',
                       'DejaVu Serif'],
        'mathtext.fontset': 'stix', 'font.size': 9, 'legend.fontsize': 8,
        'axes.grid': True, 'grid.linestyle': ':', 'grid.alpha': 0.4, 'pdf.fonttype': 42})


def new_figure(width_factor=1.0, ncols=1):
    _style()
    w, h = STY['size_cm']
    fig, ax = plt.subplots(1, ncols, figsize=(width_factor * w / 2.54, h / 2.54))
    return fig, ax


def export(fig, path, legend=True):
    if legend:
        for ax in fig.axes[:1]:
            ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f'Wrote {path}')


def fig2(O, path):
    fig, ax = new_figure()
    xmax = 0.9
    ax.plot([0, xmax], [100 * O['p_noan']] * 2, '--', color=STY['blue'], lw=STY['lw'], label='No AN')
    ax.plot(O['eta'], 100 * O['p_null'], '-s', color=STY['green'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Null-space-only Gaussian AN')
    if 'p_zou' in O:
        ax.plot(O['eta_zou'], 100 * O['p_zou'], '-d', color=STY['purple'], lw=STY['lw'],
                ms=STY['ms'], mfc='none', label='MI-minimizing AN (Zou et al.)')
    ax.plot(O['eta'], 100 * O['p_prop'], '-o', color=STY['orange'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Proposed fake-peak AN')
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, max(45, 100 * max(np.max(O['p_noan']), np.max(O['p_null']),
                                     np.max(O['p_prop'])) + 5))
    ax.set_xlabel(r'Realized AN-Power Ratio $\eta$')
    ax.set_ylabel('Target-Localization Success Rate (%)')
    export(fig, path)


def fig3(O, path):
    fig, ax = new_figure()
    xmax = 0.9
    ax.plot([0, xmax], [O['snr_noan_mean']] * 2, '--', color=STY['blue'], lw=STY['lw'],
            label='No AN (reference)')
    ax.plot(O['eta'], O['snr_null_mean'], '-s', color=STY['green'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Null-space-only Gaussian AN')
    if 'snr_zou_mean' in O:
        ax.plot(O['eta_zou'], O['snr_zou_mean'], '-d', color=STY['purple'], lw=STY['lw'],
                ms=STY['ms'], mfc='none', label='MI-minimizing AN (Zou et al.)')
    ax.plot(O['eta'], O['snr_prop_mean'], '-o', color=STY['orange'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Proposed fake-peak AN')
    ax.plot(O['eta'], O['snrz_prop_mean'], '--o', color=STY['orange'], lw=STY['lw'],
            ms=STY['ms'], mfc='white', label='Proposed, AN echo exploited')
    ax.set_xlim(0, xmax)
    ax.set_xlabel(r'Realized AN-Power Ratio $\eta$')
    ax.set_ylabel('Target Sensing SNR (dB)')
    export(fig, path)


def fig4(O, path):
    fig, ax = new_figure()
    x = O['nmse']
    ax.plot(x, 100 * O['m_noan'], '-+', color=STY['green'], lw=STY['lw'], ms=STY['ms'], label='No AN')
    ax.plot(x, 100 * O['m_null'], '-d', color=STY['purple'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Null-space-only Gaussian AN')
    ax.plot(x, 100 * O['m_nom'], '-^', color=STY['orange'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Proposed nominal')
    ax.plot(x, 100 * O['m_rob'], '-s', color=STY['blue'], lw=STY['lw'], ms=STY['ms'],
            mfc='none', label='Proposed robust')
    if 'm_rnoan' in O:
        ax.plot(x, 100 * O['m_rnoan'], '--x', color=STY['gray'], lw=STY['lw'], ms=STY['ms'],
                label='Robust no AN')
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(0, 100)
    ax.set_xticks(np.arange(x[0], x[-1] + 1, 4))
    ax.set_xlabel('CSI NMSE (dB)')
    ax.set_ylabel('Communication Outage Probability (%)')
    export(fig, path)


def pd_figure(O, cfg, path, rng=None):
    """Monte Carlo Pd (markers) against the closed-form model (lines)."""
    rng = np.random.default_rng(1) if rng is None else rng
    j = int(np.argmin(np.abs(O['eta'] - 0.5)))
    off = np.asarray(cfg['pd_offset_dB'], float)
    Ns = cfg['Nue'] + 1
    df = 2 * cfg['Nrx'] * Ns
    tau = chi2.isf(cfg['Pfa'], df)
    G = rng.exponential(1.0, (50000, cfg['Nrx'])).sum(axis=1)   # sum of Nrx unit exponentials

    def model(g_db, N):
        return np.array([np.mean(ncx2.sf(tau, df, 2 * N * 10 ** ((g_db + o) / 10) * G)) for o in off])

    fig, ax = new_figure()
    mk = ['o', 's']
    for b, N in enumerate(cfg['pd_N_list']):
        ax.plot(off, O['pd_noan_mean'][:, b], mk[b], color=STY['blue'], ms=STY['ms'], mfc='none',
                label=f'No AN, N = {N}')
        ax.plot(off, model(O['snr_noan_mean'], N), '-', color=STY['blue'], lw=STY['lw'])
        ax.plot(off, O['pd_prop_mean'][j, :, b], mk[b], color=STY['orange'], ms=STY['ms'],
                mfc='none', label=rf'Proposed ($\eta$ = {O["eta"][j]:.2f}), N = {N}')
        ax.plot(off, model(O['snr_prop_mean'][j], N), '-', color=STY['orange'], lw=STY['lw'])
    ax.set_ylim(0, 1)
    ax.set_xlabel('RCS Change Relative to Table I (dB)')
    ax.set_ylabel(rf'Detection Probability ($P_{{fa}} = 10^{{{round(np.log10(cfg["Pfa"]))}}}$)')
    export(fig, path)


def mh_figure(O, path):
    rules = [('p', 'p_noan', STY['orange'], '-o', 'Dominant peak (paper)'),
             ('p_mhv', 'p_mhv_noan', STY['green'], '-s', 'Two peaks, more votes'),
             ('p_mhr', 'p_mhr_noan', STY['purple'], '-d', 'Two peaks, random choice'),
             ('p_mha', 'p_mha_noan', STY['blue'], '-^', 'Two peaks, mechanism-aware'),
             ('p_mhs', 'p_mhs_noan', STY['gray'], '-x', 'Two guesses')]
    fig, ax = new_figure()
    xmax = 0.9
    for key, key0, col, ls, lab in rules:
        ax.plot(O['eta'], 100 * O[f'{key}_prop'], ls, color=col, lw=STY['lw'], ms=STY['ms'],
                mfc='none', label=lab)
        ax.plot([0, xmax], [100 * O[key0]] * 2, '--', color=col, lw=0.8)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 100)
    ax.set_xlabel(r'Realized AN-Power Ratio $\eta$')
    ax.set_ylabel('Target-Localization Success Rate (%)')
    export(fig, path)


def scatter_figure(S, cfg, path):
    """Estimated locations without AN and with the proposed design at eta ~ 0.5.

    Estimates are cell centers, so each marker's area is proportional to the
    number of realizations at that cell.
    """
    eta = np.mean(S['eta_prop'], axis=-1)
    j = int(np.argmin(np.abs(eta - 0.5)))
    est = [S['est_noan'].T, S['est_prop'][j].T]
    titles = ['No AN', rf'Proposed, $\eta$ = {eta[j]:.2f}']
    colors = [STY['blue'], STY['orange']]
    A = cfg['area']
    fig, axes = new_figure(width_factor=2.0, ncols=2)
    for s, ax in enumerate(axes):
        loc, cnt = np.unique(est[s], axis=0, return_counts=True)
        ax.scatter(loc[:, 0], loc[:, 1], s=8 + 400 * cnt / cnt.max(), color=colors[s], alpha=0.6,
                   label='Estimates')
        ax.plot(*cfg['PosTarget'], 'kx', ms=10, mew=1.5, label='Target')
        ax.plot(*cfg['PosFake'], 'k+', ms=10, mew=1.5, label='Virtual location')
        tx = np.asarray(cfg['PosTx'])
        rx = np.asarray(cfg['PosRx'])
        ax.plot(tx[:, 0], tx[:, 1], 'k^', ms=6, mfc='none', label='Transmit APs')
        ax.plot(rx[:, 0], rx[:, 1], 'ks', ms=6, mfc='none', label='Receive APs')
        ax.set_aspect('equal')
        ax.set_xlim(A[0] - 10, A[1] + 10)
        ax.set_ylim(A[2] - 10, A[3] + 10)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_title(titles[s], fontweight='normal')
    axes[0].legend(loc='lower left', fontsize=6, framealpha=0.85, markerscale=0.6)
    export(fig, path, legend=False)
