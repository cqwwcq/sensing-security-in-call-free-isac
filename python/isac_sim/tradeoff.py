"""What the AN power costs in communication (Reviewer 3, Comment 3).

Python port of fig.7.detection/run_comm_tradeoff.m.

For each channel realization, the proposed design is solved at
rho = cfg.rho_ref, the operating point eta ~ 0.5 of Figs. 2-3.  The no-AN
design is then given the same sensing SNR and as much communication as
possible: the largest SINR target Gamma for which the no-AN solver still
reaches the proposed design's sensing SNR is found by bisection.  Saved, per
kept realization:
  gmax_noan_dB   largest common SINR target of the no-AN alternative
  rate_noan      its sum rate, bit/s/Hz, at that target
  rate_prop      sum rate of the proposed design
  snr_prop_dB    sensing SNR of the proposed design (2)
  eta_prop       realized AN-power ratio of the proposed design
"""
import numpy as np

from .channels import gen_channels
from .detection import prop_opts, solve_design
from .metrics import sensing_metrics, user_sinr

DRIVER_DEFAULTS = {'rho_ref': 0.75, 'gamma_hi_dB': 40.0, 'bisect_iter': 7,
                   'num_runs': 100, 'base_seed': 20260826}


def one_run(cfg, geo, r):
    rng = np.random.default_rng(cfg.base_seed + r)
    ch = gen_channels(cfg, geo, rng, 'pilot')
    s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rho_ref, prop_opts(cfg, geo))
    m = sensing_metrics(cfg, geo, s.w_i, s.w_t, s.z)
    out = {'ok': bool(s.success), 'snr_prop_dB': 10 * np.log10(m.gamma_t), 'eta_prop': m.eta,
           'rate_prop': float(np.sum(np.log2(1 + user_sinr(cfg, ch.Hhat, s.u, s.w_i, s.w_t, s.z))))}

    # Bisection on the SINR target of the no-AN design at the same sensing SNR.
    lo, hi = cfg.Gamma_comm, cfg.gamma_hi_dB
    best = None
    for _ in range(cfg.bisect_iter):
        mid = (lo + hi) / 2
        s0 = solve_design(cfg, geo, ch.Hhat, cfg.P_T, 0.0, {'an_mode': 'none'}, Gamma_dB=mid)
        m0 = sensing_metrics(cfg, geo, s0.w_i, s0.w_t, s0.z)
        if s0.success and m0.gamma_t >= m.gamma_t:
            lo = mid
            best = s0
        else:
            hi = mid
    if best is None:
        # The bisection never moved up: fall back to the 10 dB target.
        best = solve_design(cfg, geo, ch.Hhat, cfg.P_T, 0.0, {'an_mode': 'none'}, Gamma_dB=lo)
        out['ok'] = out['ok'] and bool(best.success)
    out['gmax_noan_dB'] = lo
    zero = np.zeros((cfg.Ntx, cfg.Mtx), dtype=complex)
    out['rate_noan'] = float(np.sum(np.log2(1 + user_sinr(cfg, ch.Hhat, best.u, best.w_i,
                                                          best.w_t, zero))))
    return out


def collect(results, cfg):
    keep = [r for r in results if r['ok']]
    S = {'n_kept': len(keep)}
    for f in ('gmax_noan_dB', 'rate_noan', 'rate_prop', 'snr_prop_dB', 'eta_prop'):
        S[f] = np.array([r[f] for r in keep])
    return S


def print_results(S, cfg):
    print(f"Proposed (rho = {cfg.rho_ref:.2f}): eta {np.mean(S['eta_prop']):.3f}, sensing SNR "
          f"{np.mean(S['snr_prop_dB']):.2f} dB, sum rate {np.mean(S['rate_prop']):.2f} bit/s/Hz")
    print(f"No-AN alternative at the same sensing SNR: common SINR {np.mean(S['gmax_noan_dB']):.2f} dB "
          f"(+{np.mean(S['gmax_noan_dB']) - cfg.Gamma_comm:.2f} dB), sum rate "
          f"{np.mean(S['rate_noan']):.2f} bit/s/Hz")
