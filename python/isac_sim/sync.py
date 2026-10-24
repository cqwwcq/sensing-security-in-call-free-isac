"""Robustness of the effective-channel null to AP phase errors.

Python port of fig.9.outage/run_sync_sweep.m (Reviewer 2, Comment 4).

The beams are designed on the LMMSE estimates of Table I (no-AN baseline and
proposed nominal design, rho = cfg.rou).  On the actual channels, the signal
of AP k is rotated by a residual carrier phase psi_k ~ N(0, s^2), i.i.d. over
APs and cfg.num_phase draws, for each s in cfg.sync_std_deg.  Saved:
  m_out_noan, m_out_nom   outage probability Pr{min_i gamma_i < Gamma}
  c_out_*                 95 % half-widths over the channel realizations
  m_anr                   mean residual AN power at the combiner outputs
                          relative to the noise, dB
  m_anr_theory            the same from the closed form
                          E|sum_k e^{j psi_k} c_k|^2 = e^{-s^2} |sum_k c_k|^2
                          + (1 - e^{-s^2}) sum_k |c_k|^2, c_k = u_i^H H_{i,k} z_k
The sensing SNR (2) is unaffected, since it combines the transmit APs
non-coherently, and so is the attack, which works per AP.
"""
import numpy as np

from .channels import gen_channels
from .config import db2lin, dbm2w
from .detection import fake_mask, solve_design

DRIVER_DEFAULTS = {'rou': 1.3, 'sync_std_deg': [0, 1, 2, 5, 10, 20], 'num_runs': 100,
                   'num_phase': 1000, 'base_seed': 20260825}


def outage(cfg, H, rot, u, w_i, w_t, z):
    """Outage with per-AP phase rotations; also the residual AN statistics."""
    Nue = cfg.Nue
    Gamma = db2lin(cfg.Gamma_comm)
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)
    X = np.concatenate([w_i, w_t[:, :, None], z[:, :, None]], axis=2)   # Ntx x Mtx x (Nue+2)
    fail = np.zeros(rot.shape[0], bool)
    an_power = an_coh = an_inc = 0.0
    for i in range(Nue):
        # Per-AP contributions of every stream at user i's combiner output.
        g = np.einsum('m,kmn->kn', u[i].conj(), H[i])                    # Ntx x Mtx
        C = np.einsum('kn,knl->kl', g, X)                                # Ntx x (Nue+2)
        P = np.abs(rot @ C) ** 2
        sinr = P[:, i] / (P.sum(axis=1) - P[:, i] + sigma_ue)
        fail |= sinr < Gamma
        an_power += np.mean(P[:, -1]) / Nue
        an_coh += np.abs(np.sum(C[:, -1])) ** 2 / Nue
        an_inc += np.sum(np.abs(C[:, -1]) ** 2) / Nue
    return float(np.mean(fail)), an_power, an_coh, an_inc


def one_run(cfg, geo, r):
    rng = np.random.default_rng(cfg.base_seed + r)
    ch = gen_channels(cfg, geo, rng, 'pilot')
    s0 = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou, {'an_mode': 'none'})
    z0 = np.zeros((cfg.Ntx, cfg.Mtx), dtype=complex)
    s1 = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou, {
        'an_mode': 'fake', 'sym_null': cfg.sym_null_prop, 'fake_tleak': cfg.fake_tleak,
        'fake_mask': fake_mask(cfg, geo), 'fake_upper': cfg.fake_upper})
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)
    ns = len(cfg.sync_std_deg)
    o_noan = np.zeros(ns)
    o_nom = np.zeros(ns)
    anr = np.zeros(ns)
    anr_th = np.zeros(ns)
    for s, sd_deg in enumerate(cfg.sync_std_deg):
        sd = sd_deg * np.pi / 180
        rot = np.exp(1j * sd * rng.standard_normal((cfg.num_phase, cfg.Ntx)))
        o_noan[s] = outage(cfg, ch.H, rot, s0.u, s0.w_i, s0.w_t, z0)[0]
        o_nom[s], an_p, an_c, an_i = outage(cfg, ch.H, rot, s1.u, s1.w_i, s1.w_t, s1.z)
        anr[s] = an_p / sigma_ue
        anr_th[s] = (np.exp(-sd ** 2) * an_c + (1 - np.exp(-sd ** 2)) * an_i) / sigma_ue
    return {'o_noan': o_noan, 'o_nom': o_nom, 'anr': anr, 'anr_th': anr_th,
            'ok': bool(s0.success and s1.success)}


def collect(results, cfg):
    ok = np.array([r['ok'] for r in results])
    n = int(ok.sum())
    get = lambda key: np.stack([r[key] for r in results])[ok]   # noqa: E731
    S = {'sync_std_deg': np.asarray(cfg.sync_std_deg, float), 'n_kept': n}
    for key, name in (('o_noan', 'noan'), ('o_nom', 'nom')):
        x = get(key)
        S[f'm_out_{name}'] = x.mean(axis=0)
        S[f'c_out_{name}'] = 1.96 * x.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.full(x.shape[1], np.nan)
    S['m_anr'] = 10 * np.log10(get('anr').mean(axis=0))
    S['m_anr_theory'] = 10 * np.log10(get('anr_th').mean(axis=0))
    return S


def print_results(S):
    for s, sd in enumerate(S['sync_std_deg']):
        print(f"  phase std {sd:4.1f} deg: outage no AN {100 * S['m_out_noan'][s]:.2f} %, "
              f"proposed {100 * S['m_out_nom'][s]:.2f} %; residual AN/noise "
              f"{S['m_anr'][s]:.1f} dB (theory {S['m_anr_theory'][s]:.1f} dB)")
