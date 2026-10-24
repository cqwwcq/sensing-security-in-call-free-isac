"""Sweep behind Fig. 4: communication outage versus CSI NMSE.

Python port of fig.9.outage/run_outage_sweep.m.

For each NMSE and channel realization, the beams are designed on an estimate
with link NMSE nu (gen_channels 'nmse' mode) and the outage
Pr{min_i gamma_i^AN < Gamma_comm} is evaluated over cfg.num_err draws of the
actual channel H = Hhat + E with the designed combiners.  Designs:
  noan   no-AN baseline
  null   proposed nominal beams with the AN redrawn as Gaussian AN in the
         effective-channel null space at the same power
  nom    proposed nominal design (rho = cfg.rou)
  rob    proposed robust design (20)
  rnoan  robust no-AN baseline (Reviewer 1, Comment 5), if
         cfg.include_robust_noan
Saved: m_* (outage), c_* (95 % half-widths), m_sn* (sensing SNR, dB),
n_kept, n_drop, nmse.

Since ||u_i|| = 1 and E_{i,k} has i.i.d. CN(0, nmse*beta_ik) entries,
u_i^H E_{i,k} has i.i.d. CN(0, nmse*beta_ik) entries, so the error is drawn
directly on the effective channels.
"""
import numpy as np

from .channels import gen_channels
from .config import db2lin, dbm2w
from .detection import fake_mask, solve_design
from .metrics import draw_nullonly_an, effective_channel_matrix, sensing_metrics

DRIVER_DEFAULTS = {'rou': 1.3, 'nmse_db': list(range(-28, -3, 2)), 'num_runs': 300,
                   'num_err': 3000, 'base_seed': 20260824, 'include_robust_noan': True}
SNR_NAMES = {'noan': 'm_sna', 'null': 'm_snu', 'nom': 'm_sn', 'rob': 'm_sr', 'rnoan': 'm_srn'}


def designs(cfg):
    d = ['noan', 'null', 'nom', 'rob']
    if cfg.include_robust_noan:
        d.append('rnoan')
    return d


def draw_effective_errors(cfg, beta_err, rng):
    """num_err x Nue x (Ntx*Mtx) effective-channel errors u_i^H E_{i,k}."""
    D, Nue, Ntx, Mtx = cfg.num_err, cfg.Nue, cfg.Ntx, cfg.Mtx
    E = (rng.standard_normal((D, Nue, Ntx * Mtx)) + 1j * rng.standard_normal((D, Nue, Ntx * Mtx))) / np.sqrt(2)
    scale = np.repeat(np.sqrt(beta_err), Mtx, axis=0).T            # Nue x Ntx*Mtx
    return E * scale[None]


def evaluate(cfg, geo, ch, E, w_i, w_t, u, z):
    """Outage over the error draws with fixed combiners and beams."""
    Nue, Ntx, Mtx = cfg.Nue, cfg.Ntx, cfg.Mtx
    Gamma = db2lin(cfg.Gamma_comm)
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)
    B = np.concatenate([w_i, w_t[:, :, None], z[:, :, None]], axis=2).reshape(Ntx * Mtx, Nue + 2)
    G = effective_channel_matrix(ch.Hhat, u).reshape(Nue, Ntx * Mtx)
    fail = np.zeros(cfg.num_err, bool)
    for i in range(Nue):
        g = G[i][None, :] + E[:, i, :]
        P = np.abs(g @ B) ** 2
        sinr = P[:, i] / (P.sum(axis=1) - P[:, i] + sigma_ue)
        fail |= sinr < Gamma
    m = sensing_metrics(cfg, geo, w_i, w_t, z)
    return float(np.mean(fail)), float(10 * np.log10(m.gamma_t))


def one_run(cfg, geo, r):
    dlist = designs(cfg)
    nm, nd = len(cfg.nmse_db), len(dlist)
    out = np.full((nm, nd), np.nan)
    snr = np.full((nm, nd), np.nan)
    ok = np.zeros(nm, bool)
    prop = {'an_mode': 'fake', 'sym_null': cfg.sym_null_prop, 'fake_tleak': cfg.fake_tleak,
            'fake_mask': fake_mask(cfg, geo), 'fake_upper': cfg.fake_upper}
    for j, nmse_db in enumerate(cfg.nmse_db):
        # Same seed for every NMSE: the estimate differs only in its scaling.
        rng = np.random.default_rng(cfg.base_seed + r)
        ch = gen_channels(cfg, geo, rng, 'nmse', db2lin(nmse_db))
        E = draw_effective_errors(cfg, ch.beta_err, rng)
        okj = True
        nom = None
        for d, name in enumerate(dlist):
            if name == 'null':
                continue                    # filled in after 'nom'
            if name == 'noan':
                s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou, {'an_mode': 'none'})
            elif name == 'nom':
                s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou, prop)
                nom = s
            elif name == 'rob':
                s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou,
                                 dict(prop, beta_error=ch.beta_err))
            elif name == 'rnoan':
                s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, cfg.rou,
                                 {'an_mode': 'none', 'beta_error': ch.beta_err})
            okj &= s.success
            out[j, d], snr[j, d] = evaluate(cfg, geo, ch, E, s.w_i, s.w_t, s.u, s.z)
        if 'null' in dlist and nom is not None:
            d = dlist.index('null')
            p_an = np.sum(np.abs(nom.z) ** 2, axis=1)
            z = draw_nullonly_an(cfg, ch.Hhat, nom.u, nom.w_i, nom.w_t, p_an, rng)
            out[j, d], snr[j, d] = evaluate(cfg, geo, ch, E, nom.w_i, nom.w_t, nom.u, z)
        ok[j] = okj
    return {'out': out, 'snr': snr, 'ok': ok}


def collect(results, cfg):
    dlist = designs(cfg)
    out = np.stack([r['out'] for r in results])          # R x nm x nd
    snr = np.stack([r['snr'] for r in results])
    okr = np.stack([r['ok'] for r in results])            # R x nm
    R = len(results)
    S = {'nmse': np.asarray(cfg.nmse_db, float), 'n_kept': okr.sum(axis=0),
         'n_drop': R - okr.sum(axis=0)}
    for d, name in enumerate(dlist):
        m = np.zeros(len(cfg.nmse_db))
        c = np.zeros(len(cfg.nmse_db))
        s = np.zeros(len(cfg.nmse_db))
        for j in range(len(cfg.nmse_db)):
            x = out[okr[:, j], j, d]
            m[j] = np.mean(x)
            c[j] = 1.96 * np.std(x, ddof=1) / np.sqrt(x.size) if x.size > 1 else np.nan
            s[j] = np.mean(snr[okr[:, j], j, d])
        S[f'm_{name}'] = m
        S[f'c_{name}'] = c
        S[SNR_NAMES[name]] = s
    return S


def print_results(O):
    x = list(O['nmse'])
    for v in (-10, -4):
        if v not in x:
            continue
        j = x.index(v)
        line = (f"NMSE {v} dB: outage no AN {100 * O['m_noan'][j]:.2f} %, null-only "
                f"{100 * O['m_null'][j]:.2f} %, nominal {100 * O['m_nom'][j]:.2f} %, "
                f"robust {100 * O['m_rob'][j]:.2f} %")
        if 'm_rnoan' in O:
            line += f", robust no AN {100 * O['m_rnoan'][j]:.2f} %"
        print(line)
        line = (f"  sensing SNR no AN {O['m_sna'][j]:.2f} dB, nominal {O['m_sn'][j]:.2f} dB, "
                f"robust {O['m_sr'][j]:.2f} dB")
        if 'm_srn' in O:
            line += f", robust no AN {O['m_srn'][j]:.2f} dB"
        print(line)
