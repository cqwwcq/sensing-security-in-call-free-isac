"""AP-user channels and their estimates.

Python port of common/gen_channels.m.  Channels are stored as arrays of shape
Nue x Ntx x Mue x Mtx, i.e. H[i, k] is the Mue x Mtx channel from AP k to
user i (MATLAB: H(:, cols_k, i)).
"""
from types import SimpleNamespace

import numpy as np

from .config import dbm2w


def crandn(rng, *shape):
    """Circularly symmetric complex Gaussian samples with unit variance."""
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def gen_channels(cfg, geo, rng, mode='pilot', nmse=None):
    """Channels, LMMSE estimates and error variances.

    mode 'pilot': draws H_{i,k} = sqrt(beta_ik) G_ik and forms the LMMSE
    estimate from tau_p-length orthogonal pilots of per-antenna power
    P_p / M_UE (Section II), so that the error variance is epsilon_ik^2 *
    beta_ik with epsilon_ik^2 = sigma_p^2 / (tau_p p beta_ik + sigma_p^2).

    mode 'nmse': draws the estimate directly with a common link NMSE
    (linear), as in the Fig. 4 sweep: Hhat has variance (1 - nmse) beta_ik
    and the error variance is nmse * beta_ik.  The actual channel is not
    drawn (ch.H is None); the caller adds error draws with ch.beta_err.

    ch.beta_err is Ntx x Nue, the layout of beta_error in the robust solver.
    """
    Ntx, Nue, Mtx, Mue = cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Mue
    ch = SimpleNamespace()
    ch.H = np.zeros((Nue, Ntx, Mue, Mtx), dtype=complex)
    ch.Hhat = np.zeros((Nue, Ntx, Mue, Mtx), dtype=complex)
    ch.beta_err = np.zeros((Ntx, Nue))

    if mode == 'pilot':
        p = dbm2w(cfg.P_p_dBm) / Mue
        s2 = dbm2w(cfg.sigma_p_dBm)
        tp = cfg.tau_p
        for i in range(Nue):
            for k in range(Ntx):
                b = geo.beta_ue[i, k]
                H = np.sqrt(b) * crandn(rng, Mue, Mtx)
                Y = H + np.sqrt(s2 / (tp * p)) * crandn(rng, Mue, Mtx)
                ch.H[i, k] = H
                ch.Hhat[i, k] = tp * p * b / (tp * p * b + s2) * Y
                ch.beta_err[k, i] = b * s2 / (tp * p * b + s2)
    elif mode == 'nmse':
        for i in range(Nue):
            for k in range(Ntx):
                b = geo.beta_ue[i, k]
                ch.Hhat[i, k] = np.sqrt((1 - nmse) * b) * crandn(rng, Mue, Mtx)
                ch.beta_err[k, i] = nmse * b
        ch.H = None
    else:
        raise ValueError(f'Unknown channel mode "{mode}".')
    return ch
