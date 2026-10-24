"""Legitimate-side metrics.

Python port of common/user_sinr.m, common/sensing_metrics.m,
common/detection_pd.m, common/draw_nullonly_an.m and common/wilson_ci.m.
Beam layout: w_i is Ntx x Mtx x Nue, w_t and z are Ntx x Mtx, u is Nue x Mue.
"""
from types import SimpleNamespace

import numpy as np
from scipy.stats import chi2, ncx2

from .config import dbm2w
from .receiver import stream_responses


def user_sinr(cfg, H, u, w_i, w_t, z):
    """Post-combining SINR (5) of every user, linear, shape (Nue,).

    H is the designed-for or the actual channel (Nue x Ntx x Mue x Mtx).
    The interference includes the other users' streams, the sensing stream
    and the AN.
    """
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)
    q = stream_responses(H, w_i, w_t, z)              # Nue x Mue x (Nue+2)
    P = np.abs(np.einsum('im,iml->il', u.conj(), q)) ** 2
    idx = np.arange(cfg.Nue)
    sig = P[idx, idx]
    return sig / (P.sum(axis=1) - sig + sigma_ue)


def sensing_c(cfg, geo):
    """Sensing-objective weights c_k without the noise (as in the solver)."""
    return (cfg.sigmarcs ** 2 * geo.betark).sum(axis=1) / cfg.Nrx


def sensing_metrics(cfg, geo, w_i, w_t, z):
    """Sensing SNR, AN-echo term, fake-peak terms and eta.

    m.gamma_t    sensing SNR (2), linear, AN echo discarded
    m.gamma_t_z  sensing SNR when the known AN echo is also exploited,
                 gamma_t + sum_k c_k T_k (Reviewer 2, Comment 1)
    m.S, m.T, m.F  per-AP |a^H(phi_k) W_k|^2, |a^H(phi_k) z_k|^2 and
                 |a^H(phi_fake,k) z_k|^2 of Section III-B
    m.eta        realized AN-power ratio
    m.V          Nrx x Ntx x (Nue+1) echo vectors sqrt(sigma_rcs^2 g_{r,k})
                 a^H(phi_k) W_k, used by detection_pd
    """
    sigma_rx = dbm2w(cfg.sigma_rx_dBm)
    c = sensing_c(cfg, geo)
    W = np.concatenate([w_i, w_t[:, :, None]], axis=2)          # Ntx x Mtx x (Nue+1)
    aW = np.einsum('km,kml->kl', geo.a_true.conj(), W)            # Ntx x (Nue+1)
    m = SimpleNamespace()
    m.S = np.sum(np.abs(aW) ** 2, axis=1)
    m.T = np.abs(np.einsum('km,km->k', geo.a_true.conj(), z)) ** 2
    m.F = np.abs(np.einsum('km,km->k', geo.a_fake.conj(), z)) ** 2
    amp = np.sqrt(cfg.sigmarcs ** 2 * geo.betark).T               # Nrx x Ntx
    m.V = amp[:, :, None] * aW[None, :, :]
    p_tot = np.sum(np.abs(W) ** 2) + np.sum(np.abs(z) ** 2)
    m.gamma_t = np.sum(c * m.S) / sigma_rx
    m.gamma_t_z = np.sum(c * (m.S + m.T)) / sigma_rx
    m.eta = np.sum(np.abs(z) ** 2) / p_tot
    return m


def detection_pd(cfg, V, offset_dB, N_list, rng):
    """Detection probability of the GLRT for given transmit beams.

    Returns a len(offset_dB) x len(N_list) array.  V (Nrx x Ntx x Ns) holds
    the echo vectors from sensing_metrics.  The bistatic coefficients
    alpha_{r,k} ~ CN(0, 1) are drawn independently per transmit-receive pair
    (Swerling I, random echo phases), so the per-receive-AP SNR is
    gamma_r = ||sum_k alpha_{r,k} V[r,k,:]||^2 / sigma_rx^2, whose mean is the
    per-receive-AP SNR behind (2).

    Each receive AP matched-filters the Ns = Nue+1 known streams, and the
    GLRT sums the output energies over streams and receive APs.  The
    normalized statistic is chi-square with 2 Nrx Ns degrees of freedom
    under H0 and noncentral chi-square with noncentrality 2 N sum_r gamma_r
    under H1.  offset_dB scales all echoes, e.g. to emulate a smaller RCS.
    """
    Nrx, Ntx, Ns = V.shape
    sigma_rx = dbm2w(cfg.sigma_rx_dBm)
    df = 2 * Nrx * Ns
    tau = chi2.isf(cfg.Pfa, df)
    D = cfg.pd_alpha_draws
    gsum = np.zeros(D)
    for r in range(Nrx):
        alpha = (rng.standard_normal((D, Ntx)) + 1j * rng.standard_normal((D, Ntx))) / np.sqrt(2)
        cr = alpha @ V[r]                                         # D x Ns
        gsum += np.sum(np.abs(cr) ** 2, axis=1) / sigma_rx
    pd = np.zeros((len(offset_dB), len(N_list)))
    for a, off in enumerate(offset_dB):
        for b, N in enumerate(N_list):
            lam = 2 * N * 10 ** (off / 10) * gsum
            pd[a, b] = np.mean(ncx2.sf(tau, df, lam))
    return pd


def effective_channel_matrix(H, u):
    """Rows u_i^H [H_{i,1}, ..., H_{i,Ntx}]: Nue x Ntx x Mtx."""
    return np.einsum('im,ikmn->ikn', u.conj(), H)


def draw_nullonly_an(cfg, Hhat, u, w_i, w_t, p_an, rng):
    """Gaussian AN projected onto the effective-channel null space.

    Draws a complex Gaussian vector, projects the stacked vector onto N(G_eff)
    of (6), where row i of G_eff is u_i^H [Hhat_{i,1}, ..., Hhat_{i,Ntx}], and
    scales it with one common factor so that the null is preserved.  The
    factor aims at the per-AP AN powers p_an (length Ntx, W) in total, but is
    reduced if an AP would exceed its budget P_T left after w_i and w_t.
    Returns z, Ntx x Mtx.
    """
    Ntx, Mtx = cfg.Ntx, cfg.Mtx
    P_T_lin = dbm2w(cfg.P_T)
    G = effective_channel_matrix(Hhat, u).reshape(cfg.Nue, Ntx * Mtx)
    r = (rng.standard_normal(Ntx * Mtx) + 1j * rng.standard_normal(Ntx * Mtx)) / np.sqrt(2)
    zbar = r - G.conj().T @ (np.linalg.pinv(G @ G.conj().T) @ (G @ r))
    zk = zbar.reshape(Ntx, Mtx)
    pk = np.sum(np.abs(zk) ** 2, axis=1)
    p_beam = np.sum(np.abs(w_i) ** 2, axis=(1, 2)) + np.sum(np.abs(w_t) ** 2, axis=1)
    p_left = np.maximum(P_T_lin - p_beam, 0.0)
    s2 = min(np.sum(p_an) / np.sum(pk), np.min(p_left / pk))
    return np.sqrt(s2) * zk


def wilson_ci(p, n, z=1.96):
    """Wilson score interval for a binomial proportion p over n trials."""
    p = np.asarray(p, dtype=float)
    den = 1 + z ** 2 / n
    mid = (p + z ** 2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return np.maximum(mid - half, 0.0), np.minimum(mid + half, 1.0)
