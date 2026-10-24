"""Angles, steering vectors and large-scale gains of Section II.

Python port of common/compute_angle_dist.m, common/beamsteering.m and
common/setup_geometry.m.
"""
from types import SimpleNamespace

import numpy as np


def compute_angle_dist(pos1, pos2):
    """Angle (from pos1 to pos2, rad) and distance of every pair, M x N."""
    pos1 = np.atleast_2d(np.asarray(pos1, dtype=float))
    pos2 = np.atleast_2d(np.asarray(pos2, dtype=float))
    d = pos2[None, :, :] - pos1[:, None, :]
    dist = np.sqrt(np.sum(d ** 2, axis=2))
    ang = np.arctan2(d[:, :, 1], d[:, :, 0])
    return ang, dist


def beamsteering(theta, num_ant, spacing=0.5):
    """Steering vectors a(theta) with the antenna index on the last axis.

    a(theta) = [1, e^{j pi sin theta}, ..., e^{j pi (M-1) sin theta}] for
    half-wavelength spacing; the output has shape theta.shape + (num_ant,).
    """
    theta = np.asarray(theta, dtype=float)
    k = np.arange(num_ant)
    return np.exp(1j * 2 * np.pi * spacing * k * np.sin(theta)[..., None])


def setup_geometry(cfg):
    """Geometry of the scenario.

    geo.phi, geo.phi_fake   true and fake directions phi_k, phi_fake,k (Ntx,)
    geo.a_true, geo.a_fake  steering vectors a(phi_k), a(phi_fake,k), Ntx x Mtx
    geo.betark              Ntx x Nrx two-way path gains g_{r,k}
    geo.beta_ue             Nue x Ntx large-scale fading beta_{i,k}
    geo.resolvable          APs whose true and fake directions are at least
                            one Rayleigh width apart,
                            |sin(phi_k) - sin(phi_fake,k)| >= 2/Mtx
    """
    phi, d_tk = compute_angle_dist(cfg.PosTx, cfg.PosTarget)
    phi_fake, _ = compute_angle_dist(cfg.PosTx, cfg.PosFake)
    _, d_tr = compute_angle_dist(cfg.PosRx, cfg.PosTarget)
    _, d_ik = compute_angle_dist(cfg.PosUE, cfg.PosTx)

    geo = SimpleNamespace()
    geo.phi = phi[:, 0]
    geo.phi_fake = phi_fake[:, 0]
    geo.a_true = beamsteering(geo.phi, cfg.Mtx)
    geo.a_fake = beamsteering(geo.phi_fake, cfg.Mtx)
    d_tk = d_tk[:, 0]
    d_tr = d_tr[:, 0]
    # g_{r,k} = lambda^2 / ((4 pi)^3 d_{t,k}^2 d_{t,r}^2), stored as betark[k, r]
    geo.betark = cfg.lambda_c ** 2 / ((4 * np.pi) ** 3 * np.outer(d_tk ** 2, d_tr ** 2))
    geo.beta_ue = d_ik ** (-cfg.pathloss_exp)
    geo.resolvable = np.abs(np.sin(geo.phi) - np.sin(geo.phi_fake)) >= 2.0 / cfg.Mtx
    geo.d_tk = d_tk
    geo.d_tr = d_tr
    return geo
