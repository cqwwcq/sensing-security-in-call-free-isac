"""Localization step of the reconstruction attack.

Python port of common/ray_cells.m, common/localize_target.m and
common/localize_multi.m.  Search cells are indexed 0-based by (ix, iy); a
linear cell index is ix + iy * nx (the column-major order of MATLAB, so that
ties are broken among the cells in the same order).
"""
from types import SimpleNamespace

import numpy as np

from .geometry import beamsteering


def _grid(cfg, cs=None):
    cs = cfg.cell_size if cs is None else cs
    A = cfg.area
    nx = int(round((A[1] - A[0]) / cs))
    ny = int(round((A[3] - A[2]) / cs))
    return cs, A, nx, ny


def _scan(cfg, R):
    """Beampatterns B_k(theta) = a^H(theta) R_k a(theta) on the angle grid."""
    theta = np.linspace(-np.pi, np.pi, cfg.n_ang)
    At = beamsteering(theta, cfg.Mtx)                          # n_ang x Mtx
    B = np.real(np.einsum('tm,kmn,tn->kt', At.conj(), R, At))   # Ntx x n_ang
    return theta, B


def ray_cells(p0, th, A, cs, nx, ny):
    """Linear indices of the search-area cells crossed by a ray from p0."""
    L = np.hypot(A[1] - A[0], A[3] - A[2])
    t = np.arange(0.0, L + 1e-12, cs / 20)
    P = np.asarray(p0, dtype=float)[None, :] + t[:, None] * np.array([np.cos(th), np.sin(th)])
    tol = 1e-9
    inside = (P[:, 0] >= A[0] - tol) & (P[:, 0] <= A[1] + tol) \
        & (P[:, 1] >= A[2] - tol) & (P[:, 1] <= A[3] + tol)
    P = P[inside]
    # A ray from a boundary AP that points outward leaves only the AP itself.
    if P.shape[0] < 2:
        return np.zeros(0, dtype=int)
    ix = np.clip(np.floor((P[:, 0] - A[0]) / cs).astype(int), 0, nx - 1)
    iy = np.clip(np.floor((P[:, 1] - A[2]) / cs).astype(int), 0, ny - 1)
    return np.unique(ix + iy * nx)


def _mirror_rays(cfg, k, th, A, cs, nx, ny, rng):
    """Rays kept for a peak at th: the ULA cannot tell th from pi - th."""
    cand = [th, np.angle(np.exp(1j * (np.pi - th)))]
    cells = [ray_cells(cfg.PosTx[k], c, A, cs, nx, ny) for c in cand]
    valid = [c.size > 0 for c in cells]
    if abs(cand[0] - cand[1]) < 1e-9:
        valid[1] = False
    pick = [c for c in range(2) if valid[c]]
    if len(pick) > 1 and cfg.mirror_rule == 'random':
        pick = [pick[rng.integers(len(pick))]]
    return cand, cells, pick


def _argmax_random(v, rng):
    best = np.flatnonzero(v == np.max(v))
    return int(best[rng.integers(best.size)])


def _target_cell(cfg, A, cs, nx, ny):
    tx = min(int(np.floor((cfg.PosTarget[0] - A[0]) / cs)), nx - 1)
    ty = min(int(np.floor((cfg.PosTarget[1] - A[2]) / cs)), ny - 1)
    return tx, ty


def localize_target(cfg, geo, R, rng, cell_size=None):
    """Dominant-peak line fusion of the reconstruction attack.

    R is Ntx x Mtx x Mtx, the sample covariances R_{xhat_k} of (3).  For each
    AP the dominant angle of B_k(theta) defines a ray from the AP; every cell
    crossed by a ray receives one vote, and the cell with the most votes is
    the estimate (ties broken at random).  The ULA cannot tell theta from
    pi - theta, so only candidates whose ray enters the search area are kept;
    cfg.mirror_rule = 'random' picks one of them, 'both' casts both rays.

    Returns (success, err, est, dom): success if the estimate lies in the
    target's cell, err the distance between the estimated cell center and the
    target, and per-AP diagnostics dom.theta_hat, dom.fake_dominant (the
    reconstructed fake-direction response exceeds the true-direction one)
    and dom.peak_class (1 fake direction, 2 true direction, 3 elsewhere;
    within half a main lobe, |sin(theta) - sin(phi)| <= 1/Mtx).
    """
    cs, A, nx, ny = _grid(cfg, cell_size)
    Ntx, Mtx = cfg.Ntx, cfg.Mtx
    theta, B = _scan(cfg, R)
    votes = np.zeros(nx * ny)
    dom = SimpleNamespace(theta_hat=np.zeros(Ntx), fake_dominant=np.zeros(Ntx, bool),
                          peak_class=np.zeros(Ntx, int))
    for k in range(Ntx):
        th = theta[int(np.argmax(B[k]))]
        if abs(np.sin(th) - np.sin(geo.phi[k])) <= 1 / Mtx:
            dom.peak_class[k] = 2
        elif abs(np.sin(th) - np.sin(geo.phi_fake[k])) <= 1 / Mtx:
            dom.peak_class[k] = 1
        else:
            dom.peak_class[k] = 3
        cand, cells, pick = _mirror_rays(cfg, k, th, A, cs, nx, ny, rng)
        for c in pick:
            votes[cells[c]] += 1
        if pick:
            dom.theta_hat[k] = cand[pick[0]]
        bt = np.real(geo.a_true[k].conj() @ R[k] @ geo.a_true[k])
        bf = np.real(geo.a_fake[k].conj() @ R[k] @ geo.a_fake[k])
        dom.fake_dominant[k] = bf > bt

    best = _argmax_random(votes, rng)
    ix, iy = best % nx, best // nx
    est = np.array([A[0] + (ix + 0.5) * cs, A[2] + (iy + 0.5) * cs])
    tx, ty = _target_cell(cfg, A, cs, nx, ny)
    success = (ix == tx) and (iy == ty)
    err = float(np.linalg.norm(est - cfg.PosTarget))
    return success, err, est, dom


def _two_peaks(B, theta, Mtx):
    """Strongest local maximum and the strongest one at least 2/Mtx away
    from it in sin(theta) (circular neighbours, as in localize_multi.m)."""
    prv = np.roll(B, 1)
    nxt = np.roll(B, -1)
    loc = np.flatnonzero((B >= prv) & (B >= nxt))
    loc = loc[np.argsort(-B[loc], kind='stable')]
    pk = [int(loc[0])]
    s1 = np.sin(theta[loc[0]])
    for q in loc[1:]:
        if abs(np.sin(theta[q]) - s1) >= 2 / Mtx:
            pk.append(int(q))
            break
    return pk


def localize_multi(cfg, geo, R, rng):
    """Multi-hypothesis reconstruction attack (Reviewer 2, Comment 5).

    Instead of the dominant peak only, the adversary keeps the two strongest
    distinct peaks of every B_k(theta) (at least one Rayleigh width apart in
    sin(theta)), casts one ray per peak, and takes as its two hypotheses the
    most-voted cell and the most-voted cell outside a cfg.mh_exclusion-cell
    neighborhood of it.  The second hypothesis is coherent if it collects at
    least cfg.mh_coherent times the votes of the first; otherwise every rule
    keeps the first.  Each AP votes at most once per cell, and cells within
    cfg.mh_ap_guard cells of a transmit AP are not taken as hypotheses.

    Returns a namespace of booleans (chosen cell is the target's cell):
    det_votes  more votes wins (ties at random)
    det_rand   one of two coherent hypotheses at random
    det_aware  of two coherent hypotheses, the one with weaker supporting
               peaks wins (the adversary knows the fake peak is the stronger)
    det_set    the target cell is one of the two hypotheses (two guesses)
    """
    cs, A, nx, ny = _grid(cfg)
    Ntx, Mtx = cfg.Ntx, cfg.Mtx
    theta, B = _scan(cfg, R)
    votes = np.zeros(nx * ny)
    rays = []
    for k in range(Ntx):
        ap_cells = []
        for p in _two_peaks(B[k], theta, Mtx):
            cand, cells, pick = _mirror_rays(cfg, k, theta[p], A, cs, nx, ny, rng)
            for c in pick:
                ap_cells.append(cells[c])
                rays.append((cells[c], B[k, p]))
        if ap_cells:
            votes[np.unique(np.concatenate(ap_cells))] += 1

    ix_all = np.arange(nx * ny) % nx
    iy_all = np.arange(nx * ny) // nx
    v1 = votes.copy()
    for k in range(Ntx):
        ia = np.floor((cfg.PosTx[k, 0] - A[0]) / cs)
        ja = np.floor((cfg.PosTx[k, 1] - A[2]) / cs)
        v1[np.maximum(np.abs(ix_all - ia), np.abs(iy_all - ja)) <= cfg.mh_ap_guard] = -np.inf
    h1 = _argmax_random(v1, rng)
    i1, j1 = h1 % nx, h1 // nx
    far = np.maximum(np.abs(ix_all - i1), np.abs(iy_all - j1)) > cfg.mh_exclusion
    v2 = v1.copy()
    v2[~far] = -np.inf
    h2 = _argmax_random(v2, rng)
    hyp = [h1, h2]

    # Mean peak power of the rays supporting each hypothesis (NaN if none).
    score = np.full(2, np.nan)
    for h in range(2):
        pw = [pw for cells, pw in rays if np.any(cells == hyp[h])]
        if pw:
            score[h] = np.mean(pw)

    tx, ty = _target_cell(cfg, A, cs, nx, ny)
    target = tx + ty * nx
    if votes[h1] > votes[h2]:
        pick_votes = h1
    else:
        pick_votes = hyp[rng.integers(2)]
    coherent = votes[h2] >= cfg.mh_coherent * votes[h1]
    pick_rand = h1
    pick_aware = h1
    if coherent:
        pick_rand = hyp[rng.integers(2)]
        if not (score[0] == score[1]):
            # MATLAB min ignores NaN and returns the first index if both are NaN.
            w = 0 if np.all(np.isnan(score)) else int(np.nanargmin(score))
            pick_aware = hyp[w]
        else:
            pick_aware = hyp[rng.integers(2)]
    return SimpleNamespace(det_votes=pick_votes == target, det_rand=pick_rand == target,
                           det_aware=pick_aware == target, det_set=target in hyp)
