"""Adversary-side EM reconstruction of the per-AP transmit signals.

Python port of fig.7.detection/AN_Em.m, vectorized over the snapshots: every
snapshot runs its own EM iterations with its own stopping test, exactly as in
the MATLAB loops.

As in AN_Em.m, the adversary is genie-aided: it observes each AP's signal
H_{a,k} (W_{-a,k} s + z_k s_z) + noise separately, and its EM routine is
initialized with and anchored to the actual channel H_{a,k}.  The same
receiver-noise sample is added to every AP's observation of a snapshot.

The MATLAB code draws a Gamma variate in every EM iteration whose value is
never used; it is omitted here, since it affects only the random stream.
"""
import numpy as np


def an_em(Dsn, nvariance, H, noise, W, z, rng, adversary=0, max_em_iter=200):
    """Per-AP sample covariances R_{xhat_k} of (3) and the EM cap count.

    Dsn       Ns x N unit-power symbols; rows 0..Nue-1 data, row Nue sensing
    nvariance adversary noise variance sigma^2 (Covnd in AN_Em.m)
    H         Ntx x Mue x Mtx actual channels of the adversary
    noise     Mue x N receiver noise, one column per snapshot
    W         Ntx x Mtx x (Nue+1) beams [w_1,k ... w_Nue,k, w_t,k]
    z         Ntx x Mtx AN beams
    Returns R (Ntx x Mtx x Mtx, mean over the snapshots of xhat xhat^H) and
    cap_hits, the number of snapshots whose EM loop reached max_em_iter.
    """
    Ntx, Mue, Mtx = H.shape
    Nue = W.shape[2] - 1
    N = Dsn.shape[1]
    cn = Dsn[:Nue]
    sn = Dsn[Nue]
    others = [j for j in range(Nue) if j != adversary]
    cn_aux = cn[others]                                    # (Nue-1) x N
    jam = np.exp(1j * 2 * np.pi * rng.random(N))           # AN symbols, unit modulus
    Covnd = nvariance
    prior = 1e-2
    a_par = (10 + 2 * Mue) / 2
    eye = np.eye(Mtx)

    R = np.zeros((Ntx, Mtx, Mtx), dtype=complex)
    cap_hits = 0
    for k in range(Ntx):
        Hk = H[k]
        Waux = W[k][:, others]                              # Mtx x (Nue-1)
        x_other = Waux @ cn_aux                             # Mtx x N
        x_sense = np.outer(W[k][:, Nue], sn)                # Mtx x N
        Y = Hk @ (x_other + np.outer(z[k], jam) + x_sense) + noise      # Mue x N
        Y = Y.T                                             # N x Mue
        S = (np.linalg.pinv(Hk.conj().T @ Hk) @ Hk.conj().T @ Y.T).T    # N x Mtx

        # Initial convergence measure of AN_Em.m: spectral norm of the
        # matrix W_t s_t[n] + W_aux C_aux (broadcast over the snapshots).
        cvg = np.array([np.linalg.norm(x_other + x_sense[:, [n]], 2) for n in range(N)])

        Eh = Hk[None] + rng.normal(0, prior, (N, Mue, Mtx)) \
            + 1j * rng.normal(0, prior, (N, Mue, Mtx))
        Omega = np.linalg.inv(np.einsum('ni,nj->nij', S, S.conj()) + eye / prior)

        it = np.zeros(N, dtype=int)
        active = cvg >= 1e-5
        while np.any(active):
            act = np.flatnonzero(active)
            s = S[act]
            y = Y[act]
            E = Eh[act]
            Om = Omega[act]
            it[act] += 1
            r = y - np.einsum('nij,nj->ni', E, s)
            quad = np.einsum('ni,nij,nj->n', s.conj(), Om, s)
            Cl = (np.einsum('ni,ni->n', r.conj(), r) + quad * Mue) / Covnd
            b_par = (10 + 2 * Cl) / 2
            Eu = a_par / b_par
            Om = np.linalg.inv(eye / prior
                               + np.einsum('ni,nj->nij', s, s.conj()) * (Eu.real / Covnd)[:, None, None])
            rhs = Hk.conj().T[None] / prior \
                + np.einsum('ni,nj->nij', s, y.conj()) * (Eu / Covnd)[:, None, None]
            E = np.conj(np.einsum('nij,njk->nik', Om.conj(), rhs)).transpose(0, 2, 1)
            # sdhat = (Htil^H Htil)^{-1} Htil^H ytil with Htil = [Eh; chol(Omega Mue)]
            # and ytil = [y; 0], i.e. (Eh^H Eh + Mue Omega)^{-1} Eh^H y.
            A = np.einsum('nji,njk->nik', E.conj(), E) + Mue * Om
            s_new = np.linalg.solve(A, np.einsum('nji,nj->ni', E.conj(), y)[..., None])[..., 0]
            cvg_new = np.linalg.norm(s_new - s, axis=1)
            S[act] = s_new
            Eh[act] = E
            Omega[act] = Om
            cvg[act] = cvg_new
            active = (cvg >= 1e-5) & (it < max_em_iter)
        cap_hits += int(np.sum(it >= max_em_iter))
        R[k] = np.einsum('ni,nj->ij', S, S.conj()) / N
    return R, cap_hits
