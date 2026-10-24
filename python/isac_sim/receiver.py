"""MMSE combiner of eq:mmse under a selectable null set.

Python port of AN_update_receiver_mmse_null.m (both versions; the robust
version's per-user noise level is accepted through sigma_n).
"""
import numpy as np
from scipy.linalg import null_space


def stack_beams(w_i, w_t, z):
    """Ntx x Mtx x (Nue+2) array [w_1,k ... w_Nue,k, w_t,k, z_k] per AP."""
    return np.concatenate([w_i, w_t[:, :, None], z[:, :, None]], axis=2)


def stream_responses(H, w_i, w_t, z):
    """q[i, :, l] = sum_k H[i, k] x_{l,k}: Nue x Mue x (Nue+2)."""
    return np.einsum('ikmn,knl->iml', H, stack_beams(w_i, w_t, z))


def update_receiver_mmse_null(H, w_i, w_t, z, sigma_n, rx_mode=2):
    """Unit-norm combiners u (Nue x Mue) maximizing the post-combining SINR.

    The combiner of user i maximizes |u^H q_ii|^2 / (u^H R_i u) over
    {u : Q_i^H u = 0, ||u|| = 1}, where Q_i matches the null constraints
    that the transmit design imposes:
        rx_mode = 0  ->  Q_i = []              (no nulls)
        rx_mode = 1  ->  Q_i = q_iz            (AN-only null, prob:null)
        rx_mode = 2  ->  Q_i = [q_iz, q_it]    (prob:null and prob:null-t)
    sigma_n is a scalar or one effective noise level per user (robust design,
    sigma_i^2 + L_i^csi).  As in the MATLAB code, R_i contains the other
    users' streams, the sensing stream and the noise, but not the AN.
    """
    Nue = H.shape[0]
    Mue = H.shape[2]
    sigma_vec = np.broadcast_to(np.asarray(sigma_n, dtype=float), (Nue,))
    q = stream_responses(H, w_i, w_t, z)
    u = np.zeros((Nue, Mue), dtype=complex)
    for i in range(Nue):
        R = sigma_vec[i] * np.eye(Mue, dtype=complex)
        for j in range(Nue):
            if j != i:
                R += np.outer(q[i, :, j], q[i, :, j].conj())
        q_it = q[i, :, Nue]
        q_iz = q[i, :, Nue + 1]
        q_ii = q[i, :, i]
        R += np.outer(q_it, q_it.conj())

        if rx_mode == 0:
            Q = np.zeros((Mue, 0), dtype=complex)
        elif rx_mode == 1:
            Q = q_iz[:, None]
        elif rx_mode == 2:
            Q = np.stack([q_iz, q_it], axis=1)
        else:
            raise ValueError(f'Unknown rx_mode {rx_mode}.')

        Rinv_q = np.linalg.solve(R, q_ii)
        if Q.shape[1] == 0:
            u_t = Rinv_q
        else:
            Rinv_Q = np.linalg.solve(R, Q)
            gram = Q.conj().T @ Rinv_Q
            u_t = Rinv_q - Rinv_Q @ (np.linalg.pinv(gram) @ (Q.conj().T @ Rinv_q))
            if np.linalg.norm(u_t) <= 1e-10 * max(np.linalg.norm(Rinv_q), 1.0):
                # Numerically stable fallback: solve the MMSE problem directly
                # in the Euclidean null space of Q_i^H.
                N = null_space(Q.conj().T)
                if N.shape[1] == 0:
                    raise RuntimeError(
                        f'The null-constrained combiner space is empty for user {i + 1}.')
                Rr = N.conj().T @ R @ N
                u_t = N @ np.linalg.solve(Rr, N.conj().T @ q_ii)
        if np.linalg.norm(u_t) <= 0:
            u_t = np.ones(Mue, dtype=complex) / np.sqrt(Mue)
        u[i] = u_t / np.linalg.norm(u_t)
    return u
