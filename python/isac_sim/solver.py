"""Alternating CCCP-MMSE beam design of Section III.

Python port of fig.7.detection/AN_run_optimization_mi.m and
fig.9.outage/AN_run_optimization_robust.m (the robust solver differs only by
the optional CSI-error floor, passed here as beta_error).  CVX is replaced by
CVXPY with the Clarabel interior-point solver.

Numerical scaling.  The physical problem mixes powers of about 1 W with noise
powers of about 4e-13 W and two-way path gains of about 1e-13, which is fine
for CVX/SDPT3 but poorly scaled for most conic solvers.  The convex subproblem
is therefore solved in normalized units, which is an exact change of
variables: beams are divided by sqrt(P_T), the SINR and null constraints are
divided by a channel-power scale nu (so that the effective channels have
order-one entries; the noise then enters as the small constant
sigma_ue^2 / nu), and the objective weights c_k are multiplied by
P_T / sigma_rx^2, so that the objective is in units of the sensing noise.

Each CVXPY problem is built once per process and design mode and then
re-solved with new parameter values (DPP), so that the per-iteration cost is
the solver time only.
"""
import time
from types import SimpleNamespace

import cvxpy as cp
import numpy as np

from .config import db2lin, dbm2w
from .metrics import effective_channel_matrix, sensing_c
from .receiver import update_receiver_mmse_null

_SOLVED = (cp.OPTIMAL, cp.OPTIMAL_INACCURATE)
_CACHE = {}


def _build(Ntx, Mtx, Nue, use_an, use_fake, sym_null, tleak, mask, upper, robust, a_true, a_fake):
    """CVXPY problem of one CCCP transmit update, with parameters."""
    n = Ntx * Mtx
    ncol = Nue + 1 + (1 if use_an else 0)
    X = cp.Variable((n, ncol), complex=True)
    p = SimpleNamespace(X=X)

    # Block-row steering matrices: row k is a(phi_k)^H on the columns of AP k.
    Abd = np.zeros((Ntx, n), dtype=complex)
    Afd = np.zeros((Ntx, n), dtype=complex)
    Ssum = np.zeros((Ntx, n))
    for k in range(Ntx):
        Abd[k, k * Mtx:(k + 1) * Mtx] = a_true[k].conj()
        Afd[k, k * Mtx:(k + 1) * Mtx] = a_fake[k].conj()
        Ssum[k, k * Mtx:(k + 1) * Mtx] = 1.0

    p.G = cp.Parameter((Nue, n), complex=True)          # normalized u_i^H Hhat_i
    p.Hs = cp.Parameter((Nue, n), complex=True)         # conj(S_r,i) g_i / Gamma
    p.c0 = cp.Parameter(Nue, nonneg=True)               # |S_r,i|^2 / Gamma
    p.Cobj = cp.Parameter((Ntx, Nue + 1), complex=True)  # 2 c_k conj(a^H W_r,k)
    p.budget = cp.Parameter(nonneg=True)
    p.eps = cp.Parameter(nonneg=True)

    GX = p.G @ X
    absq = cp.square(cp.abs(GX))
    mask_int = np.ones((Nue, ncol))
    mask_int[np.arange(Nue), np.arange(Nue)] = 0.0
    interf = cp.sum(cp.multiply(mask_int, absq), axis=1)
    sig_lin = 2 * cp.real(cp.diag(p.Hs @ X[:, :Nue]))
    rowpow = cp.sum(cp.square(cp.abs(X)), axis=1)
    Pk = Ssum @ rowpow

    p.noise = cp.Parameter(nonneg=True)                 # sigma_ue^2 / nu
    lhs = interf + p.noise - sig_lin + p.c0
    if robust:
        p.Bt = cp.Parameter((Nue, Ntx), nonneg=True)     # normalized beta_err^T
        lhs = lhs + p.Bt @ Pk
    cons = [lhs <= 0, Pk <= p.budget]
    if use_an:
        cons.append(cp.abs(GX[:, Nue + 1]) <= p.eps)
    if sym_null:
        cons.append(cp.abs(GX[:, Nue]) <= p.eps)

    AX = Abd @ X[:, :Nue + 1]                           # a_k^H W_k, Ntx x (Nue+1)
    obj = cp.real(cp.sum(cp.multiply(p.Cobj, AX)))

    idx = np.flatnonzero(mask)
    if use_fake and idx.size > 0:
        m = idx.size
        xz = X[:, Nue + 1]
        AXm = AX[idx, :]
        Fz = Afd[idx] @ xz
        Tz = Abd[idx] @ xz
        p.rho = cp.Parameter(nonneg=True)
        p.fr = cp.Parameter(m, complex=True)            # conj(a_f^H z_r)
        p.fr2 = cp.Parameter(m, nonneg=True)            # |a_f^H z_r|^2
        lhs_f = p.rho * cp.sum(cp.square(cp.abs(AXm)), axis=1) \
            - (2 * cp.real(cp.multiply(p.fr, Fz)) - p.fr2)
        if tleak:
            lhs_f = lhs_f + cp.square(cp.abs(Tz))
        cons.append(lhs_f <= 0)
        if upper:
            # Affine lower bounds of ||a^H W_k||^2 and T_k at the current
            # iterate keep the cap convex and conservative.
            p.U1 = cp.Parameter((m, Nue + 1), complex=True)   # kappa rho conj(aW_r)
            p.U1c = cp.Parameter(m, nonneg=True)              # kappa rho |aW_r|^2
            rhs = 2 * cp.real(cp.sum(cp.multiply(p.U1, AXm), axis=1)) - p.U1c
            if tleak:
                p.U2 = cp.Parameter(m, complex=True)          # kappa conj(az_r)
                p.U2c = cp.Parameter(m, nonneg=True)          # kappa |az_r|^2
                rhs = rhs + 2 * cp.real(cp.multiply(p.U2, Tz)) - p.U2c
            cons.append(cp.square(cp.abs(Fz)) <= rhs)
    p.idx = idx
    p.prob = cp.Problem(cp.Maximize(obj), cons)
    return p


def _get_problem(cfg_key, **kw):
    if cfg_key not in _CACHE:
        _CACHE[cfg_key] = _build(**kw)
    return _CACHE[cfg_key]


def run_optimization(cfg, geo, Hhat, P_T_dBm, Gamma_comm_dB, rho, opts=None, verbose=False):
    """Beam design of problem (9), or of the robust problem (20).

    Hhat: Nue x Ntx x Mue x Mtx channel estimates.  opts (dict) may contain
    max_iter, tol, tol_floor, sym_null, fake_tleak, fake_mask, fake_upper, an_mode
    ('fake', 'none', 'nullonly') and beta_error (Ntx x Nue, robust design).

    Returns a namespace with gamma_t (linear sensing SNR (2)), w_i (Ntx x Mtx
    x Nue), w_t, z (Ntx x Mtx), u (Nue x Mue), success, iterations,
    obj_history and solve_time, plus the constraint checks of the MATLAB
    solver's final report (min_sinr_dB, max_leakage, min_bp_difference).
    """
    t0 = time.time()
    opts = dict(opts or {})
    Ntx, Nue, Mtx, Mue = cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Mue
    P_T_lin = dbm2w(P_T_dBm)
    Gamma = db2lin(Gamma_comm_dB)
    sigma_u = dbm2w(cfg.sigma_ue_dBm)
    sigma_rx = dbm2w(cfg.sigma_rx_dBm)

    max_iter = opts.get('max_iter', 30)
    tol = opts.get('tol', 4e-3)
    # Stopping rule of the MATLAB code: |f - f_prev| / (tol_floor + |f_prev|)
    # < tol, with f in physical units (about 1e-12).  With the MATLAB value
    # tol_floor = 1e-9 the floor dominates, so the rule acts as an absolute
    # tolerance of about 4e-12 and the loop typically stops after 1-2
    # iterations, as in the published runs.  tol_floor = 0 gives a true
    # relative tolerance.
    tol_floor = opts.get('tol_floor', 1e-9)
    sym_null = bool(opts.get('sym_null', False))
    tleak = bool(opts.get('fake_tleak', False))
    mask = np.asarray(opts.get('fake_mask', np.ones(Ntx, bool)), dtype=bool)
    upper = opts.get('fake_upper', np.inf)
    an_mode = opts.get('an_mode', 'fake')
    beta_err = opts.get('beta_error', None)
    robust = beta_err is not None and np.size(beta_err) > 0

    if an_mode == 'fake':
        use_an, use_fake = True, True
    elif an_mode == 'none':
        use_an, use_fake, sym_null = False, False, False
    elif an_mode == 'nullonly':
        # Null-space-only (Gaussian-AN) baseline: prob:null and prob:null-t
        # are enforced but prob:fake is omitted.  The transmit beams are
        # optimized over the reduced budget left by the prescribed AN power,
        # and the AN itself is drawn afterwards by projecting a random vector
        # onto the null space of the returned combiners.
        use_an, use_fake, sym_null = False, False, True
    else:
        raise ValueError(f'Unknown an_mode "{an_mode}".')
    if an_mode == 'nullonly':
        rx_mode = 2
    elif not use_an:
        rx_mode = 0
    elif sym_null:
        rx_mode = 2
    else:
        rx_mode = 1

    c = sensing_c(cfg, geo)
    a_true, a_fake = geo.a_true, geo.a_fake

    # --- Initialization along the dominant user-channel, true and fake directions
    p_wi = P_T_lin * 0.6 / Nue * 0.8
    p_kick = P_T_lin * 0.3 * 0.8
    P_an = P_T_lin * 0.1 * 0.8
    if not use_an:
        # The no-AN baseline starts from the same beams with the AN budget
        # reassigned to the sensing beam.
        p_kick += P_an
        P_an = 0.0
    w_i = np.zeros((Ntx, Mtx, Nue), dtype=complex)
    w_t = np.zeros((Ntx, Mtx), dtype=complex)
    z = np.zeros((Ntx, Mtx), dtype=complex)
    for k in range(Ntx):
        for i in range(Nue):
            _, _, Vh = np.linalg.svd(Hhat[i, k])
            w_i[k, :, i] = Vh[0].conj() * np.sqrt(p_wi)
        w_t[k] = a_true[k] / np.linalg.norm(a_true[k]) * np.sqrt(p_kick)
        z[k] = a_fake[k] / np.linalg.norm(a_fake[k]) * np.sqrt(P_an)

    def objective(w_i_, w_t_):
        W = np.concatenate([w_i_, w_t_[:, :, None]], axis=2)
        aW = np.einsum('km,kml->kl', a_true.conj(), W)
        return float(np.sum(c * np.sum(np.abs(aW) ** 2, axis=1)))

    key = (Ntx, Mtx, Nue, use_an, use_fake, sym_null, tleak, tuple(mask),
           bool(np.isfinite(upper)), robust, a_true.tobytes(), a_fake.tobytes())
    prob = _get_problem(key, Ntx=Ntx, Mtx=Mtx, Nue=Nue, use_an=use_an,
                        use_fake=use_fake, sym_null=sym_null, tleak=tleak,
                        mask=mask, upper=bool(np.isfinite(upper)), robust=robust,
                        a_true=a_true, a_fake=a_fake)

    scale = P_T_lin                                    # beams / sqrt(scale)
    # Channel normalization nu = P_T * max_i ||Hhat_i||_F^2 / Mue (W), the
    # power scale of the effective channels u_i^H Hhat_i x.  The SINR and null
    # constraints are expressed relative to nu instead of the noise, which is
    # about 70 dB weaker in this geometry.
    nu = scale * max(np.sum(np.abs(Hhat[i].transpose(1, 0, 2).reshape(Mue, -1)) ** 2)
                     / Mue for i in range(Nue))
    gscale = np.sqrt(scale / nu)
    prob.budget.value = P_T_lin / scale
    prob.noise.value = sigma_u / nu
    prob.eps.value = cfg.epsilon / np.sqrt(nu)
    if robust:
        prob.Bt.value = np.asarray(beta_err).T * scale / nu
    c_n = c * scale / sigma_rx

    last = objective(w_i, w_t)
    obj_history = [last]
    success = False
    u = None
    it = 0
    for it in range(1, max_iter + 1):
        if robust:
            P_ap = np.sum(np.abs(w_i) ** 2, axis=(1, 2)) + np.sum(np.abs(w_t) ** 2, axis=1) \
                + np.sum(np.abs(z) ** 2, axis=1)
            noise_eff = sigma_u + np.asarray(beta_err).T @ P_ap
        else:
            noise_eff = sigma_u
        u = update_receiver_mmse_null(Hhat, w_i, w_t, z, noise_eff, rx_mode)
        G = effective_channel_matrix(Hhat, u).reshape(Nue, Ntx * Mtx)   # physical units

        # Current iterate in normalized units.
        wi_n = w_i / np.sqrt(scale)
        wt_n = w_t / np.sqrt(scale)
        z_n = z / np.sqrt(scale)
        Gn = G * gscale
        S_r = np.array([Gn[i] @ wi_n[:, :, i].reshape(-1) for i in range(Nue)])
        prob.G.value = Gn
        prob.Hs.value = (S_r.conj()[:, None] * Gn) / Gamma
        prob.c0.value = np.abs(S_r) ** 2 / Gamma
        Wn = np.concatenate([wi_n, wt_n[:, :, None]], axis=2)
        aW_r = np.einsum('km,kml->kl', a_true.conj(), Wn)        # Ntx x (Nue+1)
        prob.Cobj.value = 2 * c_n[:, None] * aW_r.conj()
        if use_fake and prob.idx.size > 0:
            idx = prob.idx
            fr = np.einsum('km,km->k', a_fake[idx].conj(), z_n[idx])
            prob.rho.value = float(rho)
            prob.fr.value = fr.conj()
            prob.fr2.value = np.abs(fr) ** 2
            if np.isfinite(upper):
                kr = upper * rho
                prob.U1.value = kr * aW_r[idx].conj()
                prob.U1c.value = kr * np.sum(np.abs(aW_r[idx]) ** 2, axis=1)
                if tleak:
                    az = np.einsum('km,km->k', a_true[idx].conj(), z_n[idx])
                    prob.U2.value = upper * az.conj()
                    prob.U2c.value = upper * np.abs(az) ** 2

        try:
            prob.prob.solve(solver=cp.CLARABEL)
            status = prob.prob.status
        except cp.SolverError:
            status = 'solver_error'
        if status not in _SOLVED:
            if verbose:
                print(f'    [Stage 2] Iteration {it}: solver failed. Status: {status}')
            success = False
            break

        Xv = prob.X.value * np.sqrt(scale)
        Xv = Xv.reshape(Ntx, Mtx, -1)
        w_i_new = Xv[:, :, :Nue]
        w_t_new = Xv[:, :, Nue]
        z_new = Xv[:, :, Nue + 1] if use_an else np.zeros((Ntx, Mtx), dtype=complex)
        cur = objective(w_i_new, w_t_new)
        obj_history.append(cur)
        if verbose:
            print(f'    [Stage 2] Iteration {it}: sensing SNR = {10 * np.log10(cur / sigma_rx):.8f} dB')
        converged = abs(cur - last) / (tol_floor + abs(last)) < tol
        last = cur
        w_i, w_t, z = w_i_new, w_t_new, z_new
        success = True
        if converged:
            if verbose:
                print('    [Stage 2] Algorithm converged.')
            break
        if it == max_iter and verbose:
            print('    [Stage 2] Maximum number of iterations reached.')

    out = SimpleNamespace(gamma_t=last / sigma_rx, w_i=w_i, w_t=w_t, z=z, u=u,
                          success=success, iterations=it, obj_history=np.array(obj_history),
                          solve_time=time.time() - t0)
    _final_checks(cfg, geo, Hhat, out, sigma_u, verbose, Gamma_comm_dB)
    return out


def _final_checks(cfg, geo, Hhat, out, sigma_u, verbose, Gamma_comm_dB):
    """Constraint report of the MATLAB solver (SINR, AN leakage, C4)."""
    from .metrics import user_sinr
    if out.u is None:
        out.min_sinr_dB = np.nan
        out.max_leakage = np.nan
        out.min_bp_difference = np.nan
        return
    sinr = user_sinr(cfg, Hhat, out.u, out.w_i, out.w_t, out.z)
    out.min_sinr_dB = float(10 * np.log10(np.min(sinr)))
    G = effective_channel_matrix(Hhat, out.u)
    out.max_leakage = float(np.max(np.abs(np.einsum('ikn,kn->i', G, out.z))))
    Wtot = np.concatenate([out.w_i, out.w_t[:, :, None], out.z[:, :, None]], axis=2)
    bt = np.sum(np.abs(np.einsum('km,kml->kl', geo.a_true.conj(), Wtot)) ** 2, axis=1)
    bf = np.sum(np.abs(np.einsum('km,kml->kl', geo.a_fake.conj(), Wtot)) ** 2, axis=1)
    out.min_bp_difference = float(np.min(bf - bt))
    if verbose:
        print(f'    [Check] Final sensing SNR: {10 * np.log10(out.gamma_t):.7f} dB')
        print(f'    [Check] Worst-case (minimum) SINR: {out.min_sinr_dB:.7f} dB '
              f'(target {Gamma_comm_dB:.1f} dB)')
        print(f'    [Check] Worst-case (maximum) leakage: {out.max_leakage:.7e} (linear)')
        print(f'    [Check] Minimum beampattern difference (constraint C4): '
              f'{out.min_bp_difference:.7e}')
