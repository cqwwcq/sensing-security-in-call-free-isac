"""Sweep behind Figs. 2-3: localization success and sensing SNR versus eta.

Python port of fig.7.detection/run_detection_sweep.m and
fig.7.detection/make_detection_results.m.

In each channel realization, the no-AN baseline is solved once, and for every
rho in cfg.rho_list the proposed design ('fake'), the null-space-only
baseline ('nullonly', paired with the proposed design's realized eta) and, if
available, the MI-minimizing baseline are solved on the LMMSE estimates.  The
reconstruction attack and the localization then run on the actual channels.

Saved arrays (realizations on the last axis), besides the variables of the
published sweep (det_*, err_*, eta_*, snr_*):
  snrz_*     sensing SNR (dB) with the known AN echo exploited
  fdom_*     fraction of APs whose reconstructed fake-direction response
             exceeds the true-direction response (Reviewer 1, Comment 2)
  pkf_*, pkt_*  fraction of APs whose reconstructed dominant peak lies at the
             fake / true direction; fdom_*, pkf_* and pkt_* count only APs
             with resolvable directions
  mhv_*, mhr_*, mha_*, mhs_*  success of the multi-hypothesis attack with the
             votes / random / mechanism-aware rule and with two guesses
  fhit_*     the dominant-peak estimate lies in the virtual location's cell
  detc_*     dominant-peak success with the search cells of cfg.cell_list
  est_*      estimated location of the dominant-peak attack
  pd_*       detection probability, nOffset x nN (per rho for prop/null/zou)
  emcap_*    snapshots whose EM loop reached its iteration cap

The MI-minimizing baseline needs a module isac_sim/zou.py with a function
run_optimization_zou(cfg, geo, Hhat, rho, opts) returning the same namespace
as solver.run_optimization.  It is not in the repository (neither is the
MATLAB AN_run_optimization_zou.m); without it the baseline is skipped.
"""
import numpy as np

from .attack import an_em
from .channels import gen_channels
from .config import dbm2w
from .localization import localize_multi, localize_target
from .metrics import detection_pd, draw_nullonly_an, sensing_metrics, wilson_ci
from .solver import run_optimization

try:
    from .zou import run_optimization_zou   # noqa: F401
    HAS_ZOU = True
except ImportError:
    HAS_ZOU = False

DRIVER_DEFAULTS = {'num_runs': 500, 'base_seed': 20260823}
FIELDS = ['det', 'err', 'eta', 'snr', 'snrz', 'fdom', 'pkf', 'pkt',
          'mhv', 'mhr', 'mha', 'mhs', 'fhit', 'emcap']


def fake_mask(cfg, geo):
    """APs at which the fake-peak constraint (9g) is imposed."""
    return geo.resolvable.copy() if cfg.fake_resolvable_only else np.ones(cfg.Ntx, bool)


def solve_design(cfg, geo, Hhat, P_T_dBm, rho, opts, Gamma_dB=None):
    o = dict(cfg.solver_opts)
    o.update(opts)
    g = cfg.Gamma_comm if Gamma_dB is None else Gamma_dB
    return run_optimization(cfg, geo, Hhat, P_T_dBm, g, rho, o)


def prop_opts(cfg, geo):
    return {'an_mode': 'fake', 'sym_null': cfg.sym_null_prop, 'fake_tleak': cfg.fake_tleak,
            'fake_mask': fake_mask(cfg, geo), 'fake_upper': cfg.fake_upper}


def evaluate_design(cfg, geo, ch, Dsn, noise_adv, w_i, w_t, z, rng):
    """Legitimate metrics on the designed beams, attack on the actual channels."""
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)
    a = cfg.adversary - 1
    W = np.concatenate([w_i, w_t[:, :, None]], axis=2)          # Ntx x Mtx x (Nue+1)
    Rk, cap = an_em(Dsn, sigma_ue, ch.H[a], noise_adv, W, z, rng, adversary=a)
    success, err, est, dom = localize_target(cfg, geo, Rk, rng)
    mh = localize_multi(cfg, geo, Rk, rng)
    m = sensing_metrics(cfg, geo, w_i, w_t, z)
    res = geo.resolvable
    e = {'det': float(success), 'err': err, 'eta': m.eta,
         'snr': 10 * np.log10(m.gamma_t), 'snrz': 10 * np.log10(m.gamma_t_z),
         'fdom': float(np.mean(dom.fake_dominant[res])),
         'pkf': float(np.mean(dom.peak_class[res] == 1)),
         'pkt': float(np.mean(dom.peak_class[res] == 2)),
         'mhv': float(mh.det_votes), 'mhr': float(mh.det_rand),
         'mha': float(mh.det_aware), 'mhs': float(mh.det_set),
         'emcap': float(cap)}
    corner = np.array([cfg.area[0], cfg.area[2]])
    fcell = np.floor((cfg.PosFake - corner) / cfg.cell_size)
    e['fhit'] = float(np.array_equal(np.floor((est - corner) / cfg.cell_size), fcell))
    e['est'] = est
    # Dominant-peak attack with coarser search cells (Reviewer 3, Comment 4).
    e['detc'] = np.array([float(localize_target(cfg, geo, Rk, rng, cell_size=cs)[0])
                          for cs in cfg.cell_list])
    e['pd'] = detection_pd(cfg, m.V, cfg.pd_offset_dB, cfg.pd_N_list, rng)
    return e


def _blank(cfg):
    e = {f: np.nan for f in FIELDS}
    e['est'] = np.full(2, np.nan)
    e['detc'] = np.full(len(cfg.cell_list), np.nan)
    e['pd'] = np.full((len(cfg.pd_offset_dB), len(cfg.pd_N_list)), np.nan)
    return e


def one_run(cfg, geo, r):
    """One channel realization of the sweep (one_run of the MATLAB driver)."""
    rng = np.random.default_rng(cfg.base_seed + r)
    ch = gen_channels(cfg, geo, rng, 'pilot')
    Ns = cfg.Nue + 1
    N = cfg.N_snapshot
    P_T_lin = dbm2w(cfg.P_T)
    sigma_ue = dbm2w(cfg.sigma_ue_dBm)

    # Unit-power QPSK symbols (rows 0..Nue-1: data, row Nue: sensing) and
    # adversary receiver noise, independent across snapshots.
    Dsn = (np.sign(rng.standard_normal((Ns, N))) + 1j * np.sign(rng.standard_normal((Ns, N)))) / np.sqrt(2)
    noise_adv = np.sqrt(sigma_ue / 2) * (rng.standard_normal((cfg.Mue, N))
                                         + 1j * rng.standard_normal((cfg.Mue, N)))

    def evaluate(w_i, w_t, z):
        return evaluate_design(cfg, geo, ch, Dsn, noise_adv, w_i, w_t, z, rng)

    out = {'ok': True, 'ok_zou': HAS_ZOU}
    s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, 0.0, {'an_mode': 'none'})
    out['ok'] &= s.success
    out['noan'] = evaluate(s.w_i, s.w_t, np.zeros_like(s.z))

    prop, null, zou = [], [], []
    for rho in cfg.rho_list:
        # --- Proposed fake-peak design ---
        s = solve_design(cfg, geo, ch.Hhat, cfg.P_T, rho, prop_opts(cfg, geo))
        out['ok'] &= s.success
        prop.append(evaluate(s.w_i, s.w_t, s.z))

        # --- Null-space-only Gaussian AN, paired at the same eta ---
        eta = prop[-1]['eta']
        P_beam_dBm = 10 * np.log10(max(1 - eta, 1e-6) * P_T_lin * 1000)
        s = solve_design(cfg, geo, ch.Hhat, P_beam_dBm, rho, {'an_mode': 'nullonly'})
        out['ok'] &= s.success
        z = draw_nullonly_an(cfg, ch.Hhat, s.u, s.w_i, s.w_t,
                             eta * P_T_lin * np.ones(cfg.Ntx), rng)
        null.append(evaluate(s.w_i, s.w_t, z))

        # --- MI-minimizing baseline (adapted from Zou et al.) ---
        if HAS_ZOU:
            s = run_optimization_zou(cfg, geo, ch.Hhat, rho, dict(cfg.solver_opts))
            out['ok_zou'] &= s.success
            zou.append(evaluate(s.w_i, s.w_t, s.z))
        else:
            zou.append(_blank(cfg))
    out['prop'], out['null'], out['zou'] = prop, null, zou
    return out


def collect(results, cfg):
    """Stack the per-run results; realizations on the last axis."""
    ok = np.array([r['ok'] for r in results])
    ok_zou = ok & np.array([r['ok_zou'] for r in results])
    S = {}
    for name, keep in (('noan', ok), ('null', ok), ('prop', ok), ('zou', ok_zou)):
        if name == 'zou' and not HAS_ZOU:
            continue
        rs = [r for r, k in zip(results, keep) if k]
        if name == 'noan':
            for f in FIELDS:
                S[f'{f}_{name}'] = np.array([r[name][f] for r in rs])
            S[f'est_{name}'] = np.stack([r[name]['est'] for r in rs], axis=-1)
            S[f'detc_{name}'] = np.stack([r[name]['detc'] for r in rs], axis=-1)
            S[f'pd_{name}'] = np.stack([r[name]['pd'] for r in rs], axis=-1)
            for f in ('eta', 'snrz', 'fdom'):
                S.pop(f'{f}_{name}')
        else:
            for f in FIELDS:
                S[f'{f}_{name}'] = np.array([[e[f] for e in r[name]] for r in rs]).T
            S[f'est_{name}'] = np.stack([np.stack([e['est'] for e in r[name]]) for r in rs], axis=-1)
            S[f'detc_{name}'] = np.stack([np.stack([e['detc'] for e in r[name]]) for r in rs], axis=-1)
            S[f'pd_{name}'] = np.stack([np.stack([e['pd'] for e in r[name]]) for r in rs], axis=-1)
    S['num_ok'] = int(ok.sum())
    S['num_ok_zou'] = int(ok_zou.sum())
    return S


def make_results(S, cfg):
    """Aggregate a sweep into the plotted values (make_detection_results.m)."""
    has_zou = 'det_zou' in S
    O = {}
    n = S['det_noan'].shape[-1]
    O['p_noan'] = np.mean(S['det_noan'])
    O['lo_noan'], O['hi_noan'] = wilson_ci(O['p_noan'], n)
    for d in ('null', 'prop'):
        O[f'p_{d}'] = np.mean(S[f'det_{d}'], axis=-1)
        O[f'lo_{d}'], O[f'hi_{d}'] = wilson_ci(O[f'p_{d}'], n)
    O['eta'] = np.mean(S['eta_prop'], axis=-1)
    O['eta_null'] = np.mean(S['eta_null'], axis=-1)
    O['snr_noan_mean'] = np.mean(S['snr_noan'])
    O['snr_null_mean'] = np.mean(S['snr_null'], axis=-1)
    O['snr_prop_mean'] = np.mean(S['snr_prop'], axis=-1)
    O['snrz_prop_mean'] = np.mean(S['snrz_prop'], axis=-1)
    O['fdom_prop_mean'] = np.mean(S['fdom_prop'], axis=-1)
    O['pkf_prop_mean'] = np.mean(S['pkf_prop'], axis=-1)
    O['pkt_prop_mean'] = np.mean(S['pkt_prop'], axis=-1)
    O['pkt_noan_mean'] = np.mean(S['pkt_noan'])
    O['p_detc_noan'] = np.mean(S['detc_noan'], axis=-1)        # nCell
    O['p_detc_prop'] = np.mean(S['detc_prop'], axis=-1)        # nRho x nCell
    O['p_detc_null'] = np.mean(S['detc_null'], axis=-1)
    for rule in ('mhv', 'mhr', 'mha', 'mhs', 'fhit'):
        O[f'p_{rule}_noan'] = np.mean(S[f'{rule}_noan'])
        O[f'p_{rule}_null'] = np.mean(S[f'{rule}_null'], axis=-1)
        O[f'p_{rule}_prop'] = np.mean(S[f'{rule}_prop'], axis=-1)
    O['pd_noan_mean'] = np.mean(S['pd_noan'], axis=-1)          # nOffset x nN
    O['pd_null_mean'] = np.mean(S['pd_null'], axis=-1)          # nRho x nOffset x nN
    O['pd_prop_mean'] = np.mean(S['pd_prop'], axis=-1)
    O['emcap_total'] = float(np.sum(S['emcap_noan']) + np.sum(S['emcap_prop'])
                             + np.sum(S['emcap_null']))
    if has_zou:
        nz = S['det_zou'].shape[-1]
        O['p_zou'] = np.mean(S['det_zou'], axis=-1)
        O['lo_zou'], O['hi_zou'] = wilson_ci(O['p_zou'], nz)
        O['eta_zou'] = np.mean(S['eta_zou'], axis=-1)
        O['snr_zou_mean'] = np.mean(S['snr_zou'], axis=-1)
        O['pd_zou_mean'] = np.mean(S['pd_zou'], axis=-1)
    return O


def print_results(O, cfg):
    """Printout of the values quoted in the paper and the response letter."""
    rho_list = np.asarray(cfg['rho_list'])
    j = int(np.argmin(np.abs(O['eta'] - 0.5)))
    print(f"No AN: success {100 * O['p_noan']:.2f} %, sensing SNR {O['snr_noan_mean']:.2f} dB")
    print(f"eta = {O['eta'][j]:.3f} (rho = {rho_list[j]:.2f}):")
    print(f"  proposed : success {100 * O['p_prop'][j]:.2f} %, SNR {O['snr_prop_mean'][j]:.2f} dB, "
          f"SNR with AN echo {O['snrz_prop_mean'][j]:.2f} dB, fake peak dominant at "
          f"{100 * O['fdom_prop_mean'][j]:.0f} % of APs")
    pkf, pkt = O['pkf_prop_mean'][j], O['pkt_prop_mean'][j]
    print(f"  proposed : dominant peak at fake / true / elsewhere: {100 * pkf:.0f} / "
          f"{100 * pkt:.0f} / {100 * max(0, 1 - pkf - pkt):.0f} % of APs "
          f"(no AN: {100 * O['pkt_noan_mean']:.0f} % at true)")
    for c, cs in enumerate(cfg['cell_list']):
        print(f"  {cs:2.0f} m cells: success no AN {100 * O['p_detc_noan'][c]:.1f} %, "
              f"proposed {100 * O['p_detc_prop'][j, c]:.1f} %, "
              f"null-only {100 * O['p_detc_null'][j, c]:.1f} %")
    print(f"  estimates in the virtual location's cell: no AN {100 * O['p_fhit_noan']:.1f} %, "
          f"proposed {100 * O['p_fhit_prop'][j]:.1f} %")
    print('  multi-hypothesis attack, success with votes / random / aware rule / two guesses:')
    print(f"    no AN   : {100 * O['p_mhv_noan']:5.1f} / {100 * O['p_mhr_noan']:5.1f} / "
          f"{100 * O['p_mha_noan']:5.1f} / {100 * O['p_mhs_noan']:5.1f} %")
    j1 = int(np.argmin(np.abs(rho_list - 1)))
    for jj in sorted({j, j1}):
        print(f"    proposed (rho = {rho_list[jj]:.2f}, eta = {O['eta'][jj]:.3f}): "
              f"{100 * O['p_mhv_prop'][jj]:5.1f} / {100 * O['p_mhr_prop'][jj]:5.1f} / "
              f"{100 * O['p_mha_prop'][jj]:5.1f} / {100 * O['p_mhs_prop'][jj]:5.1f} %")
    print(f"  null-only: success {100 * O['p_null'][j]:.2f} %, SNR {O['snr_null_mean'][j]:.2f} dB "
          f"(its own eta {O['eta_null'][j]:.3f})")
    if 'p_zou' in O:
        jz = int(np.argmin(np.abs(O['eta_zou'] - 0.5)))
        print(f"  MI-min   : success {100 * O['p_zou'][jz]:.2f} % at eta {O['eta_zou'][jz]:.3f}, "
              f"SNR {O['snr_zou_mean'][jz]:.2f} dB")
    offs = list(cfg['pd_offset_dB'])
    if 0 in offs:
        i0 = offs.index(0)
        for b, N in enumerate(cfg['pd_N_list']):
            print(f"  Pd (Pfa = {cfg['Pfa']:g}, N = {N}): no AN {O['pd_noan_mean'][i0, b]:.4f}, "
                  f"proposed {O['pd_prop_mean'][j, i0, b]:.4f}")
    print(f"  EM iteration cap reached in {O['emcap_total']:.0f} snapshots in total")
