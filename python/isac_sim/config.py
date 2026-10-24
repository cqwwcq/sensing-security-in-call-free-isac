"""Default simulation parameters shared by all sweep drivers.

Python port of common/sim_config.m and common/apply_overrides.m.

Values marked [Table I] are taken from Table I of paper_draft.tex.  Values
marked [ASSUMED] are not stated in the paper and were filled in from the
earlier multi-static ISAC paper in this repository (path-loss exponent, RCS,
UE placement, AP ordering, search area) or chosen here.  Replace them with
the values of the original driver scripts to reproduce the published numbers
exactly.
"""
import ast
import copy
from types import SimpleNamespace

import numpy as np


def sim_config():
    cfg = SimpleNamespace()

    # --- Network size and antennas [Table I] ---
    cfg.Ntx = 8
    cfg.Nrx = 4
    cfg.Nue = 4
    cfg.Mtx = 16
    cfg.Mue = 16
    cfg.Mrx = 1           # [ASSUMED] the solver's c_k has no receive array gain

    # --- Power, QoS, noise [Table I] ---
    cfg.P_T = 30.0                # dBm per transmit AP
    cfg.Gamma_comm = 10.0         # dB
    cfg.sigma_ue_dBm = -94.0      # user noise power
    cfg.sigma_rx_dBm = -94.0      # sensing-receiver noise power
    cfg.sigma_p_dBm = -94.0       # pilot noise power
    cfg.epsilon = 1e-10           # null-constraint tolerance used in the sweeps

    # --- Carrier and pilots [Table I] ---
    cfg.fc = 1.9e9
    cfg.lambda_c = 3e8 / cfg.fc
    cfg.tau_p = 64
    cfg.P_p_dBm = 20.0            # total pilot power per user

    # --- Geometry [Table I], ordering of the earlier paper [ASSUMED] ---
    cfg.PosTx = np.array([[-200, -200], [200, 200], [-200, 200], [200, -200],
                          [0, -200], [0, 200], [-200, 0], [200, 0]], dtype=float)
    cfg.PosRx = np.array([[150, 150], [-150, -150], [-150, 150], [150, -150]],
                         dtype=float)
    cfg.PosTarget = np.array([-75.0, 75.0])
    # Relocated so that all TX APs resolve the true and fake directions
    # (Reviewer 2, Comment 5); published: [5, -25], resolvable at 5 of 8 APs.
    cfg.PosFake = np.array([175.0, 35.0])
    cfg.PosUE = np.array([[120, 120], [-120, -120], [-120, 120], [120, -120]],
                         dtype=float)                       # [ASSUMED]
    cfg.adversary = 1             # 1-based, as in AN_Em.m: user 1 is the adversary

    # --- Large-scale fading and target [ASSUMED, earlier paper] ---
    cfg.pathloss_exp = 3.0        # beta_ik = d_ik^(-pathloss_exp)
    cfg.sigmarcs = np.sqrt(10.0)  # RCS variance 10 dBsm (solver uses sigmarcs^2)

    # --- Reconstruction attack and localization ---
    cfg.N_snapshot = 128          # [Table I]
    cfg.cell_size = 10.0          # m [Table I]
    cfg.cell_list = [10.0, 20.0, 50.0]   # coarser cells, more favorable adversary
    cfg.n_ang = 3604              # angle-scan points over [-pi, pi]
    cfg.area = [-200.0, 200.0, -200.0, 200.0]   # [xmin xmax ymin ymax] [ASSUMED]
    cfg.mirror_rule = 'random'    # ULA front-back ambiguity: 'random' or 'both'
    cfg.mh_exclusion = 3          # multi-hypothesis attack: second hypothesis lies
                                  # more than this many cells from the first
    cfg.mh_coherent = 0.6         # ... and competes only with at least this
                                  # fraction of the first hypothesis' votes
    cfg.mh_ap_guard = 2           # no hypotheses within this many cells of a TX AP

    # --- Design sweep ---
    cfg.rho_list = [0.05, 0.15, 0.3, 0.5, 0.75, 1.05, 1.5, 2.2, 3.5]
    cfg.sym_null_prop = True      # enforce (9f) in the proposed design, as the
                                  # paper states; the published sweep used False
    cfg.fake_tleak = True         # revised (9g) with T_k on the right-hand side
                                  # (Reviewer 1, Comment 2); published: False
    cfg.fake_resolvable_only = True   # impose (9g) only at APs with resolvable
                                  # true and fake directions; published: False
    cfg.fake_upper = np.inf       # cap kappa on the fake peak; e.g. 1.2 with
                                  # rho ~ 1 gives power-matched peaks
    cfg.solver_opts = {'max_iter': 30, 'tol': 4e-3}

    # --- Detection probability (Reviewer 3, Comment 2) ---
    cfg.Pfa = 1e-6
    cfg.pd_N_list = [5, 128]              # integration lengths (N >= Nue+1)
    cfg.pd_offset_dB = list(range(-30, 1, 2))   # SNR offsets for Pd-vs-RCS curves
    cfg.pd_alpha_draws = 2000             # Swerling-I echo draws per realization

    # --- Monte Carlo ---
    cfg.workers = 0               # parallel processes; 0 = all CPU cores, 1 = serial
    return cfg


def apply_overrides(cfg, overrides):
    """Replace attributes of cfg by the entries of the dict overrides.

    Lets a driver run with different settings without editing any file.
    Unknown keys are an error, as in apply_overrides.m.
    """
    cfg = copy.deepcopy(cfg)
    for key, val in (overrides or {}).items():
        if not hasattr(cfg, key):
            raise KeyError(f'apply_overrides: unknown setting "{key}".')
        old = getattr(cfg, key)
        if isinstance(old, np.ndarray):
            val = np.asarray(val, dtype=old.dtype)
        setattr(cfg, key, val)
    return cfg


def parse_overrides(items):
    """Turn command-line ['num_runs=1', 'PosFake=[5,-25]'] into a dict."""
    out = {}
    for item in items or []:
        if '=' not in item:
            raise ValueError(f'Expected key=value, got "{item}".')
        key, text = item.split('=', 1)
        text = text.strip()
        if text.lower() in ('inf', '+inf'):
            val = np.inf
        elif text.lower() in ('true', 'false'):
            val = text.lower() == 'true'
        else:
            try:
                val = ast.literal_eval(text)
            except (ValueError, SyntaxError):
                val = text
        out[key.strip()] = val
    return out


def cfg_to_dict(cfg):
    """Flat dict of the configuration, for saving next to the results."""
    out = {}
    for key, val in vars(cfg).items():
        if isinstance(val, dict):
            for k2, v2 in val.items():
                out[f'{key}.{k2}'] = v2
        else:
            out[key] = val
    return out


def db2lin(x_db):
    return 10.0 ** (np.asarray(x_db) / 10.0)


def dbm2w(x_dbm):
    return 10.0 ** (np.asarray(x_dbm) / 10.0) / 1000.0
