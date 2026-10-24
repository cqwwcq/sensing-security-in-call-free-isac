# Python port of the simulation code

This folder is a Python rewrite of the MATLAB code in `common/`, `fig.7.detection/` and
`fig.9.outage/`, so that the revision simulations can run without MATLAB or CVX. Every
MATLAB function has a Python counterpart with the same algorithm, settings and defaults.

## Setup

```bash
cd python
pip install -r requirements.txt     # numpy, scipy, cvxpy (with Clarabel), matplotlib
```

Python 3.10 or later. CVX is replaced by CVXPY with the Clarabel interior-point solver,
which `pip install cvxpy` installs as well.

## Running the revision simulations

These are the steps of `RUN_SIMULATIONS.md`, in Python. Every driver takes
`--set key=value ...` to override any setting of `isac_sim/config.py` or of the driver
itself, `--out` for the output file, and `--workers` for the number of parallel processes
(default: all CPU cores). Results go to `results/` as `.npz` files.

Each finished channel realization is saved to a checkpoint folder next to the output file
(`results/<name>.parts_<hash>/`). If a sweep is interrupted, start the same command again
and it resumes where it stopped. The hash changes with the settings, so different settings
never mix.

```bash
# Step 1: smoke test (about a minute)
python run_detection_sweep.py --out results/smoke_detection.npz --set num_runs=1 "rho_list=[0.75, 1.05]"

# Step 2: check against the published results (published design and virtual location)
python run_detection_sweep.py --out results/validate_published.npz --set num_runs=100 \
    "PosFake=[5, -25]" sym_null_prop=false fake_tleak=false fake_resolvable_only=false
python make_detection_results.py results/validate_published.npz results/validate_published_results.npz

# Step 3: detection sweep, Figs. 2-3 and the letter figures
python run_detection_sweep.py                    # 500 realizations
python make_detection_results.py                 # printout of the quoted numbers
python make_plots.py fig2 fig3 pd mh scatter     # add --paper to replace the paper's PDFs
# power-matched operating point for Reviewer 2, Comment 5
python run_detection_sweep.py --out results/detection_sweep_matched.npz \
    --set num_runs=200 "rho_list=[0.9, 1.0, 1.05, 1.1]" fake_upper=1.2
python make_detection_results.py results/detection_sweep_matched.npz results/detection_results_matched.npz

# Step 4: outage sweep, Fig. 4 and Reviewer 1, Comment 5
python run_outage_sweep.py
python make_plots.py fig4

# Step 5: synchronization sweep, Reviewer 2, Comment 4
python run_sync_sweep.py

# Step 6: communication tradeoff, Reviewer 3, Comment 3
python run_comm_tradeoff.py
python run_comm_tradeoff.py --out results/comm_tradeoff_rho030.npz --set rho_ref=0.3
python run_comm_tradeoff.py --out results/comm_tradeoff_rho050.npz --set rho_ref=0.5
```

Load a result in Python with `isac_sim.runner.load_npz(path)`, which returns a dict of
arrays plus the saved configuration under `'cfg'`. Realizations are on the last axis, as
in the MATLAB `.mat` files.

## MATLAB to Python map

| MATLAB | Python |
|---|---|
| `common/sim_config.m`, `apply_overrides.m` | `isac_sim/config.py` |
| `common/compute_angle_dist.m`, `beamsteering.m`, `setup_geometry.m` | `isac_sim/geometry.py` |
| `common/gen_channels.m` | `isac_sim/channels.py` |
| `AN_update_receiver_mmse_null.m` (both folders) | `isac_sim/receiver.py` |
| `AN_run_optimization_mi.m`, `AN_run_optimization_robust.m` | `isac_sim/solver.py` (`run_optimization`; robust when `beta_error` is given) |
| `fig.7.detection/AN_Em.m` | `isac_sim/attack.py` |
| `common/ray_cells.m`, `localize_target.m`, `localize_multi.m` | `isac_sim/localization.py` |
| `common/user_sinr.m`, `sensing_metrics.m`, `detection_pd.m`, `draw_nullonly_an.m`, `wilson_ci.m` | `isac_sim/metrics.py` |
| `run_detection_sweep.m`, `make_detection_results.m` | `isac_sim/detection.py`, `run_detection_sweep.py`, `make_detection_results.py` |
| `run_outage_sweep.m` | `isac_sim/outage.py`, `run_outage_sweep.py` |
| `run_sync_sweep.m` | `isac_sim/sync.py`, `run_sync_sweep.py` |
| `run_comm_tradeoff.m` | `isac_sim/tradeoff.py`, `run_comm_tradeoff.py` |
| `paper_fig_style.m`, `plot_*.m` | `isac_sim/plots.py`, `make_plots.py` |

Array layouts: channels are `Nue x Ntx x Mue x Mtx` (`H[i, k]` is the channel from AP `k`
to user `i`), beams `w_i` are `Ntx x Mtx x Nue`, `w_t` and `z` are `Ntx x Mtx`, combiners
`u` are `Nue x Mue`. Indices are 0-based in code; `cfg.adversary` keeps the 1-based value
of the MATLAB settings.

## Differences from the MATLAB code

- **Random numbers.** NumPy's generator replaces MATLAB's, seeded with the same
  `base_seed + r` per realization, so individual realizations differ from the MATLAB
  runs. Averages agree within Monte Carlo error.
- **Solver.** CVXPY with Clarabel replaces CVX with SDPT3 (`cvx_precision high`). The
  transmit subproblem is the same convex problem; it is solved in normalized units
  (beams divided by `sqrt(P_T)`, SINR and null constraints divided by a channel-power
  scale), an exact change of variables that avoids the badly scaled raw units (noise
  about 4e-13 W, link SNRs about 70 dB). Each CVXPY problem is compiled once per process
  and design mode and then re-solved with new parameters.
- **CCCP stopping rule.** Reproduced exactly: `|f - f_prev| / (1e-9 + |f_prev|) < tol`
  with `f` in physical units. Because `f` is about `1e-12`, the `1e-9` term dominates and
  the rule acts as an absolute tolerance, so the loop usually stops after one or two
  iterations, in the MATLAB code and the published runs as well. Pass `tol_floor=0` in
  `solver_opts` for a true relative tolerance, e.g.
  `--set "solver_opts={'max_iter': 30, 'tol': 4e-3, 'tol_floor': 0}"`.
- **EM attack.** Vectorized over the snapshots, with each snapshot's own stopping test as
  in the MATLAB loops. The Gamma variate that `AN_Em.m` draws in every iteration is never
  used and is omitted. The least-squares step uses `(Eh^H Eh + Mue*Omega)^{-1} Eh^H y`,
  which equals the MATLAB form with `Htil = [Eh; chol(Mue*Omega)]`.
- **MI-minimizing baseline.** As in MATLAB, it needs a solver that is not in the
  repository (`AN_run_optimization_zou.m`). Add `isac_sim/zou.py` with
  `run_optimization_zou(cfg, geo, Hhat, rho, opts)` returning the namespace of
  `solver.run_optimization` to include it; otherwise it is skipped.
- **Outputs.** `.npz` instead of `.mat`, with the same variable names, plus `emcap_*`,
  the number of snapshots in which the EM loop reached its 200-iteration cap. Figures
  default to `results/figures/`; `make_plots.py --paper` writes Figs. 2-4 into the paper
  folder, as the MATLAB plot scripts do.

The settings marked `[ASSUMED]` in `isac_sim/config.py` are the same guesses as in
`common/sim_config.m`. Replace them with the values of the original driver scripts if
they are known.
