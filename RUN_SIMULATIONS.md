# Running the revision simulations in MATLAB

These steps produce every number and figure that the reply letter and the revised paper
need. Run them in order. Steps 1 and 2 are short checks; do not skip them, since a
problem found there saves hours later.

A Python port of all the MATLAB code is in `python/`, with the same steps as Python
commands in `python/README.md`. It needs no MATLAB or CVX.

## Before you start

- MATLAB R2020a or later (the plot scripts use `exportgraphics`).
- Toolboxes: Statistics and Machine Learning (required), Parallel Computing (optional;
  without it, add `'use_parfor', false` to every override struct below).
- CVX: run `cvx_setup` once and then `savepath`, so that CVX stays on the path.
- Get the latest `main` of this repository.
- If you still have the MI-minimizing baseline solver, save it as
  `fig.7.detection/AN_run_optimization_zou.m` (same arguments and first seven outputs as
  `AN_run_optimization_mi.m`). Without it, that baseline is skipped.
- Every sweep saves its results only when it finishes, so keep the computer awake
  (disable sleep) during the long runs.
- The full runs overwrite `detection_sweep.mat`, `detection_results.mat`,
  `outage_results.mat` and the figure PDFs in `figure/`. The published versions stay in the
  git history.

Every driver takes an optional second argument that overrides settings of
`common/sim_config.m` without editing any file, e.g.
`run_detection_sweep('', struct('num_runs', 1))`. The first argument is the output file
(`''` for the default).

## Step 0: start a log

```matlab
cd('<path to the repository>')
diary('matlab_run_log.txt')        % keeps all printouts; I need them for the letter
```

## Step 1: smoke test (a few minutes)

```matlab
cd fig.7.detection
run_detection_sweep('smoke_detection.mat', ...
    struct('num_runs', 1, 'rho_list', [0.75 1.05], 'use_parfor', false));
```

Check that:

- CVX reports no errors and the run ends with `Saved ... (1 of 1 runs kept ...)`;
- for the proposed design, the solver's line `Minimum beampattern difference
  (constraint C4)` is positive.

If either fails, stop and send me the printout. Afterwards, delete `smoke_detection.mat`.

If you use `parfor` in the later steps and a run fails with an undefined CVX function,
start the pool and run `pctRunOnAll cvx_setup` once, then rerun the step.

## Step 2: check the rewritten scripts against the published results (about 30 min)

The sweep scripts were rewritten from scratch, and some settings in `sim_config.m` are
marked [ASSUMED]. This run uses the published design and virtual location to check that
the new scripts reproduce your published numbers.

```matlab
run_detection_sweep('validate_published.mat', struct('num_runs', 100, ...
    'PosFake', [5 -25], 'sym_null_prop', false, 'fake_tleak', false, ...
    'fake_resolvable_only', false));
make_detection_results('validate_published.mat', 'validate_published_results.mat');
```

Compare the printout with the paper: no-AN success about 35.8 %, proposed about 2.8 % at
η ≈ 0.50, sensing SNR about 14.7 dB (no AN) and 11.1 dB (proposed). With 100 runs, expect
deviations of up to a few percentage points. If the no-AN success is far from 36 %, the
[ASSUMED] settings differ from your original ones: stop and send me the printout together
with your original `Parameters.m` / `LinesCont.m` if you have them.

## Step 3: detection sweep, Figs. 2-3 and three letter figures (several hours)

```matlab
run_detection_sweep();             % 500 realizations -> detection_sweep.mat
make_detection_results();          % -> detection_results.mat and a printout
plot_localization_success();       % Fig. 2 -> figure/localization_success_tvt.pdf
plot_sensing_snr();                % Fig. 3 -> figure/sensing_snr_final.pdf
plot_detection_pd();               % letter, Reviewer 3 C2 -> detection_pd_response.pdf
plot_multi_hypothesis();           % letter, Reviewer 2 C5 -> multi_hypothesis_response.pdf
plot_location_scatter();           % letter Fig. R1, Reviewer 1 C3 -> location_scatter_response.pdf
```

Power-matched operating point for Reviewer 2, Comment 5 (about an hour):

```matlab
run_detection_sweep('detection_sweep_matched.mat', struct('num_runs', 200, ...
    'rho_list', [0.9 1.0 1.05 1.1], 'fake_upper', 1.2));
make_detection_results('detection_sweep_matched.mat', 'detection_results_matched.mat');
```

## Step 4: outage sweep, Fig. 4 and Reviewer 1, Comment 5 (several hours)

```matlab
cd ../fig.9.outage
run_outage_sweep();                % -> outage_results.mat
plot_outage();                     % Fig. 4 -> figure/outage_vs_nmse.pdf and a printout
```

## Step 5: synchronization sweep, Reviewer 2, Comment 4 (under an hour)

```matlab
run_sync_sweep();                  % -> sync_results.mat and a printout
```

## Step 6: communication tradeoff, Reviewer 3, Comment 3 (1-2 hours per run)

```matlab
cd ../fig.7.detection
run_comm_tradeoff();                                               % rho = 0.75 (eta ~ 0.5)
run_comm_tradeoff('comm_tradeoff_rho030.mat', struct('rho_ref', 0.3));   % Table R1 rows
run_comm_tradeoff('comm_tradeoff_rho050.mat', struct('rho_ref', 0.5));
```

## Step 7: finish

```matlab
diary off
```

Commit and push these files to `main` (or upload them):

| Folder | Files |
|---|---|
| repository root | `matlab_run_log.txt` |
| `fig.7.detection/` | `validate_published.mat`, `validate_published_results.mat`, `detection_sweep.mat`, `detection_results.mat`, `detection_sweep_matched.mat`, `detection_results_matched.mat`, `comm_tradeoff*.mat`, `detection_pd_response.pdf`, `multi_hypothesis_response.pdf`, `location_scatter_response.pdf` |
| `fig.9.outage/` | `outage_results.mat`, `sync_results.mat` |
| `Fake-Peak .../figure/` | `localization_success_tvt.pdf`, `sensing_snr_final.pdf`, `outage_vs_nmse.pdf` |

With these, I fill in the [X] values and the published numbers in the reply letter, insert
Fig. R1 and Table R1, and then revise the manuscript.

## Approximate run times

Based on the published runs (1.7 h for the detection sweep and 4.4 h for the outage sweep
with a parallel pool), expect roughly 2-4 h for step 3, 5-6 h for step 4, under 1 h for
step 5, and 1-2 h per run in step 6. Times scale with the number of CPU cores.
