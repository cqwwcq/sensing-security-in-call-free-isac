# Paper figures ↔ MATLAB code

Paper: `Fake-Peak Artificial Noise for Sensing-Location Privacy in Cell-Free Massive MIMO ISAC/paper_draft.tex`
(TVT correspondence, Figs. 1–4).

The folder names `fig.7.detection` and `fig.9.outage` follow the figure numbering of an
earlier, longer version of the paper (the logs say "Fig. 9/10 sweep" and save to
`figure/detection_rate_final.pdf`). They do **not** match the current figure numbers.

Every MATLAB file listed here has a Python counterpart in `python/`; see the map in
`python/README.md`.

## Summary

| Paper fig. | PDF in `figure/` | Code folder | Data file used for the plotted values | Values checked against the text |
|---|---|---|---|---|
| Fig. 1 | — (TikZ, receive-processing asymmetry) | none; drawn in `paper_draft.tex` | — | — |
| Fig. 2 | `localization_success_tvt.pdf` | `fig.7.detection/` | `detection_results.mat` (aggregated), `detection_sweep.mat` (per realization) | ✅ |
| Fig. 3 | `sensing_snr_final.pdf` | `fig.7.detection/` (same sweep as Fig. 2) | `detection_sweep.mat` (`snr_*`) | ✅ |
| Fig. 4 | `outage_vs_nmse.pdf` | `fig.9.outage/` | `outage_results.mat` | ✅ |

## Fig. 2: localization success rate vs. η, and Fig. 3: sensing SNR vs. η

Both figures come from one sweep over ρ ∈ {0.05, 0.15, 0.3, 0.5, 0.75, 1.05, 1.5, 2.2, 3.5}
(500 channel realizations, 461 kept; 432 kept for the MI-minimizing baseline;
seed 20260823; N = 128 snapshots; 10 m grid; 3604-point angle scan).

| File | Role |
|---|---|
| `AN_run_optimization_mi.m` | Transmit/receive design, problem (9) solved by the CCCP update (12) alternating with the MMSE update (15). `opts.an_mode`: `'fake'` = proposed, `'none'` = no-AN baseline, `'nullonly'` = null-space-only Gaussian-AN baseline. Uses CVX. |
| `AN_update_receiver_mmse_null.m` | Constrained MMSE combiner (15) with Q_iᴴu_i = 0. `rx_mode` 0/1/2 = no null / AN null / AN + sensing-beam null. |
| `AN_Em.m` | Reconstruction attack: per-AP EM estimate x̂_k[n] and the sample covariance R_x̂k of (3); user 1 is the adversary. |
| `detection_sweep.mat` | Per-realization results: `det_*` (success 0/1), `err_*` (localization error), `eta_*` (realized η), `snr_*` (sensing SNR, dB) for `noan`, `null`, `prop`, `zou`. |
| `detection_sweep.mat.bak` | Same content. In `detection_sweep.mat` the four `*_zou` datasets lack the `MATLAB_class` attribute, which MATLAB may need to load them; the `.bak` copy has it. |
| `detection_results.mat` | Fig. 2 plot data: `p_*` success rate, `lo_*`/`hi_*` 95 % Wilson interval, `eta`, `eta_zou`, `cfg`. |
| `det_plot.log` | ⚠️ From a **different, older run** (800 realizations, 744 kept, no MI-minimizing baseline). Its numbers (36.8 %, 3.2 %, 26.5 %) do not match the paper. |

Recomputed from the data files (η ≈ 0.50 operating point is ρ = 0.75, or ρ = 0.5 for the MI baseline):

| Quantity | Paper text | Data |
|---|---|---|
| No-AN success | 35.8 % | 35.79 % (`p_noan`) |
| Proposed / null-only success at η ≈ 0.50 | 2.8 % / 26.0 % | 2.82 % / 26.03 % (η = 0.495) |
| MI-minimizing success at η ≈ 0.50 | 23.6 % | 23.61 % (η = 0.507) |
| Proposed reaches 0 % | η ≈ 0.58 | η = 0.576 |
| MI-minimizing near 0 % | η ≈ 0.78 | 0.69 % at η = 0.780 |
| No-AN sensing SNR | 14.7 dB | 14.72 dB |
| Proposed / null-only / MI sensing SNR at η ≈ 0.50 | 11.1 / 11.4 / 11.5 dB | 11.07 / 11.40 / 11.46 dB |

The sensing SNR in Fig. 3 is the mean of per-realization dB values, not the dB value of the mean linear SNR.

## Fig. 4: communication outage vs. CSI NMSE

300 realizations × 3000 CSI-error draws per NMSE point, NMSE ∈ {−28, −26, …, −4} dB,
ρ = 1.3, seed 20260824.

| File | Role |
|---|---|
| `AN_run_optimization_robust.m` | Same solver as above plus the robust design (20): passing `opts.beta_error` (N_tx × N_UE error variances) adds the CSI-error floor L_i^csi of (19) to the SINR denominator. Modes `'fake'`, `'none'`, `'nullonly'`. |
| `AN_update_receiver_mmse_null.m` | Same as in `fig.7.detection`, extended to take one effective noise level per user (σ_i² + L_i^csi) for the robust design. |
| `beamsteering.m`, `compute_angle_dist.m` | a(φ) and φ_k = ∠(p_t − p_k). |
| `outage_results.mat` | Fig. 4 plot data: `m_noan`, `m_null`, `m_nom`, `m_rob` (outage probability), `c_*` (95 % half-widths), `m_sna`, `m_snu`, `m_sn`, `m_sr` (sensing SNR, dB), `n_kept`, `nmse`, `cfg`. |
| `outage_plot_new.log` | Printout of the plotted values; matches the paper. |
| `outage_sim.log`, `outage_plot.log` | ⚠️ From an **older run** (196 realizations × 4000 draws, nominal and robust only). Not the paper version. The log says it saved `outage_sweep.mat`, which is not in the repo. |

| Quantity | Paper text | Data |
|---|---|---|
| Outage at −10 dB: no-AN / nominal / robust | 75.1 % / 18.6 % / 1.6 % | 75.12 % / 18.64 % / 1.64 % |
| Robust sensing-SNR cost at −10 dB | 0.19 dB | 9.31 − 9.12 = 0.19 dB |
| Outage at −4 dB: robust / others | 31.6 % / > 99 % | 31.64 % / 99.28–99.96 % |

The four curves are no-AN, null-space-only, proposed nominal and proposed robust.

## Regenerating the figures

The original sweep drivers, plot scripts and helpers were not uploaded. They have been
rewritten for the revision from the paper, the solver interfaces and the `cfg` structs
saved in the result files. They have **not** been run with CVX yet, so their numbers
will not match the published ones exactly until the [ASSUMED] settings below are
replaced with the original values.

Run from MATLAB, in this order (step-by-step commands, checks and run times are in
[`RUN_SIMULATIONS.md`](RUN_SIMULATIONS.md)):

| Step | Script | Output |
|---|---|---|
| 1 | `fig.7.detection/run_detection_sweep.m` | `detection_sweep.mat` |
| 2 | `fig.7.detection/make_detection_results.m` | `detection_results.mat`, printout of the quoted values |
| 3 | `fig.7.detection/plot_localization_success.m` | `figure/localization_success_tvt.pdf` (Fig. 2) |
| 4 | `fig.7.detection/plot_sensing_snr.m` | `figure/sensing_snr_final.pdf` (Fig. 3) |
| 5 | `fig.7.detection/plot_detection_pd.m` | `fig.7.detection/detection_pd_response.pdf` (response letter, R3-2) |
| 5b | `fig.7.detection/plot_multi_hypothesis.m` | `fig.7.detection/multi_hypothesis_response.pdf` (response letter, R2-5) |
| 5c | `fig.7.detection/plot_location_scatter.m` | `fig.7.detection/location_scatter_response.pdf` (response letter Fig. R1, R1-3) |
| 6 | `fig.9.outage/run_outage_sweep.m` | `outage_results.mat` |
| 7 | `fig.9.outage/plot_outage.m` | `figure/outage_vs_nmse.pdf` (Fig. 4), printout |
| 8 | `fig.9.outage/run_sync_sweep.m` | `sync_results.mat`, printout (response letter, R2-4) |
| 9 | `fig.7.detection/run_comm_tradeoff.m` | `comm_tradeoff.mat`, printout (response letter, R3-3) |

Shared helpers are in `common/`: `sim_config.m` (all parameters), `setup_geometry.m`,
`gen_channels.m` (LMMSE pilots or fixed NMSE), `localize_target.m` (dominant-peak line
fusion with cell voting), `localize_multi.m` (multi-hypothesis attack), `ray_cells.m`,
`draw_nullonly_an.m`, `sensing_metrics.m`, `detection_pd.m`,
`wilson_ci.m`, `paper_fig_style.m`, and copies of `beamsteering.m` / `compute_angle_dist.m`.

Additions for the revision, beyond the published sweep:

- **AN echo exploited** (Reviewer 2, Comment 1): `snrz_*`, plotted dashed in Fig. 3.
- **Revised fake-peak constraint** (Reviewer 1, Comment 2): both solvers take
  `opts.fake_tleak` (adds T_k = |a^H(φ_k) z_k|² to the right-hand side of (9g), so ρ > 1
  guarantees the fake-direction response exceeds the true-direction one) and
  `opts.fake_mask` (APs at which (9g) is imposed). `sim_config.m` turns both on
  (`cfg.fake_tleak`, `cfg.fake_resolvable_only`); the mask keeps the APs whose true and fake
  directions are at least one Rayleigh width apart, |sin φ_k − sin φ_fake,k| ≥ 2/M_tx. With the
  relocated virtual location `cfg.PosFake = [175 35]` all 8 APs pass (min |Δ sin| = 0.163), so
  the mask is all-true; with the published (5, −25), the APs at (−200, 200), (200, −200) and
  (0, −200) fail, being nearly collinear with the target and the fake location. Set both
  options to `false` and `cfg.PosFake = [5 -25]` to reproduce the published sweep.
- **Fake-peak dominance** (Reviewer 1, Comment 2): `fdom_*`, the fraction of APs whose
  reconstructed fake-direction response exceeds the true-direction response, and `pkf_*` /
  `pkt_*`, the fraction whose reconstructed dominant peak lies at the fake / true direction,
  all over the APs with resolvable directions.
- **Multi-hypothesis adversary** (Reviewer 2, Comment 5; Reviewer 1, Comment 3):
  `localize_multi.m` keeps the two strongest distinct peaks per AP, forms the two most-voted,
  separated location hypotheses, and resolves them by votes, at random, or mechanism-aware
  (weaker peaks win); `mhv_*`, `mhr_*`, `mha_*` hold the three success indicators and `mhs_*`
  the two-guess success. `est_*` stores the dominant-peak estimates for the scatter plot and
  `fhit_*` whether they land in the virtual location's cell. Each AP votes at most once per
  cell, and cells next to a TX AP are not taken as hypotheses (`cfg.mh_ap_guard`); a second
  hypothesis competes only with at least `cfg.mh_coherent` of the first one's votes. The fake
  hypothesis is supported only by APs with resolvable directions, so with the published
  `PosFake = (5, −25)` (5 of 8 resolvable) a vote-counting adversary recovers the target; this
  is why the default is now (175, 35), which all 8 APs resolve.
- **Fake-peak cap** (Reviewer 2, Comment 5): `opts.fake_upper` (κ) in both solvers caps the
  fake-direction response at κ(ρ‖a^H W_k‖² + T_k), so ρ ≈ 1 gives power-matched true and fake
  peaks. Off by default (`cfg.fake_upper = Inf`); for the matched operating point, rerun with
  e.g. `cfg.fake_upper = 1.2` and `cfg.rho_list` around 1.
- **AP synchronization** (Reviewer 2, Comment 4): `run_sync_sweep.m` rotates each AP's signal
  by a residual phase ψ_k ~ N(0, s²), s ∈ {0, 1, 2, 5, 10, 20}°, on the actual channels and
  reports the outage of the no-AN and proposed designs and the residual AN power at the
  combiner outputs, next to the closed form (1 − e^{−s²}) Σ_k |u_i^H H_{i,k} z_k|² (plus the
  e^{−s²}-weighted coherent residual from CSI errors).
- **Communication cost of the AN** (Reviewer 3, Comment 3): `run_comm_tradeoff.m` gives the
  no-AN design the proposed design's sensing SNR at ρ = 0.75 (η ≈ 0.5) and finds by bisection
  the largest common SINR target it still meets; it reports that target and the sum rates
  (via `common/user_sinr.m`). Tradeoff rows for other η need `cfg.rho_ref` changed and a rerun.
- **More favorable adversary** (Reviewer 3, Comment 4): the detection sweep also evaluates the
  dominant-peak attack with the search cells of `cfg.cell_list` = 10/20/50 m (`detc_*`).
- **Detection probability** (Reviewer 3, Comment 2): `pd_*`, GLRT over the Nue+1 known
  streams with Swerling-I echo draws on the designed beams, for N ∈ {5, 128} and RCS
  offsets −30…0 dB. `plot_detection_pd.m` overlays the closed-form model quoted in the
  response letter, as a check of those numbers.
- **Robust no-AN baseline** (Reviewer 1, Comment 5): `m_rnoan` / `m_srn` in the outage sweep.
- **Fixes from the points below**: (9f) is enforced in the proposed design
  (`cfg.sym_null_prop = true`; set `false` to match the published sweep), and the adversary's
  noise is drawn independently per snapshot.

Settings marked [ASSUMED] in `common/sim_config.m` are not stated in the paper: the path-loss
exponent 3 and RCS 10 dBsm (from the earlier multi-static paper in this repo), the UE positions,
the AP ordering, the 400 m × 400 m search area, and how the ULA front-back ambiguity is resolved
(`cfg.mirror_rule`). The no-AN localization success (35.8 % in the paper) is the most sensitive
check of these.

Still missing:

- **MI-minimizing (Zou et al.) baseline.** No solver mode implements it:
  `AN_run_optimization_mi.m` only has `'fake'`, `'none'` and `'nullonly'`, despite its name.
  `run_detection_sweep.m` calls `AN_run_optimization_zou.m` (same arguments and first seven
  outputs as `AN_run_optimization_mi.m`) if it is on the path, and skips the baseline otherwise.
- MATLAB toolboxes needed: CVX for the solvers; Statistics (`gamrnd`, `normrnd`, `chi2inv`,
  `ncx2cdf`); Parallel Computing for the parallel sweeps (`cfg.use_parfor`).

## Points to check against the paper

1. **Sensing-beam null (9f).** In `'fake'` mode `SYM_NULL` defaults to `false`, so (9f) is
   enforced only if the driver passes `opts.sym_null = true`. The paper states that the
   proposed design enforces (9f). The constraint check printed at the end of the solver
   (`[Check] Null-space leakage (constraint C3)`) covers the AN only.
2. **Per-AP separation is built into the attack.** `AN_Em.m` forms each AP's observation
   `H_{a,k}(W_k s + z_k s_z) + noise` separately, so the adversary gets perfectly separated
   per-AP signals. This is the assumption Reviewer 1 (comment 1) and Reviewer 2 (comment 2) question.
3. **Adversary CSI.** `AN_Em.m` computes its initial least-squares estimate with the true
   `H_{a,k}` (`pinv(H'H)H'y`) and uses the true `H_{a,k}` as the prior mean of the EM channel
   update, so the attack effectively has perfect CSI. The paper describes the adversary using
   its own estimate Ĥ_{a,k}^adv.
4. **Noise in the attack.** `AN_Em.m` adds its `sigma2` argument to every snapshot. In the
   earlier `LinesCont.m` caller, `sigma2` was one M_UE × 1 noise vector drawn per realization,
   i.e. the same noise sample on all 128 snapshots and all APs. Check that the current driver
   passes independent noise per snapshot.
5. **Null-space-only baseline.** The beams are optimized over the budget left after the AN
   power, and the AN is drawn afterwards by projecting a random vector onto the null space of
   the returned combiners. That matches the paper's description, but is worth stating in the
   setup since η is then prescribed rather than optimized.
