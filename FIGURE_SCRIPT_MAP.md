# Paper figures ↔ MATLAB code

Paper: `Fake-Peak Artificial Noise for Sensing-Location Privacy in Cell-Free Massive MIMO ISAC/paper_draft.tex`
(TVT correspondence, Figs. 1–4).

The folder names `fig.7.detection` and `fig.9.outage` follow the figure numbering of an
earlier, longer version of the paper (the logs say "Fig. 9/10 sweep" and save to
`figure/detection_rate_final.pdf`). They do **not** match the current figure numbers.

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

## Missing to regenerate the figures

None of the figures can be regenerated from the repo alone. The uploaded files are the
solver functions and the saved results; these pieces are not included:

- **Sweep drivers** that create `detection_sweep.mat` / `detection_results.mat` and
  `outage_results.mat` (channel and position generation, LMMSE estimation, CSI-error
  draws, the ρ and NMSE loops, outage evaluation).
- **Plot scripts** that write `localization_success_tvt.pdf`, `sensing_snr_final.pdf` and
  `outage_vs_nmse.pdf`.
- **MI-minimizing (Zou et al.) baseline.** No solver mode implements it:
  `AN_run_optimization_mi.m` only has `'fake'`, `'none'` and `'nullonly'`, despite its name.
- **Helpers used by Figs. 2–3**: geometry and large-scale fading (formerly `Parameters.m`),
  symbol generation (formerly `signalgeneration.m`), dominant-peak line fusion and grid
  voting (formerly `LinesCont.m`), and `beamsteering.m` / `compute_angle_dist.m` (present only
  in `fig.9.outage/`).
- MATLAB toolboxes needed: CVX for the solvers; Statistics (`gamrnd`, `normrnd`) for `AN_Em.m`;
  Parallel Computing for the parallel sweeps.

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
