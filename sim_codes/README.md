# Simulation code ↔ paper figure map

The folder names (`fig.2`, `fig3`, …, `fig.8`) follow the figure numbering of an
earlier, longer version of the paper. They do **not** match the figure numbers in
the current TVT correspondence draft (`paper_draft.tex`, Figs. 1–4).

All shared functions (`Parameters.m`, `run_optimization.m`, `AN_run_optimization.m`,
`update_receiver_mmse.m`, `AN_update_receiver_mmse.m`, `Em.m`, `AN_Em.m`,
`beamsteering.m`, `compute_angle_dist.m`, `signalgeneration.m`,
`generate_projected_AN.m`) are byte-identical copies across the folders that contain them.

## 1. Figures in the current draft (`paper_draft.tex`)

| Fig. | File | Produced by | Status |
|---|---|---|---|
| 1 | TikZ (receive-processing asymmetry) | `paper_draft.tex`, lines 99–121 | No MATLAB needed |
| 2 | `figure/localization_success_tvt.pdf` (success rate vs. η; No AN / null-space-only Gaussian AN / MI-minimizing AN / proposed) | **Not in `sim_codes`** | Closest predecessor: `fig5/` (detection rate vs. AN-power ratio, Gaussian vs. fake-peak). Missing: MI-minimizing (Zou et al.) baseline, N_UE = 4, P_T = 30 dBm, LMMSE CSI, 10 m grid |
| 3 | `figure/sensing_snr_final.pdf` (sensing SNR vs. η, same four schemes) | **Not in `sim_codes`** | Closest predecessor: `SNRvsPowerRatio.pdf` (single curve, 14.7 dB at η = 0, matching the no-AN value in the text). γ_t is only printed by `AN_run_optimization.m`, never saved |
| 4 | `figure/outage_vs_nmse.pdf` (outage vs. CSI NMSE; No AN / null-space-only / proposed nominal / proposed robust) | **Not in `sim_codes`** | Nothing here models pilot-based LMMSE estimation, an NMSE sweep, or the robust design (17)–(20); all designs use perfect `H` |

## 2. Legacy figure folders

| Folder | What it shows | Data generator | Plot script | PDF in paper folder |
|---|---|---|---|---|
| `fig.2` | Per-AP directional MI I(Θ;Θᵖ), without vs. with fake-peak AN (bar chart) | `run_MI_direction_detection_compare.m` (ρ = 1.3) → `result.MI.mat` | `picture_MI.m` (values typed in); `createfigure.m` is an unused auto-generated helper; `figure2.fig` | none |
| `fig3` | System layout | — | `picture_systemmodel.m` | `system_layout.pdf` |
| `fig4` | Normalized transmit beampattern of AP 7, conventional vs. proposed | `system_simulation.m` (5 runs; `ki=1` no-AN → `*MC2.mat`, `ki=2` proposed ρ = 0.5977 → `*MC0.mat`) | `beampattern_picture.m`; `untitled.fig` | none (old builds used `figure/beampattern_stacked.pdf`, `fig3a.pdf`, `figure4.pdf`, not in repo) |
| `fig5` | Detection (= localization-success) rate vs. AN-power ratio; fake-peak vs. Gaussian AN, M_tx = 16 and 8 | `system_simulation.m` (ρ = 0…1.5) → `LinesCont.m` (`ki=1` Gaussian AN via `generate_projected_AN.m`, `ki=2` fake-peak) | `picture_gaussAN_VS_fakeAN.m` (values typed in, smoothed); `figure5.fig` | `guassANvsFakeAN.pdf` |
| `fig.6` | Detection rate vs. AN-power ratio for M_UE = 16/12/8/4 | `system_simulation.m` (ρ = 0…2.5, 100 runs, `Mue` changed by hand) → `LinesCont.m` | `picture_diffMue.m` (values typed in, smoothed) | `diffMue.pdf` |
| `fig.7` | CDF of the minimum user SINR: no AN / Gaussian AN / fake-peak AN | `system_simulation.m` (1000 runs; ρ = 0 and 0.5977 with `run_optimization.m`, ρ = 1.3 with `AN_run_optimization.m`) → `sinr_*.mat` | `picture_sinr_cdf.m`; `untitled.fig` | `SINR.pdf` |
| `fig.8` | CDF of localization error, 5 schemes | **missing** (only `final_Poscont.mat`, 1000×2×5×10; `LinesCont.m` saves a different layout) | `picture_CDF.m`; `untitled.fig` | `loc-error-cdf.pdf` |
| — | Sensing SNR vs. AN-power ratio | **missing** | **missing** | `SNRvsPowerRatio.pdf` |
| — | CCCP convergence | **missing** | **missing** | `convergence.pdf` |

## 3. Functions ↔ equations in `paper_draft.tex`

| Function | Paper element | Notes |
|---|---|---|
| `Parameters.m` | System model, Table I geometry, g_{r,k} in (2) (`betark`), H_{i,k} | β_{i,k}: COST-231 Hata three-slope (h_BS = 15 m, h_UE = 1.65 m) + correlated 8 dB shadowing. `sigmarcs = sqrt(10)` → RCS σ² = 10 m² (= ξ² in (2)) |
| `beamsteering.m`, `compute_angle_dist.m` | a(φ) = [1, e^{jπ sin φ}, …]ᵀ, φ_k = ∠(p_t − p_k) | Global azimuth from `atan2` |
| `signalgeneration.m` | s_i[n], s_t[n], N = 128 snapshots | 16-QAM data symbols have power 1/N_UE, not 1 |
| `run_optimization.m` | No-AN design (problem (9) without z, (9e), (9g)) | Per-AP budget (1−ρ)P_T; with ρ > 0 it is the power-reserved design for the Gaussian-AN baseline |
| `AN_run_optimization.m` | Proposed CCCP transmit update (12) and its alternation with the MMSE receive update | Objective, (9d) and linearized (9g) as in the paper; ε_tol = 4·10⁻³, 30 iterations max; initialization as in Sec. III-B |
| `update_receiver_mmse.m`, `AN_update_receiver_mmse.m` | Receive combiner update, (14)–(16) | Unconstrained MMSE R⁻¹q normalized; no Q_iᴴu_i = 0 constraint |
| `generate_projected_AN.m` | Null-space-only Gaussian AN baseline | Per-AP projection onto null{u_iᴴH_{i,k}}, power ρP_T |
| `Em.m`, `AN_Em.m` | Reconstruction attack: x̂_k[n] and R_x̂k in (3) | User 1 is the adversary; its own stream is removed |
| `LinesCont.m` | Dominant-peak direction lines, grid voting, success = estimate in target cell | 3604-point scan over [0, 2π]; 20 m cells; 10 seeds × 100 runs |

## 4. Where code and current draft differ

| Item | `paper_draft.tex` | `sim_codes` |
|---|---|---|
| N_UE | 4 | 2 |
| P_T | 30 dBm | 50 dBm |
| Grid cell Δ_g | 10 m | 20 m |
| CSI | LMMSE estimates from τ_p = 64 pilots, P_p = 20 dBm | Perfect `H` everywhere |
| AN nulling | Network-wide Σ_k ĝᴴ_{i,k} z_k = 0, (6)/(9e) | Per-AP \|u_iᴴH_{i,k}z_k\| ≤ 10⁻¹⁰ for every k |
| Sensing-beam nulling (9f) | Present | Absent |
| MMSE combiner (15) | Constrained (Q_iᴴu_i = 0) | Unconstrained |
| Robust design (17)–(20) | Present | Absent |
| MI-minimizing baseline [Zou et al.] | Present | Absent |
| Reported γ_t | Excludes z | `gamma_t_opt` includes the AN echo; the value without z is only printed (`val_noz`) |
| Receive-AP antennas / receive gain | Absorbed into ξ² | Not modeled (`Rxbeamsteering` is unused) |

Matching items: transmit/receive AP and target/fake coordinates, f_c = 1.9 GHz,
noise powers (−94 dBm), M_tx = M_UE = 16, N_tx = 8, N_rx = 4, Γ_comm = 10 dB, N = 128.

## 5. Issues found while tracing the code

1. **Per-AP separation is built into the attack.** `AN_Em.m` forms
   `y = H_{a,k}(W_k s + z_k s_z) + n` for each AP k separately, so the adversary gets perfectly
   separated per-AP observations. This is exactly what Reviewer 1 (comment 1) and Reviewer 2
   (comment 2) question. The adversary's channel prior is the true `H` plus N(0, 10⁻⁴) error,
   which amounts to near-perfect adversary CSI.
2. **Front-back ambiguity in the attack.** Because a(θ) depends only on sin θ, θ and π − θ give
   identical beampatterns. `LinesCont.m` scans [0, 2π] and takes the first maximum, so some APs
   draw the mirrored line. With ideal noiseless pencil beams toward the target, 3 of 8 APs
   (AP 3, 5, 6) get the mirrored direction. This likely contributes to the low no-AN success
   rate that Reviewer 3 (comment 4) points out. Arrays oriented so the area lies on one side, or
   an adversary that tests both mirror candidates, would remove it.
3. **Noise is not white over snapshots.** `sigma2` is one M_UE×1 vector drawn once per seed and
   added unchanged to all 128 snapshots, all APs and all runs.
4. **Wrong file in `fig4/beampattern_picture.m`.** Both `W1` and `W2` load `WMC0.mat`, so the
   "conventional" curve is the proposed design without z. It should load `WMC2.mat`.
5. **Saved data doesn't match the figures.** `fig5/*.mat` hold a single run
   (`num_runs = 1`, 16 ρ values), while `LinesCont.m` loops over 100 runs. These files are test
   output, not the data behind `guassANvsFakeAN.pdf`.
6. **Hard-coded Windows paths.** `LinesCont.m` loads from `C:\Users\pc\Desktop\代码\AN_simulation\fig.5\…`
   or `…\fig.6\…`, and `beampattern_picture.m` from `D:\研究生毕业论文\…`.
7. **Plot scripts don't read the simulation output.** `picture_MI.m`,
   `picture_gaussAN_VS_fakeAN.m` and `picture_diffMue.m` contain typed-in numbers, and the latter
   two apply `smoothdata(…, 'gaussian', 4)` before plotting.
