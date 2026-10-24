function run_sync_sweep(out_file, overrides)
%RUN_SYNC_SWEEP Robustness of the effective-channel null to AP phase errors.
%
%   Reviewer 2, Comment 4.  The beams are designed on the LMMSE estimates of
%   Table I (no-AN baseline and proposed nominal design, rho = cfg.rou).
%   On the actual channels, the signal of AP k is rotated by a residual
%   carrier phase psi_k ~ N(0, s^2), i.i.d. over APs and cfg.num_phase
%   draws, for each s in cfg.sync_std_deg.  Timing and frequency offsets
%   within a frame are modeled through the same per-AP phase.  Outputs, in
%   sync_results.mat in this folder:
%     m_out_noan, m_out_nom   outage probability Pr{min_i gamma_i < Gamma}
%     c_out_*                 95% half-widths over the channel realizations
%     m_anr                   mean residual AN power at the combiner outputs
%                             relative to the noise, dB
%     m_anr_theory            the same from the closed form
%                             E|sum_k e^{j psi_k} c_k|^2 = e^{-s^2} |sum_k c_k|^2
%                             + (1 - e^{-s^2}) sum_k |c_k|^2,
%                             with c_k = u_i^H H_{i,k} z_k
%   The sensing SNR (2) is unaffected, since it combines the transmit APs
%   non-coherently, and so is the attack, which works per AP.
%
%   RUN_SYNC_SWEEP(out_file, overrides) writes to out_file ('' for the
%   default) and replaces any cfg setting by the fields of the struct
%   overrides, e.g. struct('num_runs', 1, 'use_parfor', false).

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'common'));
if nargin < 1 || isempty(out_file)
    out_file = fullfile(here, 'sync_results.mat');
end
if nargin < 2
    overrides = struct();
end

cfg = sim_config();
cfg.rou = 1.3;
cfg.sync_std_deg = [0 1 2 5 10 20];
cfg.num_runs = 100;
cfg.num_phase = 1000;
cfg.base_seed = 20260825;
cfg = apply_overrides(cfg, overrides);
geo = setup_geometry(cfg);

R = cfg.num_runs;
ns = numel(cfg.sync_std_deg);
out_noan = nan(R, ns); out_nom = nan(R, ns);
anr = nan(R, ns); anr_th = nan(R, ns);
okr = false(R, 1);
t0 = tic;
if cfg.use_parfor
    parfor r = 1:R
        [out_noan(r, :), out_nom(r, :), anr(r, :), anr_th(r, :), okr(r)] = one_run(cfg, geo, r);
    end
else
    for r = 1:R
        [out_noan(r, :), out_nom(r, :), anr(r, :), anr_th(r, :), okr(r)] = one_run(cfg, geo, r);
    end
end
cfg.elapsed = toc(t0);

n = nnz(okr);
S.cfg = cfg;
S.sync_std_deg = cfg.sync_std_deg;
S.n_kept = n;
S.m_out_noan = mean(out_noan(okr, :), 1);
S.m_out_nom = mean(out_nom(okr, :), 1);
S.c_out_noan = 1.96 * std(out_noan(okr, :), 0, 1) / sqrt(n);
S.c_out_nom = 1.96 * std(out_nom(okr, :), 0, 1) / sqrt(n);
S.m_anr = 10*log10(mean(anr(okr, :), 1));
S.m_anr_theory = 10*log10(mean(anr_th(okr, :), 1));
save(out_file, '-struct', 'S');
fprintf('Saved %s (%d of %d runs kept, %.0f s).\n', out_file, n, R, cfg.elapsed);
for s = 1:ns
    fprintf('  phase std %4.1f deg: outage no AN %.2f %%, proposed %.2f %%; residual AN/noise %.1f dB (theory %.1f dB)\n', ...
        cfg.sync_std_deg(s), 100*S.m_out_noan(s), 100*S.m_out_nom(s), S.m_anr(s), S.m_anr_theory(s));
end
end

function [o_noan, o_nom, anr, anr_th, ok] = one_run(cfg, geo, r)
rng(cfg.base_seed + r);
ch = gen_channels(cfg, geo, 'pilot');
[w0, t0, u0, ~, ok0] = solve(cfg, geo, ch, struct('an_mode', 'none'));
z0 = zeros(cfg.Mtx, cfg.Ntx);
[w1, t1, u1, z1, ok1] = solve(cfg, geo, ch, struct('an_mode', 'fake', ...
    'sym_null', cfg.sym_null_prop, 'fake_tleak', cfg.fake_tleak, ...
    'fake_mask', fake_mask(cfg, geo), 'fake_upper', cfg.fake_upper));
ok = ok0 && ok1;

ns = numel(cfg.sync_std_deg);
o_noan = zeros(1, ns); o_nom = zeros(1, ns);
anr = zeros(1, ns); anr_th = zeros(1, ns);
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;
for s = 1:ns
    sd = cfg.sync_std_deg(s) * pi/180;
    rot = exp(1j * sd * randn(cfg.num_phase, cfg.Ntx));     % draws x Ntx
    o_noan(s) = outage(cfg, ch.H, rot, u0, w0, t0, z0);
    [o_nom(s), res] = outage(cfg, ch.H, rot, u1, w1, t1, z1);
    anr(s) = res.an_power / sigma_ue;
    anr_th(s) = (exp(-sd^2) * res.an_coh + (1 - exp(-sd^2)) * res.an_inc) / sigma_ue;
end
end

function [p_out, res] = outage(cfg, H, rot, u, w_i, w_t, z)
% Outage with per-AP phase rotations; also the residual AN statistics.
Nue = cfg.Nue; Ntx = cfg.Ntx; Mtx = cfg.Mtx;
Gamma = 10^(cfg.Gamma_comm/10);
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;
fail = false(size(rot, 1), 1);
res.an_power = 0; res.an_coh = 0; res.an_inc = 0;
for i = 1:Nue
    % Per-AP contributions of every stream at user i's combiner output.
    C = zeros(Ntx, Nue + 2);
    for k = 1:Ntx
        g = u{i}' * H(:, (k-1)*Mtx + (1:Mtx), i);
        C(k, :) = g * [w_i(:, :, k), w_t(:, k), z(:, k)];
    end
    P = abs(rot * C).^2;                                 % draws x streams
    sinr = P(:, i) ./ (sum(P, 2) - P(:, i) + sigma_ue);
    fail = fail | (sinr < Gamma);
    res.an_power = res.an_power + mean(P(:, end)) / Nue;
    res.an_coh = res.an_coh + abs(sum(C(:, end)))^2 / Nue;
    res.an_inc = res.an_inc + sum(abs(C(:, end)).^2) / Nue;
end
p_out = mean(fail);
end

function [w_i, w_t, u, z, ok] = solve(cfg, geo, ch, opts)
o = cfg.solver_opts;
for f = fieldnames(opts).'
    o.(f{1}) = opts.(f{1});
end
[~, w_i, w_t, u, z, ~, ok] = AN_run_optimization_robust( ...
    cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Nrx, cfg.Mue, cfg.P_T, cfg.Gamma_comm, ...
    cfg.sigma_ue_dBm, cfg.sigma_rx_dBm, cfg.epsilon, ch.Hhat, ...
    geo.Txbeamsteering, geo.Fake_Txbeamsteering, cfg.sigmarcs, geo.betark, cfg.rou, o);
end

function mask = fake_mask(cfg, geo)
% APs at which the fake-peak constraint (9g) is imposed.
if cfg.fake_resolvable_only
    mask = geo.resolvable;
else
    mask = true(1, cfg.Ntx);
end
end
