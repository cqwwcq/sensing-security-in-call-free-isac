function run_comm_tradeoff(out_file)
%RUN_COMM_TRADEOFF What the AN power costs in communication (Reviewer 3, C3).
%
%   For each channel realization, the proposed design is solved at
%   rho = cfg.rho_ref, the operating point eta ~ 0.5 of Figs. 2-3.  The
%   no-AN design is then given the same sensing SNR and as much
%   communication as possible: the largest SINR target Gamma for which the
%   no-AN solver still reaches the proposed design's sensing SNR is found by
%   bisection.  The gap between that Gamma and the 10 dB requirement is the
%   communication the AN power could buy instead.  Outputs, in
%   comm_tradeoff.mat in this folder, per kept realization:
%     gmax_noan_dB   largest common SINR target of the no-AN alternative
%     rate_noan      its sum rate, bit/s/Hz, at that target
%     rate_prop      sum rate of the proposed design
%     snr_prop_dB    sensing SNR of the proposed design (2)
%     eta_prop       realized AN-power ratio of the proposed design

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'common'));
if nargin < 1
    out_file = fullfile(here, 'comm_tradeoff.mat');
end

cfg = sim_config();
cfg.rho_ref = 0.75;
cfg.gamma_hi_dB = cfg.Gamma_comm + 30;   % upper end of the bisection
cfg.bisect_iter = 7;                     % 30 dB / 2^7, about 0.23 dB
cfg.num_runs = 100;
cfg.base_seed = 20260826;
geo = setup_geometry(cfg);

R = cfg.num_runs;
res = cell(R, 1);
t0 = tic;
if cfg.use_parfor
    parfor r = 1:R
        res{r} = one_run(cfg, geo, r);
    end
else
    for r = 1:R
        res{r} = one_run(cfg, geo, r);
    end
end
cfg.elapsed = toc(t0);
res = [res{:}];
ok = [res.ok];
res = res(ok);

S.cfg = cfg;
S.n_kept = nnz(ok);
for f = {'gmax_noan_dB', 'rate_noan', 'rate_prop', 'snr_prop_dB', 'eta_prop'}
    S.(f{1}) = [res.(f{1})];
end
save(out_file, '-struct', 'S');
fprintf('Saved %s (%d of %d runs kept, %.0f s).\n', out_file, S.n_kept, R, cfg.elapsed);
fprintf('Proposed (rho = %.2f): eta %.3f, sensing SNR %.2f dB, sum rate %.2f bit/s/Hz\n', ...
    cfg.rho_ref, mean(S.eta_prop), mean(S.snr_prop_dB), mean(S.rate_prop));
fprintf('No-AN alternative at the same sensing SNR: common SINR %.2f dB (+%.2f dB), sum rate %.2f bit/s/Hz\n', ...
    mean(S.gmax_noan_dB), mean(S.gmax_noan_dB) - cfg.Gamma_comm, mean(S.rate_noan));
end

function out = one_run(cfg, geo, r)
rng(cfg.base_seed + r);
ch = gen_channels(cfg, geo, 'pilot');

opts = struct('an_mode', 'fake', 'sym_null', cfg.sym_null_prop, ...
    'fake_tleak', cfg.fake_tleak, 'fake_upper', cfg.fake_upper);
if cfg.fake_resolvable_only
    opts.fake_mask = geo.resolvable;
end
[w_i, w_t, u, z, ok] = solve(cfg, geo, ch, cfg.Gamma_comm, cfg.rho_ref, opts);
m = sensing_metrics(cfg, geo, w_i, w_t, z);
out.ok = ok;
out.snr_prop_dB = 10*log10(m.gamma_t);
out.eta_prop = m.eta;
out.rate_prop = sum(log2(1 + user_sinr(cfg, ch.Hhat, u, w_i, w_t, z)));

% Bisection on the SINR target of the no-AN design at the same sensing SNR.
lo = cfg.Gamma_comm;
hi = cfg.gamma_hi_dB;
best = [];
for it = 1:cfg.bisect_iter
    mid = (lo + hi) / 2;
    [wi0, wt0, u0, z0, ok0] = solve(cfg, geo, ch, mid, 0, struct('an_mode', 'none'));
    m0 = sensing_metrics(cfg, geo, wi0, wt0, z0);
    if ok0 && m0.gamma_t >= m.gamma_t
        lo = mid;
        best = struct('w_i', wi0, 'w_t', wt0, 'u', {u0});
    else
        hi = mid;
    end
end
if isempty(best)
    % The bisection never moved up: fall back to the 10 dB target.
    [wi0, wt0, u0, ~, ok0] = solve(cfg, geo, ch, lo, 0, struct('an_mode', 'none'));
    out.ok = out.ok && ok0;
    best = struct('w_i', wi0, 'w_t', wt0, 'u', {u0});
end
out.gmax_noan_dB = lo;
zero = zeros(cfg.Mtx, cfg.Ntx);
out.rate_noan = sum(log2(1 + user_sinr(cfg, ch.Hhat, best.u, best.w_i, best.w_t, zero)));
end

function [w_i, w_t, u, z, ok] = solve(cfg, geo, ch, Gamma_dB, rho, opts)
o = cfg.solver_opts;
for f = fieldnames(opts).'
    o.(f{1}) = opts.(f{1});
end
[~, w_i, w_t, u, z, ~, ok] = AN_run_optimization_mi( ...
    cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Nrx, cfg.Mue, cfg.P_T, Gamma_dB, ...
    cfg.sigma_ue_dBm, cfg.sigma_rx_dBm, cfg.epsilon, ch.Hhat, ...
    geo.Txbeamsteering, geo.Fake_Txbeamsteering, cfg.sigmarcs, geo.betark, rho, o);
end
