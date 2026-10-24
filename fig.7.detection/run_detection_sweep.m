function run_detection_sweep(out_file)
%RUN_DETECTION_SWEEP Sweep behind Figs. 2-3: localization success and SNR vs eta.
%
%   RUN_DETECTION_SWEEP() runs cfg.num_runs channel realizations.  In each,
%   the no-AN baseline is solved once, and for every rho in cfg.rho_list the
%   proposed design ('fake'), the null-space-only baseline ('nullonly',
%   paired with the proposed design's realized eta) and, if available, the
%   MI-minimizing baseline are solved on the LMMSE estimates.  The
%   reconstruction attack (AN_Em.m) and the localization then run on the
%   actual channels.  Results go to detection_sweep.mat in this folder.
%
%   New outputs for the revision, besides the variables of the published
%   sweep (det_*, err_*, eta_*, snr_*):
%     snrz_*     sensing SNR (dB) with the known AN echo exploited
%     fdom_*     fraction of APs whose reconstructed fake-direction response
%                exceeds the true-direction response (Reviewer 1, Comment 2)
%     pkf_*, pkt_*  fraction of APs whose reconstructed dominant peak lies at
%                the fake / true direction (the rest lie elsewhere); fdom_*,
%                pkf_* and pkt_* count only APs with resolvable directions
%     mhv_*, mhr_*, mha_*, mhs_*  success of the multi-hypothesis attack
%                with the votes / random / mechanism-aware rule, and with
%                two guesses (Reviewer 2, Comment 5), see LOCALIZE_MULTI
%     fhit_*     the dominant-peak estimate lies in the virtual location's cell
%     est_*      estimated location of the dominant-peak attack, 2 x runs
%                (or nRho x 2 x runs), for the scatter plot of Reviewer 1,
%                Comment 3
%     pd_*       detection probability, nOffset x nN x runs (or nRho x
%                nOffset x nN x runs), see DETECTION_PD
%
%   The MI-minimizing baseline needs AN_run_optimization_zou.m on the path,
%   with the same arguments and first seven outputs as
%   AN_run_optimization_mi.m.  It is not in the repository; without it the
%   baseline is skipped.
%
%   Requires CVX, the Statistics and Machine Learning Toolbox, and, for
%   cfg.use_parfor, the Parallel Computing Toolbox with CVX set up on the
%   workers.

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'common'));
if nargin < 1
    out_file = fullfile(here, 'detection_sweep.mat');
end

cfg = sim_config();
cfg.num_runs = 500;
cfg.base_seed = 20260823;
assert(cfg.N_snapshot == 128, 'AN_Em.m is hard-coded for 128 snapshots.');
assert(cfg.adversary == 1, 'AN_Em.m treats user 1 as the adversary.');
geo = setup_geometry(cfg);

has_zou = exist('AN_run_optimization_zou', 'file') == 2;
if ~has_zou
    warning('AN_run_optimization_zou.m not found: skipping the MI-minimizing baseline.');
end

R = cfg.num_runs;
res = cell(R, 1);
t0 = tic;
if cfg.use_parfor
    parfor r = 1:R
        res{r} = one_run(cfg, geo, r, has_zou);
    end
else
    for r = 1:R
        res{r} = one_run(cfg, geo, r, has_zou);
    end
end
cfg.elapsed = toc(t0);

res = [res{:}];
ok = [res.ok];
ok_zou = ok & [res.ok_zou];
cfg.num_ok = nnz(ok);
cfg.num_ok_zou = nnz(ok_zou);

S = struct('cfg', cfg);
for name = {'noan', 'null', 'prop'}
    S = collect(S, res(ok), name{1});
end
if has_zou
    S = collect(S, res(ok_zou), 'zou');
end
save(out_file, '-struct', 'S', '-v7.3');
fprintf('Saved %s (%d of %d runs kept, %.0f s).\n', out_file, cfg.num_ok, R, cfg.elapsed);
end

function S = collect(S, res, name)
% Stack the per-run results of one design along the last dimension.
x = [res.(name)];
fields = {'det', 'err', 'eta', 'snr', 'snrz', 'fdom', 'pkf', 'pkt', ...
    'mhv', 'mhr', 'mha', 'mhs', 'fhit'};
for f = fields
    S.([f{1} '_' name]) = cat(2, x.(f{1}));
end
S.(['pd_' name]) = cat(ndims(x(1).pd) + 1, x.pd);
S.(['est_' name]) = cat(3, x.est);
if strcmp(name, 'noan')
    S = rmfield(S, {'eta_noan', 'snrz_noan', 'fdom_noan'});
end
end

function out = one_run(cfg, geo, r, has_zou)
rng(cfg.base_seed + r);
ch = gen_channels(cfg, geo, 'pilot');
nr = numel(cfg.rho_list);
Ns = cfg.Nue + 1;
P_T_lin = 10^(cfg.P_T/10) / 1000;
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;

% Unit-power QPSK symbols (rows 1..Nue: data, row Nue+1: sensing) and
% adversary receiver noise, independent across snapshots.
Dsn = (sign(randn(Ns, cfg.N_snapshot)) + 1j*sign(randn(Ns, cfg.N_snapshot))) / sqrt(2);
noise_adv = sqrt(sigma_ue/2) * (randn(cfg.Mue, cfg.N_snapshot) + 1j*randn(cfg.Mue, cfg.N_snapshot));

solve = @(P_T_dBm, rho, opts) solve_design(cfg, geo, ch, P_T_dBm, rho, opts);
evaluate = @(w_i, w_t, z) evaluate_design(cfg, geo, ch, Dsn, noise_adv, w_i, w_t, z);

out.ok = true;
out.ok_zou = has_zou;

% --- No-AN baseline ---
[w_i, w_t, ~, z, okd] = solve(cfg.P_T, 0, struct('an_mode', 'none'));
out.ok = out.ok && okd;
out.noan = evaluate(w_i, w_t, zeros(size(z)));

blank = struct('det', nan, 'err', nan, 'eta', nan, 'snr', nan, 'snrz', nan, ...
    'fdom', nan, 'pkf', nan, 'pkt', nan, ...
    'mhv', nan, 'mhr', nan, 'mha', nan, 'mhs', nan, 'fhit', nan, 'est', nan(1, 2), ...
    'pd', nan(numel(cfg.pd_offset_dB), numel(cfg.pd_N_list)));
prop = repmat(blank, nr, 1); null = prop; zou = prop;

for j = 1:nr
    rho = cfg.rho_list(j);

    % --- Proposed fake-peak design ---
    opts = struct('an_mode', 'fake', 'sym_null', cfg.sym_null_prop, ...
        'fake_tleak', cfg.fake_tleak, 'fake_mask', fake_mask(cfg, geo), ...
        'fake_upper', cfg.fake_upper);
    [w_i, w_t, ~, z, okd] = solve(cfg.P_T, rho, opts);
    out.ok = out.ok && okd;
    prop(j) = evaluate(w_i, w_t, z);

    % --- Null-space-only Gaussian AN, paired at the same eta ---
    eta = prop(j).eta;
    P_beam_dBm = 10*log10(max(1 - eta, 1e-6) * P_T_lin * 1000);
    [w_i, w_t, u, ~, okd] = solve(P_beam_dBm, rho, struct('an_mode', 'nullonly'));
    out.ok = out.ok && okd;
    z = draw_nullonly_an(cfg, ch.Hhat, u, w_i, w_t, eta * P_T_lin * ones(1, cfg.Ntx));
    null(j) = evaluate(w_i, w_t, z);

    % --- MI-minimizing baseline (adapted from Zou et al.) ---
    if has_zou
        [w_i, w_t, ~, z, okd] = solve_zou(cfg, geo, ch, rho);
        out.ok_zou = out.ok_zou && okd;
        zou(j) = evaluate(w_i, w_t, z);
    end
end

out.prop = stack(prop);
out.null = stack(null);
out.zou = stack(zou);
end

function mask = fake_mask(cfg, geo)
% APs at which the fake-peak constraint (9g) is imposed.
if cfg.fake_resolvable_only
    mask = geo.resolvable;
else
    mask = true(1, cfg.Ntx);
end
end

function s = stack(x)
% Per-rho results as column vectors; pd as nRho x nOffset x nN.
s.det = [x.det].';
s.err = [x.err].';
s.eta = [x.eta].';
s.snr = [x.snr].';
s.snrz = [x.snrz].';
s.fdom = [x.fdom].';
s.pkf = [x.pkf].';
s.pkt = [x.pkt].';
s.mhv = [x.mhv].';
s.mhr = [x.mhr].';
s.mha = [x.mha].';
s.mhs = [x.mhs].';
s.fhit = [x.fhit].';
s.est = cat(1, x.est);
s.pd = permute(cat(3, x.pd), [3 1 2]);
end

function [w_i, w_t, u, z, ok] = solve_design(cfg, geo, ch, P_T_dBm, rho, opts)
o = cfg.solver_opts;
for f = fieldnames(opts).'
    o.(f{1}) = opts.(f{1});
end
[~, w_i, w_t, u, z, ~, ok] = AN_run_optimization_mi( ...
    cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Nrx, cfg.Mue, P_T_dBm, cfg.Gamma_comm, ...
    cfg.sigma_ue_dBm, cfg.sigma_rx_dBm, cfg.epsilon, ch.Hhat, ...
    geo.Txbeamsteering, geo.Fake_Txbeamsteering, cfg.sigmarcs, geo.betark, rho, o);
end

function [w_i, w_t, u, z, ok] = solve_zou(cfg, geo, ch, rho)
[~, w_i, w_t, u, z, ~, ok] = AN_run_optimization_zou( ...
    cfg.Ntx, cfg.Nue, cfg.Mtx, cfg.Nrx, cfg.Mue, cfg.P_T, cfg.Gamma_comm, ...
    cfg.sigma_ue_dBm, cfg.sigma_rx_dBm, cfg.epsilon, ch.Hhat, ...
    geo.Txbeamsteering, geo.Fake_Txbeamsteering, cfg.sigmarcs, geo.betark, rho, cfg.solver_opts);
end

function e = evaluate_design(cfg, geo, ch, Dsn, noise_adv, w_i, w_t, z)
% Legitimate metrics on the designed beams, attack on the actual channels.
% Fields are created in the order of the blank struct in ONE_RUN, so that
% the results can be assigned into its struct arrays.
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;
W = cat(2, w_i, reshape(w_t, cfg.Mtx, 1, cfg.Ntx));     % Mtx x (Nue+1) x Ntx
Rall = AN_Em(Dsn, sigma_ue, cfg.Ntx, cfg.Nue, cfg.Mue, ch.H(:, :, cfg.adversary), ...
    noise_adv, W, z, cfg.Mtx);
Rk = reshape(mean(Rall, 3), cfg.Mtx, cfg.Mtx, cfg.Ntx);
[success, err, est, dom] = localize_target(cfg, geo, Rk);
mh = localize_multi(cfg, geo, Rk);

m = sensing_metrics(cfg, geo, w_i, w_t, z);
e.det = double(success);
e.err = err;
e.eta = m.eta;
e.snr = 10*log10(m.gamma_t);
e.snrz = 10*log10(m.gamma_t_z);
% Dominance diagnostics over the APs with resolvable true and fake directions.
res = geo.resolvable;
e.fdom = mean(dom.fake_dominant(res));
e.pkf = mean(dom.peak_class(res) == 1);
e.pkt = mean(dom.peak_class(res) == 2);
e.mhv = double(mh.det_votes);
e.mhr = double(mh.det_rand);
e.mha = double(mh.det_aware);
e.mhs = double(mh.det_set);
fcell = floor((cfg.PosFake - cfg.area([1 3])) / cfg.cell_size);
e.fhit = double(isequal(floor((est - cfg.area([1 3])) / cfg.cell_size), fcell));
e.est = est;
e.pd = detection_pd(cfg, m.V, cfg.pd_offset_dB, cfg.pd_N_list);
end
