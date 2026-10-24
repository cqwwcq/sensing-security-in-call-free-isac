function run_outage_sweep(out_file)
%RUN_OUTAGE_SWEEP Sweep behind Fig. 4: communication outage vs CSI NMSE.
%
%   For each NMSE and channel realization, the beams are designed on an
%   estimate with link NMSE nu (GEN_CHANNELS 'nmse' mode) and the outage
%   Pr{min_i gamma_i^AN < Gamma_comm} is evaluated over cfg.num_err draws of
%   the actual channel H = Hhat + E with the designed combiners.  Designs:
%     noan   no-AN baseline
%     null   proposed nominal beams with the AN redrawn as Gaussian AN in
%            the effective-channel null space at the same power (this
%            reproduces the published m_snu = m_sn)
%     nom    proposed nominal design (rho = cfg.rou)
%     rob    proposed robust design (20)
%     rnoan  robust no-AN baseline (Reviewer 1, Comment 5), if
%            cfg.include_robust_noan is true
%   Results go to outage_results.mat in this folder, with the variables of
%   the published file (m_*, c_*, m_sn*, n_kept, n_drop, nmse, cfg).
%
%   Since ||u_i|| = 1 and E_{i,k} has i.i.d. CN(0, nmse*beta_ik) entries,
%   u_i^H E_{i,k} has i.i.d. CN(0, nmse*beta_ik) entries, so the error is
%   drawn directly on the effective channels.

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'common'));
if nargin < 1
    out_file = fullfile(here, 'outage_results.mat');
end

cfg = sim_config();
cfg.rou = 1.3;
cfg.nmse_db = -28:2:-4;
cfg.num_runs = 300;
cfg.num_err = 3000;
cfg.base_seed = 20260824;
cfg.include_robust_noan = true;
geo = setup_geometry(cfg);

designs = {'noan', 'null', 'nom', 'rob'};
if cfg.include_robust_noan
    designs{end+1} = 'rnoan';
end
nd = numel(designs);
nm = numel(cfg.nmse_db);
R = cfg.num_runs;

out = zeros(R, nm, nd);       % per-run outage probability
snr = zeros(R, nm, nd);       % sensing SNR (dB)
okr = false(R, nm);
t0 = tic;
if cfg.use_parfor
    parfor r = 1:R
        [out(r, :, :), snr(r, :, :), okr(r, :)] = one_run(cfg, geo, r, designs);
    end
else
    for r = 1:R
        [out(r, :, :), snr(r, :, :), okr(r, :)] = one_run(cfg, geo, r, designs);
    end
end
cfg.elapsed = toc(t0);

S.cfg = cfg;
S.nmse = cfg.nmse_db;
S.n_kept = sum(okr, 1);
S.n_drop = R - S.n_kept;
snr_names = struct('noan', 'm_sna', 'null', 'm_snu', 'nom', 'm_sn', 'rob', 'm_sr', 'rnoan', 'm_srn');
out_names = struct('noan', 'noan', 'null', 'null', 'nom', 'nom', 'rob', 'rob', 'rnoan', 'rnoan');
for d = 1:nd
    m = zeros(1, nm); c = zeros(1, nm); s = zeros(1, nm);
    for j = 1:nm
        x = out(okr(:, j), j, d);
        m(j) = mean(x);
        c(j) = 1.96 * std(x) / sqrt(numel(x));
        s(j) = mean(snr(okr(:, j), j, d));
    end
    S.(['m_' out_names.(designs{d})]) = m;
    S.(['c_' out_names.(designs{d})]) = c;
    S.(snr_names.(designs{d})) = s;
end
save(out_file, '-struct', 'S');
fprintf('Saved %s (%.0f s).\n', out_file, cfg.elapsed);
end

function [out, snr, ok] = one_run(cfg, geo, r, designs)
nm = numel(cfg.nmse_db);
nd = numel(designs);
out = nan(nm, nd);
snr = nan(nm, nd);
ok = false(1, nm);
for j = 1:nm
    % Same seed for every NMSE: the estimate differs only in its scaling.
    rng(cfg.base_seed + r);
    nu = 10^(cfg.nmse_db(j)/10);
    ch = gen_channels(cfg, geo, 'nmse', nu);
    E = draw_effective_errors(cfg, ch.beta_err);
    okj = true;
    nom = [];
    for d = 1:nd
        switch designs{d}
            case 'noan'
                [w_i, w_t, u, z, okd] = solve(cfg, geo, ch, struct('an_mode', 'none'));
            case 'nom'
                [w_i, w_t, u, z, okd] = solve(cfg, geo, ch, struct('an_mode', 'fake', ...
                    'sym_null', cfg.sym_null_prop, 'fake_tleak', cfg.fake_tleak, ...
                    'fake_mask', fake_mask(cfg, geo)));
                nom = struct('w_i', w_i, 'w_t', w_t, 'u', {u}, 'z', z);
            case 'null'
                continue;                 % filled in after 'nom'
            case 'rob'
                [w_i, w_t, u, z, okd] = solve(cfg, geo, ch, struct('an_mode', 'fake', ...
                    'sym_null', cfg.sym_null_prop, 'fake_tleak', cfg.fake_tleak, ...
                    'fake_mask', fake_mask(cfg, geo), 'beta_error', ch.beta_err));
            case 'rnoan'
                [w_i, w_t, u, z, okd] = solve(cfg, geo, ch, ...
                    struct('an_mode', 'none', 'beta_error', ch.beta_err));
        end
        okj = okj && okd;
        [out(j, d), snr(j, d)] = evaluate(cfg, geo, ch, E, w_i, w_t, u, z);
    end
    d = find(strcmp(designs, 'null'));
    if ~isempty(d) && ~isempty(nom)
        p_an = sum(abs(nom.z).^2, 1);
        z = draw_nullonly_an(cfg, ch.Hhat, nom.u, nom.w_i, nom.w_t, p_an);
        [out(j, d), snr(j, d)] = evaluate(cfg, geo, ch, E, nom.w_i, nom.w_t, nom.u, z);
    end
    ok(j) = okj;
end
% Shape the outputs as the (r, :, :) slices of the parfor arrays.
out = reshape(out, 1, nm, nd);
snr = reshape(snr, 1, nm, nd);
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

function E = draw_effective_errors(cfg, beta_err)
% num_err x (Ntx*Mtx) x Nue effective-channel errors u_i^H E_{i,k}.
Ntx = cfg.Ntx; Mtx = cfg.Mtx; Nue = cfg.Nue; D = cfg.num_err;
E = (randn(D, Ntx*Mtx, Nue) + 1j*randn(D, Ntx*Mtx, Nue)) / sqrt(2);
for i = 1:Nue
    scale = kron(sqrt(beta_err(:, i)).', ones(1, Mtx));    % 1 x Ntx*Mtx
    E(:, :, i) = E(:, :, i) .* scale;
end
end

function [p_out, snr_dB] = evaluate(cfg, geo, ch, E, w_i, w_t, u, z)
% Outage over the error draws with fixed combiners and beams.
Nue = cfg.Nue; Ntx = cfg.Ntx; Mtx = cfg.Mtx;
Gamma = 10^(cfg.Gamma_comm/10);
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;
B = [reshape(permute(w_i, [1 3 2]), Ntx*Mtx, Nue), w_t(:), z(:)];   % stacked beams
fail = false(cfg.num_err, 1);
for i = 1:Nue
    g = u{i}' * ch.Hhat(:, :, i) + E(:, :, i);       % num_err x Ntx*Mtx
    P = abs(g * B).^2;                               % num_err x (Nue+2)
    sinr = P(:, i) ./ (sum(P, 2) - P(:, i) + sigma_ue);
    fail = fail | (sinr < Gamma);
end
p_out = mean(fail);
m = sensing_metrics(cfg, geo, w_i, w_t, z);
snr_dB = 10*log10(m.gamma_t);
end
