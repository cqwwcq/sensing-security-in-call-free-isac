function ch = gen_channels(cfg, geo, mode, nmse)
%GEN_CHANNELS AP-user channels and their estimates in the solver layout.
%
%   ch = GEN_CHANNELS(cfg, geo, 'pilot') draws H_{i,k} = sqrt(beta_ik) G_ik
%   and forms the LMMSE estimate from tau_p-length orthogonal pilots of
%   per-antenna power P_p / M_UE (Section II), so that the error variance is
%   epsilon_ik^2 * beta_ik with epsilon_ik^2 = sigma_p^2 / (tau_p p beta_ik
%   + sigma_p^2).
%
%   ch = GEN_CHANNELS(cfg, geo, 'nmse', nmse) draws the estimate directly
%   with a common link NMSE (linear), as in the Fig. 4 sweep: Hhat has
%   variance (1 - nmse) beta_ik and the error E variance nmse * beta_ik.
%   The actual channel is not drawn in this mode; the caller adds error
%   draws with ch.beta_err.
%
%   Channels are stored as Mue x (Ntx*Mtx) x Nue, i.e. H(:, cols_k, i) is
%   the Mue x Mtx channel from AP k to user i.  ch.beta_err is Ntx x Nue,
%   the layout of opts.beta_error in AN_run_optimization_robust.m.

Ntx = cfg.Ntx; Nue = cfg.Nue; Mtx = cfg.Mtx; Mue = cfg.Mue;
ch.H = zeros(Mue, Ntx*Mtx, Nue);
ch.Hhat = zeros(Mue, Ntx*Mtx, Nue);
ch.beta_err = zeros(Ntx, Nue);

switch mode
    case 'pilot'
        p = 10^(cfg.P_p_dBm/10) / 1000 / Mue;
        s2 = 10^(cfg.sigma_p_dBm/10) / 1000;
        for i = 1:Nue
            for k = 1:Ntx
                cols = (k-1)*Mtx + (1:Mtx);
                b = geo.beta_ue(i, k);
                H = sqrt(b) * crandn(Mue, Mtx);
                Y = H + sqrt(s2 / (cfg.tau_p * p)) * crandn(Mue, Mtx);
                ch.H(:, cols, i) = H;
                ch.Hhat(:, cols, i) = cfg.tau_p*p*b / (cfg.tau_p*p*b + s2) * Y;
                ch.beta_err(k, i) = b * s2 / (cfg.tau_p*p*b + s2);
            end
        end
    case 'nmse'
        for i = 1:Nue
            for k = 1:Ntx
                cols = (k-1)*Mtx + (1:Mtx);
                b = geo.beta_ue(i, k);
                ch.Hhat(:, cols, i) = sqrt((1 - nmse) * b) * crandn(Mue, Mtx);
                ch.beta_err(k, i) = nmse * b;
            end
        end
        ch.H = [];
    otherwise
        error('Unknown channel mode "%s".', mode);
end
end

function X = crandn(m, n)
X = (randn(m, n) + 1j*randn(m, n)) / sqrt(2);
end
