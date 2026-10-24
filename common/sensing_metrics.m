function m = sensing_metrics(cfg, geo, w_i, w_t, z)
%SENSING_METRICS Sensing SNR, AN-echo term, fake-peak terms and eta.
%
%   m.gamma_t     sensing SNR (2), linear, AN echo discarded
%   m.gamma_t_z   sensing SNR when the known AN echo is also exploited,
%                 gamma_t + sum_k c_k T_k (Reviewer 2, Comment 1)
%   m.S, m.T, m.F per-AP |a^H(phi_k) W_k|^2, |a^H(phi_k) z_k|^2 and
%                 |a^H(phi_fake,k) z_k|^2 of Section III-B
%   m.eta         realized AN-power ratio
%   m.V           Nrx x Ntx x (Nue+1) echo vectors sqrt(sigma_rcs^2 g_{r,k})
%                 a^H(phi_k) W_k, used by DETECTION_PD

Ntx = cfg.Ntx; Nrx = cfg.Nrx; Nue = cfg.Nue;
sigma_rx = 10^(cfg.sigma_rx_dBm/10) / 1000;
c = sum(cfg.sigmarcs^2 * geo.betark, 2).' / Nrx;     % 1 x Ntx, as in the solver

m.S = zeros(1, Ntx); m.T = zeros(1, Ntx); m.F = zeros(1, Ntx);
m.V = zeros(Nrx, Ntx, Nue + 1);
p_tot = 0;
for k = 1:Ntx
    Wk = [w_i(:, :, k), w_t(:, k)];
    aW = geo.a_true(:, k)' * Wk;                      % 1 x (Nue+1)
    m.S(k) = norm(aW)^2;
    m.T(k) = abs(geo.a_true(:, k)' * z(:, k))^2;
    m.F(k) = abs(geo.a_fake(:, k)' * z(:, k))^2;
    for r = 1:Nrx
        m.V(r, k, :) = sqrt(cfg.sigmarcs^2 * geo.betark(k, r)) * aW;
    end
    p_tot = p_tot + norm(Wk, 'fro')^2 + norm(z(:, k))^2;
end
m.gamma_t = sum(c .* m.S) / sigma_rx;
m.gamma_t_z = sum(c .* (m.S + m.T)) / sigma_rx;
m.eta = sum(abs(z(:)).^2) / p_tot;
end
