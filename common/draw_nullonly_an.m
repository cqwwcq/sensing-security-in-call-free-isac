function z = draw_nullonly_an(cfg, Hhat, u, w_i, w_t, p_an)
%DRAW_NULLONLY_AN Gaussian AN projected onto the effective-channel null space.
%
%   z = DRAW_NULLONLY_AN(cfg, Hhat, u, w_i, w_t, p_an) draws a complex
%   Gaussian vector, projects the stacked vector onto N(G_eff) of (6), where
%   row i of G_eff is u_i^H [Hhat_{i,1}, ..., Hhat_{i,Ntx}], and scales it
%   with one common factor so that the null is preserved.  The factor aims
%   at the per-AP AN powers p_an (1 x Ntx, W) in total, but is reduced if an
%   AP would exceed its budget P_T left after w_i and w_t.  z is Mtx x Ntx.

Ntx = cfg.Ntx; Mtx = cfg.Mtx; Nue = cfg.Nue;
P_T_lin = 10^(cfg.P_T/10) / 1000;

G = zeros(Nue, Ntx*Mtx);
for i = 1:Nue
    G(i, :) = u{i}' * Hhat(:, :, i);
end
r = (randn(Ntx*Mtx, 1) + 1j*randn(Ntx*Mtx, 1)) / sqrt(2);
zbar = r - G' * (pinv(G*G') * (G*r));
zk = reshape(zbar, Mtx, Ntx);
pk = sum(abs(zk).^2, 1);

p_beam = zeros(1, Ntx);
for k = 1:Ntx
    p_beam(k) = norm(w_i(:, :, k), 'fro')^2 + norm(w_t(:, k))^2;
end
p_left = max(P_T_lin - p_beam, 0);

s2 = min(sum(p_an) / sum(pk), min(p_left ./ pk));
z = sqrt(s2) * zk;
end
