function sinr = user_sinr(cfg, H, u, w_i, w_t, z)
%USER_SINR Post-combining SINR (5) of every user, linear, 1 x Nue.
%
%   H is Mue x (Ntx*Mtx) x Nue in the solver layout (designed-for or actual
%   channels), u the cell of unit-norm combiners, and w_i, w_t, z the beams.
%   The interference includes the other users' streams, the sensing stream
%   and the AN.

Nue = cfg.Nue; Ntx = cfg.Ntx; Mtx = cfg.Mtx;
sigma_ue = 10^(cfg.sigma_ue_dBm/10) / 1000;
B = [reshape(permute(w_i, [1 3 2]), Ntx*Mtx, Nue), w_t(:), z(:)];
sinr = zeros(1, Nue);
for i = 1:Nue
    P = abs(u{i}' * H(:, :, i) * B).^2;
    sinr(i) = P(i) / (sum(P) - P(i) + sigma_ue);
end
end
