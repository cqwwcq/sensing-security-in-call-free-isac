function geo = setup_geometry(cfg)
%SETUP_GEOMETRY Angles, steering vectors and large-scale gains of Section II.
%
%   geo.phi, geo.phi_fake      true and fake directions phi_k, phi_fake,k
%   geo.Txbeamsteering         a(phi_k) as a 1 x Ntx x Mtx array, the layout
%   geo.Fake_Txbeamsteering    expected by AN_run_optimization_*.m
%   geo.betark                 Ntx x Nrx two-way path gains g_{r,k}
%   geo.beta_ue                Nue x Ntx large-scale fading beta_{i,k}

[phi, d_tk] = compute_angle_dist(cfg.PosTx, cfg.PosTarget);   % Ntx x 1
phi_fake = compute_angle_dist(cfg.PosTx, cfg.PosFake);
[~, d_tr] = compute_angle_dist(cfg.PosRx, cfg.PosTarget);     % Nrx x 1
[~, d_ik] = compute_angle_dist(cfg.PosUE, cfg.PosTx);         % Nue x Ntx

geo.phi = phi(:).';
geo.phi_fake = phi_fake(:).';
geo.Txbeamsteering = beamsteering(geo.phi, cfg.Mtx);
geo.Fake_Txbeamsteering = beamsteering(geo.phi_fake, cfg.Mtx);
geo.a_true = reshape(geo.Txbeamsteering, cfg.Ntx, cfg.Mtx).';      % Mtx x Ntx
geo.a_fake = reshape(geo.Fake_Txbeamsteering, cfg.Ntx, cfg.Mtx).';

% g_{r,k} = lambda^2 / ((4 pi)^3 d_{t,k}^2 d_{t,r}^2), stored as betark(k, r)
geo.betark = cfg.lambda_c^2 ./ ((4*pi)^3 * (d_tk(:).^2) * (d_tr(:).'.^2));
geo.beta_ue = d_ik .^ (-cfg.pathloss_exp);
geo.d_tk = d_tk(:).';
geo.d_tr = d_tr(:).';
end
