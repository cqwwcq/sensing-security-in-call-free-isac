function cfg = sim_config()
%SIM_CONFIG Default simulation parameters shared by all sweep drivers.
%
%   Values marked [Table I] are taken from Table I of paper_draft.tex.
%   Values marked [ASSUMED] are not stated in the paper and were filled in
%   from the earlier multi-static ISAC paper in this repository (path-loss
%   exponent, RCS, UE placement, AP ordering, search area) or chosen here.
%   Replace them with the values of the original driver scripts to
%   reproduce the published numbers exactly.

% --- Network size and antennas [Table I] ---
cfg.Ntx = 8;
cfg.Nrx = 4;
cfg.Nue = 4;
cfg.Mtx = 16;
cfg.Mue = 16;
cfg.Mrx = 1;          % [ASSUMED] the solver's c_k has no receive array gain

% --- Power, QoS, noise [Table I] ---
cfg.P_T = 30;                 % dBm per transmit AP
cfg.Gamma_comm = 10;          % dB
cfg.sigma_ue_dBm = -94;       % user noise power
cfg.sigma_rx_dBm = -94;       % sensing-receiver noise power
cfg.sigma_p_dBm = -94;        % pilot noise power
cfg.epsilon = 1e-10;          % null-constraint tolerance used in the sweeps

% --- Carrier and pilots [Table I] ---
cfg.fc = 1.9e9;
cfg.lambda_c = 3e8 / cfg.fc;
cfg.tau_p = 64;
cfg.P_p_dBm = 20;             % total pilot power per user

% --- Geometry [Table I], ordering of the earlier paper [ASSUMED] ---
cfg.PosTx = [-200 -200; 200 200; -200 200; 200 -200; ...
             0 -200; 0 200; -200 0; 200 0];
cfg.PosRx = [150 150; -150 -150; -150 150; 150 -150];
cfg.PosTarget = [-75 75];
cfg.PosFake = [5 -25];
cfg.PosUE = [120 120; -120 -120; -120 120; 120 -120];   % [ASSUMED]
cfg.adversary = 1;            % AN_Em.m treats user 1 as the adversary

% --- Large-scale fading and target [ASSUMED, earlier paper] ---
cfg.pathloss_exp = 3;         % beta_ik = d_ik^(-pathloss_exp)
cfg.sigmarcs = sqrt(10);      % RCS variance 10 dBsm (solver uses sigmarcs^2)

% --- Reconstruction attack and localization ---
cfg.N_snapshot = 128;         % hard-coded in AN_Em.m [Table I]
cfg.cell_size = 10;           % m [Table I]
cfg.n_ang = 3604;             % angle-scan points over [-pi, pi)
cfg.area = [-200 200 -200 200];  % search area [xmin xmax ymin ymax] [ASSUMED]
cfg.mirror_rule = 'random';   % ULA front-back ambiguity: 'random' or 'both'
cfg.mh_exclusion = 3;         % multi-hypothesis attack: second hypothesis lies
                              % more than this many cells from the first

% --- Design sweep ---
cfg.rho_list = [0.05 0.15 0.3 0.5 0.75 1.05 1.5 2.2 3.5];
cfg.sym_null_prop = true;     % enforce (9f) in the proposed design, as the
                              % paper states; the published sweep used false
cfg.fake_tleak = true;        % revised (9g) with T_k on the right-hand side
                              % (Reviewer 1, Comment 2); published: false
cfg.fake_resolvable_only = true;  % impose (9g) only at APs with resolvable
                              % true and fake directions; published: false
cfg.fake_upper = Inf;          % cap kappa on the fake peak; e.g. 1.2 with
                              % rho ~ 1 gives power-matched peaks
cfg.solver_opts = struct('max_iter', 30, 'tol', 4e-3);

% --- Detection probability (Reviewer 3, Comment 2) ---
cfg.Pfa = 1e-6;
cfg.pd_N_list = [5 128];          % integration lengths (N >= Nue+1)
cfg.pd_offset_dB = -30:2:0;       % SNR offsets for Pd-vs-RCS curves
cfg.pd_alpha_draws = 2000;        % Swerling-I echo draws per realization

% --- Monte Carlo ---
cfg.use_parfor = true;
end
