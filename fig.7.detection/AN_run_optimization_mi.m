function [gamma_t_opt, w_i_opt, w_t_opt,u_opt,z_opt, solve_time, success,power_z,miu_z,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense,ITET_NUM,obj_history] = AN_run_optimization_mi(...
    Ntx, Nue, Mtx, Nrx,Mue,P_T_dBm, Gamma_comm_dB, sigma_noise_dBm, sigma_ap_dBm,epsilon, ...
    H_all, Txbeamsteering, Fake_Txbeamsteering,sigmarcs, betark,rou,varargin)
    
    tic; % Start the solver timer
    
    % --- 1. Convert parameters to linear scale ---
    P_T_lin = 10^(P_T_dBm / 10) / 1000;      % W
    Gamma_comm_lin = 10^(Gamma_comm_dB / 10);
    sigma_noise_lin = 10^(sigma_noise_dBm / 10) / 1000; % 
    sigma_ap_linear = 10^(sigma_ap_dBm / 10) / 1000; % 
    epsilon_lin = epsilon; % Assume epsilon is already in linear scale

    MAX_ITER_MAIN = 30;
    tol_main = 4e-3;
    SYM_NULL = false;   % eq:sym-null on the sensing beam (Section III-B)
    % Revised prob:fake: add the AN response T_k in the true direction to
    % the right-hand side, so that the fake-direction response dominates
    % the true-direction one for any rho > 1 (Reviewer 1, Comment 2).
    FAKE_TLEAK = false;
    % APs at which prob:fake is imposed; by default all.  Restricting it to
    % APs whose true and fake directions are resolvable avoids spending
    % sensing power where no distinguishable fake peak can be formed.
    FAKE_MASK = true(1, Ntx);
    % an_mode selects the design being solved:
    %   'fake' : proposed design, prob:null + prob:null-t + prob:fake
    %   'none' : no-AN baseline, sensing objective under QoS and power only
    AN_MODE = 'fake';
    if ~isempty(varargin)
        opts = varargin{1};
        if isfield(opts, 'max_iter')
            MAX_ITER_MAIN = opts.max_iter;
        end
        if isfield(opts, 'tol')
            tol_main = opts.tol;
        end
        if isfield(opts, 'sym_null')
            SYM_NULL = opts.sym_null;
        end
        if isfield(opts, 'fake_tleak')
            FAKE_TLEAK = opts.fake_tleak;
        end
        if isfield(opts, 'fake_mask')
            FAKE_MASK = logical(opts.fake_mask);
        end
        if isfield(opts, 'an_mode')
            AN_MODE = opts.an_mode;
        end
    end
    switch AN_MODE
        case 'fake'
            USE_AN = true;  USE_FAKE = true;
        case 'none'
            USE_AN = false; USE_FAKE = false; SYM_NULL = false;
        case 'nullonly'
            % Null-space-only (Gaussian-AN) baseline: prob:null and
            % prob:null-t are enforced but prob:fake is omitted.  Since the
            % sensing objective does not involve z and the joint null makes
            % its user interference vanish, the AN direction is immaterial
            % to the legitimate metrics.  The transmit beams are therefore
            % optimized here over the reduced budget left by the prescribed
            % AN power, and the AN itself is drawn afterwards by projecting
            % a random vector onto the null space of the returned combiners.
            USE_AN = false; USE_FAKE = false; SYM_NULL = true;
        otherwise
            error('Unknown an_mode "%s".', AN_MODE);
    end
    % Receiver constraint set consistent with the nulls actually imposed.
    if strcmp(AN_MODE, 'nullonly')
        RX_MODE = 2;              % Q_i = [q_iz, q_it], with q_iz = 0 here
    elseif ~USE_AN
        RX_MODE = 0;              % unconstrained MMSE
    elseif SYM_NULL
        RX_MODE = 2;              % Q_i = [q_iz, q_it]
    else
        RX_MODE = 1;              % Q_i = q_iz
    end
    obj_history = nan(1, MAX_ITER_MAIN + 1);


    % --- 2. Precompute all channels and matrices ---
    % This maps the simulation parameters to the optimization formulation.
    
    h_ik_row_list = cell(Nue, Ntx);     % (i, k) stores the 1xMtx row vector h_ik
    a_phi_k_H_list = cell(Ntx, 1);      % (k) stores the 1xMtx row vector a_phi_k^H
    a_phi_k_H_list_fake = cell(Ntx,1);  % (k) stores the 1xMtx row vector a_phi_k_fake^H
    a_theta_ik_H_list = cell(Nue, Ntx); % (i, k) stores the 1xMtx row vector a_theta_ik^H
    for k = 1:Ntx
        % 2a. Sensing-target vector a_phi_k^H
        a_phi_k = squeeze(Txbeamsteering(:,k,:)); % Mtx x 1 row vector
        a_phi_k_H_list{k} = a_phi_k'; % 1 x Mtx column vector
        a_phi_k_fake = squeeze(Fake_Txbeamsteering(:,k,:));
        a_phi_k_H_list_fake{k} = a_phi_k_fake'; % 1 x Mtx column vector
        for i = 1:Nue
            % 2b. Communication channel h_ik and Q_ik
            % Extract the channel from AP k to user i.
            H_list_raw{i, k} = H_all(:, (1+(k-1)*Mtx):k*Mtx, i); % (MuexMtx)    
          
            
            
        end
    end
    
     
    % 2d. Sensing-objective weights c_k
    c_k_vec = zeros(1, Ntx);
    for k = 1:Ntx
    % c_k_vec(k) = (1 / (Nrx * sigma_ap_linear)) * sum(sigmarcs^2 * betark(k, :));
        c_k_vec(k) = (1 / (Nrx )) * sum(sigmarcs^2 * betark(k, :));

    end
 

    % --- 3. Stage 1: Find a feasible point (P-Feas-CCCP) ---
    % fprintf('    [Stage 1] Searching for a feasible point...\n');
    % p_init_total = P_T_lin * 0.001; 
    % p_init_per_beam = p_init_total / (Ntx * (Nue + 2));% Initial power per beam, including AN
    % scale = sqrt(p_init_per_beam); % Assumed power level of each beam
    % (Mtx x Nue x Ntx)
    % w_i_r = (randn(Mtx, Nue, Ntx) + 1j*randn(Mtx, Nue, Ntx)) * scale / sqrt(Mtx);
    % (Mtx x Ntx)
    % w_t_r = (randn(Mtx, Ntx) + 1j*randn(Mtx, Ntx)) * scale / sqrt(Mtx);
     p_wi_per_link=P_T_lin * 0.6/Nue*0.8;
     p_kickstart = P_T_lin * 0.3*0.8;
     P_an=P_T_lin * 0.1*0.8;
    if ~USE_AN
        % The no-AN baseline starts from the same beams with the AN budget
        % reassigned to the sensing beam, so that the initial per-AP power
        % matches the AN designs.
        p_kickstart = p_kickstart + P_an;
        P_an = 0;
    end
    for k = 1:Ntx
         for i = 1:Nue
            H_mat  = H_list_raw{i, k};       
            % Use conjugate transpose to obtain the beam direction (Mtx x 1).
            % If h is a row vector, h' is a column vector.
            [~, ~, V] = svd(H_mat);         
            dominant_dir = V(:, 1); % Size: Mtx x 1
            % Normalize and assign.
            w_i_r(:, i, k) = dominant_dir * sqrt(p_wi_per_link);
        end
        target_vec = a_phi_k_H_list{k}'; 
        dir_norm = target_vec / norm(target_vec);
        w_t_r(:, k) = dir_norm * sqrt(p_kickstart );
        fake_target_vec = a_phi_k_H_list_fake{k}'; 
        z_r(:, k) = fake_target_vec / norm (fake_target_vec)* sqrt(P_an);
    end
    % z_r = (randn(Mtx, Ntx) + 1j*randn(Mtx, Ntx)) * scale / sqrt(Mtx);
    % u_opt = AN_update_receiver_mmse(Ntx, Nue, Mue, Mtx, H_list_raw, w_i_r, w_t_r, z_r, sigma_noise_lin);
    % 
    % h_ik_row_list = cell(Nue, Ntx);
    %     for k = 1:Ntx
    %         for i = 1:Nue
    %             % h_eff = u_i' * H_raw
    %             h_ik_row_list{i, k} = u_opt{i}' * H_list_raw{i, k};
    %         end
    %     end
    % MAX_ITER_PHASE_1 = 20;
    % found_feasible = false;
    % 
    % for r = 1:MAX_ITER_PHASE_1
    % 
    %     cvx_begin quiet
    %         % Define CVX variables
    %         variable w_i_vars(Mtx, Nue, Ntx) complex;
    %         variable w_t_vars(Mtx, Ntx) complex;
    %         variable s_var; % Slack variable
    %         variable z_vars(Mtx, Ntx) complex;
    % 
    %         minimize(s_var);
    % 
    %         subject to
    %             for k = 1:Ntx
    %                 % (C2) Total power constraint (per AP)
    % 
    %                 power_expr_k = sum(sum_square_abs(w_i_vars(:, :, k))) + sum_square_abs(w_t_vars(:, k))+sum_square_abs(z_vars(:, k));
    %                 power_expr_k <= P_T_lin;
    %                 % Null-space constraint (C3): aggregate AN lies in the null space of U
    %                 for i = 1:Nue
    %                     an_ik =h_ik_row_list{i, k} * z_vars(:, k);
    %                     abs(an_ik)<= epsilon_lin;
    %                 end
    %                % an_sense_leakage = a_phi_k_H_list{k} * z_vars(:, k);
    %                %  abs(an_sense_leakage) <= epsilon_lin; 
    %                 % (C4) Beampattern constraint
    %                 w_total_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k), z_vars(:,k)];
    %                 w_total_k_r = [squeeze(w_i_r(:,:,k)), w_t_r(:,k), z_r(:,k)];
    %                 W_ISAC_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k)];         
    %                 expression grad_fake_term
    %                 grad_fake_term =2 * real( a_phi_k_H_list_fake{k} * w_total_k_var * w_total_k_r'*a_phi_k_H_list_fake{k}');
    %                 const_fake_term = real( a_phi_k_H_list_fake{k} * w_total_k_r * w_total_k_r'*a_phi_k_H_list_fake{k}' );
    %                 f=pow_pos( norm( a_phi_k_H_list{k} * w_total_k_var, 'fro' ), 2 );
    %                 1.1*f-grad_fake_term+const_fake_term<=0;
    % 
    % 
    %             end
    %             for i = 1:Nue         
    %                 % (C1'_feas) QoS constraint
    %                 f_i_interference = 0;
    %                 for l = 1:Nue
    %                     if l ~= i
    %                         interf_sum_l = 0;
    %                         for k = 1:Ntx
    %                             interf_sum_l = interf_sum_l +h_ik_row_list{i, k} * w_i_vars(:, l, k);
    %                         end
    %                 f_i_interference=f_i_interference+ square_abs(interf_sum_l);
    %                     end
    %                 end
    %                 sense_interf_sum = 0;
    %                 for k = 1:Ntx
    %                     sense_interf_sum = sense_interf_sum + h_ik_row_list{i, k} * w_t_vars(:, k);
    %                 end
    %                 f_i_interference = f_i_interference + square_abs(sense_interf_sum);
    %                 an_interf_sum = 0;
    %                 for k = 1:Ntx
    %                     an_interf_sum = an_interf_sum + h_ik_row_list{i, k} * z_vars(:, k);
    %                 end
    %                 f_i_interference = f_i_interference + square_abs(an_interf_sum);
    % 
    %                 S_val_r = 0;
    %                 for k = 1:Ntx
    %                     S_val_r = S_val_r + h_ik_row_list{i, k} * w_i_r(:, i, k);
    %                 end
    %                 S_expr_var = 0;
    %                 for k = 1:Ntx
    %                     S_expr_var = S_expr_var + h_ik_row_list{i, k} * w_i_vars(:, i, k);
    %                 end
    %                 Signal_linear = 2 * real(conj(S_val_r) * S_expr_var) - abs(S_val_r)^2;
    %                 (Gamma_comm_lin * (f_i_interference + sigma_noise_lin) - Signal_linear) / sigma_noise_lin <= s_var;
    % 
    %             end
    % 
    %     cvx_end
    % 
    %     if ~contains(cvx_status, 'Solved')
    % 
    %         fprintf('    [Stage 1] Iteration %d: solver failed. Status: %s\n', r, cvx_status);
    %         break;
    %     end
    % 
    %     fprintf('    [Stage 1] Iteration %d: s = %.6e\n', r, s_var);
    % 
    %     if s_var <= 0
    %         fprintf('    [Stage 1] Feasible solution found.\n');
    %         found_feasible = true;
    %         w_i_r = w_i_vars;
    %         w_t_r = w_t_vars;
    %         z_r = z_vars;
    %         break;
    %     end
    %     w_i_r = w_i_vars;
    %     w_t_r = w_t_vars;
    %     z_r = z_vars;
    %     end
    % 
    % 
    % 
    % if ~found_feasible
    %     fprintf('    [Stage 1] No feasible solution found. The algorithm terminates.\n');
    %     gamma_t_opt = NaN;
    %     w_i_opt = {};
    %     w_t_opt = {};
    %     z_opt = {};
    %     solve_time = toc;
    %     success = false;
    %     return; 
    % end
    % --- 4. Stage 2: Maximize sensing SINR via alternating optimization (P-CCCP) ---


    
    initial_obj_val = 0;
    for k = 1:Ntx
        w_ISAC_k_init = [squeeze(w_i_r(:,:,k)), w_t_r(:,k)];
        initial_obj_val = initial_obj_val + c_k_vec(k) * ...
            real(a_phi_k_H_list{k} * w_ISAC_k_init * w_ISAC_k_init' * a_phi_k_H_list{k}');
    end
    obj_history(1) = initial_obj_val;
    last_obj_val = initial_obj_val;
    fprintf('    [Stage 2] Solving the main optimization problem...\n');

    for p = 1:MAX_ITER_MAIN
        u_opt = AN_update_receiver_mmse_null(Ntx, Nue, Mue, Mtx, H_list_raw, w_i_r, w_t_r, z_r, sigma_noise_lin, RX_MODE);
        % u_op=tupdate_receiver_mmse(Ntx, Nue, Mue, Mtx, H_list_raw, w_i_r, w_t_r, sigma_noise_lin);
        h_ik_row_list = cell(Nue, Ntx);
        for k = 1:Ntx
            for i = 1:Nue
                % h_eff = u_i' * H_raw
                h_ik_row_list{i, k} = u_opt{i}' * H_list_raw{i, k};
            end
        end
        cvx_begin quiet
        cvx_precision high
            variable w_i_vars(Mtx, Nue, Ntx) complex;
            variable w_t_vars(Mtx, Ntx) complex;
            variable z_vars(Mtx, Ntx) complex;
            % --- Build the objective function ---
 
            obj_linear_expr=0;
            for k = 1:Ntx
                A_k_true_matrix=a_phi_k_H_list{k}' * a_phi_k_H_list{k};
                w_total_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k), z_vars(:,k)];
                w_total_k_r = [squeeze(w_i_r(:,:,k)), w_t_r(:,k), z_r(:,k)];
                w_ISAC_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k)];
                w_ISAC_k_r = [squeeze(w_i_r(:,:,k)), w_t_r(:,k)];
                reward_k = c_k_vec(k) * (2 * real(trace( (A_k_true_matrix * w_ISAC_k_r)' * w_ISAC_k_var ))-pow_pos( norm( a_phi_k_H_list{k} * w_ISAC_k_r, 'fro' ), 2 ));                
                % reward_k = c_k_vec(k) * (2 * real(trace( (A_k_true_matrix * w_total_k_r)' * w_total_k_var ))-pow_pos( norm( a_phi_k_H_list{k} * w_total_k_r, 'fro' ), 2 ));                

            obj_linear_expr = obj_linear_expr + reward_k;
            end
        
            maximize( obj_linear_expr );

            subject to
                for k = 1:Ntx
                    % (C2) Total power constraint (per AP)
                   
                    power_expr_k = sum(sum_square_abs(w_i_vars(:, :, k))) + sum_square_abs(w_t_vars(:, k))+ ...
                               sum_square_abs(z_vars(:, k));
                    power_expr_k <= P_T_lin;

  
                    % (C4) Beampattern constraint
                    w_total_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k), z_vars(:,k)];
                    w_total_k_r = [squeeze(w_i_r(:,:,k)), w_t_r(:,k), z_r(:,k)];
                    W_ISAC_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k)];
                    if USE_FAKE && FAKE_MASK(k)
                        expression grad_fake_term
                        grad_fake_term =2 * real( a_phi_k_H_list_fake{k} *  z_vars(:,k) * z_r(:,k)'*a_phi_k_H_list_fake{k}');
                        const_fake_term = real( a_phi_k_H_list_fake{k} * z_r(:,k) * z_r(:,k)'*a_phi_k_H_list_fake{k}' );
                        f=pow_pos( norm( a_phi_k_H_list{k} * W_ISAC_k_var, 'fro' ), 2 );
                        if FAKE_TLEAK
                            t_leak = pow_pos(abs(a_phi_k_H_list{k} * z_vars(:,k)), 2);
                            rou*f+t_leak-grad_fake_term+const_fake_term<=0;
                        else
                            rou*f-grad_fake_term+const_fake_term<=0;
                        end
                    end
                    if ~USE_AN
                        % No-AN baseline: the AN beam is removed entirely.
                        z_vars(:,k) == 0;
                    end


                end
                
                for i = 1:Nue
                    % (C3) Null-space constraint
                    if USE_AN
                        an_joint_i = 0;
                        for k = 1:Ntx
                            an_joint_i = an_joint_i + ...
                                h_ik_row_list{i,k} * z_vars(:,k);
                        end
                        abs(an_joint_i) <= epsilon_lin;
                    end
                    % (C3b) Symmetric null: place the sensing beam under the
                    % same joint effective-channel condition as the AN.
                    if SYM_NULL
                        sense_joint_i = 0;
                        for k = 1:Ntx
                            sense_joint_i = sense_joint_i + ...
                                h_ik_row_list{i,k} * w_t_vars(:,k);
                        end
                        abs(sense_joint_i) <= epsilon_lin;
                    end
                    % (C1'_feas) QoS constraint
                    f_i_interference = 0;
                    for l = 1:Nue
                        if l ~= i
                            interf_sum_l = 0;
                            for k = 1:Ntx
                                interf_sum_l = interf_sum_l +h_ik_row_list{i, k} * w_i_vars(:, l, k);
                            end
                            f_i_interference=f_i_interference+ square_abs(interf_sum_l);
                        end
                    end
                    sense_interf_sum = 0;
                    for k = 1:Ntx
                        sense_interf_sum = sense_interf_sum + h_ik_row_list{i, k} * w_t_vars(:, k);
                    end
                    f_i_interference = f_i_interference + square_abs(sense_interf_sum);
                    an_interf_sum = 0;
                    for k = 1:Ntx
                        an_interf_sum = an_interf_sum + h_ik_row_list{i, k} * z_vars(:, k);
                    end
                    f_i_interference = f_i_interference + square_abs(an_interf_sum);
                    S_val_r = 0;
                    for k = 1:Ntx
                        S_val_r = S_val_r + h_ik_row_list{i, k} * w_i_r(:, i, k);
                    end
                    S_expr_var = 0;
                    for k = 1:Ntx
                        S_expr_var = S_expr_var + h_ik_row_list{i, k} * w_i_vars(:, i, k);
                    end
                    Signal_linear = 2 * real(conj(S_val_r) * S_expr_var) - abs(S_val_r)^2;
                    (Gamma_comm_lin * (f_i_interference + sigma_noise_lin) - Signal_linear) / sigma_noise_lin <= 0;
                end
               
        cvx_end

        
        


        current_obj_val = 0;
        if ~contains(cvx_status, 'Solved')
            
                fprintf('    [Stage 2] Iteration %d: solver failed. Status: %s\n', p, cvx_status);
                success = false;   
                break;
          
        end
   
        
        % --- Check convergence ---
        current_obj_val = 0;
        for k = 1:Ntx
            w_ISAC_k_current = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k)];
            current_obj_val = current_obj_val + c_k_vec(k) * ...
                real(a_phi_k_H_list{k} * w_ISAC_k_current * w_ISAC_k_current' * a_phi_k_H_list{k}');
        end
        obj_history(p + 1) = current_obj_val;
      


        
        fprintf('    [Stage 2] Iteration %d: sensing SNR = %.8f dB\n', p, 10*log10(current_obj_val/ sigma_ap_linear));
        
        if abs(current_obj_val - last_obj_val) / (1e-9 + abs(last_obj_val)) < tol_main
            fprintf('    [Stage 2] Algorithm converged.\n');
            w_i_r = w_i_vars; 
            w_t_r = w_t_vars; 
            z_r = z_vars;
            last_obj_val = current_obj_val;
            success = true;
            break;
        end
        
        last_obj_val = current_obj_val;
        w_i_r = w_i_vars;
        w_t_r = w_t_vars;
        z_r = z_vars;
        
        if p == MAX_ITER_MAIN
            fprintf('    [Stage 2] Maximum number of iterations reached.\n');
        end
        success = true;

    end
 
    ITET_NUM=p;   
    obj_history = obj_history(~isnan(obj_history));

  

    % --- 5. Final results ---
    gamma_t_opt = last_obj_val/ sigma_ap_linear;
    w_i_opt = w_i_r;
    w_t_opt = w_t_r;
    z_opt = z_r;
    solve_time = toc; % Stop the timer
   
      % --- 6. Final constraint verification (detailed report) ---
    fprintf('    [Check] Verifying the constraints of the final solution...\n');
    

    
    % Initialize worst-case tracking variables.
    min_sinr_dB = Inf;
    max_leakage_lin = -Inf;
    max_power_lin = -Inf;
    min_beampattern_difference= Inf;
    interf_UE=0
    interf_sense=0;
    % Check (C1) communication SINR and (C3) sensing leakage.
    for i = 1:Nue
        % Verify SINR using the optimized U.
        u_vec = u_opt{i};
        sig_sum = 0;
        for k = 1:Ntx
            sig_sum = sig_sum + H_list_raw{i, k} * w_i_opt(:, i, k);
        end
        signal_p = abs(u_vec' * sig_sum)^2;
        
        interf_p = 0;
        for l = 1:Nue
            if l ~= i
                inf_sum_l = 0;
                for k = 1:Ntx
                    inf_sum_l = inf_sum_l + H_list_raw{i, k} * w_i_opt(:, l, k);
                end
                interf_p = interf_p + abs(u_vec' * inf_sum_l)^2;
            end
        end
        sense_sum = 0;
        for k = 1:Ntx
            sense_sum = sense_sum + H_list_raw{i, k} * w_t_opt(:, k);
        end
        interf_UE=interf_UE+interf_p;
        interf_sense=interf_sense+abs(u_vec' * sense_sum)^2;
        interf_p = interf_p + abs(u_vec' * sense_sum)^2;
        an_sum = 0;
        for k = 1:Ntx
            an_sum = an_sum + H_list_raw{i, k} * z_opt(:, k);
        end
        an_interference = abs(u_vec' * an_sum)^2;

        final_sinr_lin = signal_p / (interf_p + sigma_noise_lin+an_interference);
        final_sinr_dB = 10*log10(final_sinr_lin);
        min_sinr_dB = min(min_sinr_dB, final_sinr_dB);
    end
    interf_UE=interf_UE/Nue;    
    interf_sense=interf_sense/Nue;

    % % (C3) Null-space leakage
    % % Track the maximum leakage
    for i = 1:Nue
        joint_leakage_i = 0;
    
        for k = 1:Ntx
            h_ik_row = u_opt{i}' * H_list_raw{i,k};
    
            joint_leakage_i = joint_leakage_i + ...
                h_ik_row * z_opt(:,k);
        end
    
        max_leakage_lin = max( ...
            max_leakage_lin, abs(joint_leakage_i));
    end

    for k = 1:Ntx
       
        
       % (C4) Beampattern check
       W_IASC_k_opt = [squeeze(w_i_opt(:,:,k)), w_t_opt(:,k)];
       W_total_k_opt = [W_IASC_k_opt,z_opt(:,k)];
       beampattern_true=real(a_phi_k_H_list{k} * W_total_k_opt * W_total_k_opt'*a_phi_k_H_list{k}');
       beampattern_fake=real(a_phi_k_H_list_fake{k} * W_total_k_opt * W_total_k_opt'*a_phi_k_H_list_fake{k}');
       beampattern_difference=beampattern_fake-beampattern_true;
       min_beampattern_difference= min(min_beampattern_difference,beampattern_difference);
       fprintf('True-target beampattern: %.7e\n', beampattern_true);
       fprintf('Fake-target beampattern: %.7e\n',beampattern_fake);

    end

    % Check (C2) power.
    power_z=0;
    miu_z=0;
    power_total=0;
    power_1=0;
    power_2=0;
    for k = 1:Ntx
        power_k_i1 = sum(sum_square_abs(w_i_opt(:, 1, k)));
        power_k_i2 = sum(sum_square_abs(w_i_opt(:, 2, k)));
        power_k_i3 = sum(sum_square_abs(w_i_opt(:, 3, k)));
        power_k_i4 = sum(sum_square_abs(w_i_opt(:, 4, k)));
        power_k_t  = sum_square_abs(w_t_opt(:, k));
        power_k_z  = sum_square_abs(z_opt(:, k));
        fprintf('    [Check] --- AP %.1f power (constraint C2) ---\n',k);
        fprintf('    [Check] User 1 communication power: %.7f \n,User 2 communication power: %.7f \n,User 3 communication power: %.7f \n,User 4 communication power: %.7f \n,Sensing power: %.7f \n,AN power %.7f \n', power_k_i1 ,power_k_i2,power_k_i3,power_k_i4,power_k_t,power_k_z );
        % fprintf('    [Check] User 1 communication power: %.7f \n,User 2 communication power: %.7f \n,Sensing power: %.7f \n,AN power %.7f \n', power_k_i1 ,power_k_i2,power_k_t,power_k_z );
        fprintf('    [Check]            Target power: %.3f \n', P_T_lin);
        power_total_k=power_k_i1+power_k_i2+power_k_i3+power_k_i4+power_k_t+power_k_z;
        miu_k_z=power_k_z/power_total_k;
        miu_z=miu_z+miu_k_z;
        power_z=power_z+power_k_z;
        power_total=power_total+power_total_k;
        power_1=power_1+power_k_i1;
        power_2=power_2+power_k_i2;
    end
    power_1=power_1/Ntx;
    power_2=power_2/Ntx;
    power_total=power_total/Ntx;
    power_z=power_z/Ntx;
    miu_z=miu_z/Ntx;
    val_noz = 0;
    for k = 1:Ntx
        W_ISAC_k_opt = [squeeze(w_i_opt(:,:,k)), w_t_opt(:,k)];
        val_noz = val_noz + c_k_vec(k) * real(a_phi_k_H_list{k} * W_ISAC_k_opt * W_ISAC_k_opt'*a_phi_k_H_list{k}');
    end
    val_noz=val_noz / sigma_ap_linear;
    % --- Print a unified verification summary ---
    fprintf('    [Check] Final sensing SNR: %.7f dB (%.7e linear)\n', 10*log10(gamma_t_opt), gamma_t_opt);
    fprintf('    [Check] Sensing SNR excluding z contribution: %.7f dB (%.7e linear)\n',10*log10(val_noz), val_noz);
    fprintf('    [Check] --- Communication SINR (constraint C1) ---\n');
    fprintf('    [Check] Worst-case (minimum) SINR: %.7f dB\n', min_sinr_dB);
    fprintf('    [Check]            Target SINR: %.7f dB\n', Gamma_comm_dB);
   
    
 
    fprintf('    [Check] --- Null-space leakage (constraint C3) ---\n');

    fprintf('    [Check] Worst-case (maximum) leakage: %.7e (linear)\n', max_leakage_lin);

    fprintf('    [Check] Minimum beampattern difference (constraint C4): %.7e\n',min_beampattern_difference);
    % 
    % 
    
end
