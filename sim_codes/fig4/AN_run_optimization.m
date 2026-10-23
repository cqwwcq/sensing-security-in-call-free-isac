function [gamma_t_opt, w_i_opt, w_t_opt,u_opt,z_opt, solve_time, success,power_z,miu_z,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense] = AN_run_optimization(...
    Ntx, Nue, Mtx, Nrx,Mue,P_T_dBm, Gamma_comm_dB, sigma_noise_dBm, sigma_ap_dBm,epsilon, ...
    H_all, Txbeamsteering, Fake_Txbeamsteering,sigmarcs, betark,rou)
    
    tic; % 启动求解计时器
    
    % --- 1. 转换参数为线性尺度 ---
    P_T_lin = 10^(P_T_dBm / 10) / 1000;      % W
    Gamma_comm_lin = 10^(Gamma_comm_dB / 10);
    sigma_noise_lin = 10^(sigma_noise_dBm / 10) / 1000; % 
    sigma_ap_linear = 10^(sigma_ap_dBm / 10) / 1000; % 
    epsilon_lin = epsilon; % 假设 epsilon 已经是线性尺度


    % --- 2. 预计算所有信道和矩阵 ---
    % (这是将您的参数映射到优化问题公式的关键步骤)
    
    h_ik_row_list = cell(Nue, Ntx);     % (i, k) 存储 1xMtx 行向量 h_ik
    a_phi_k_H_list = cell(Ntx, 1);        % (k) 存储 1xMtx 行向量 a_phi_k^H
    a_phi_k_H_list_fake = cell(Ntx,1);   %(k) 存储 1xMtx 行向量a_phi_k_fake^H
    a_theta_ik_H_list = cell(Nue, Ntx); % (i, k) 存储 1xMtx 行向量 a_theta_ik^H
    for k = 1:Ntx
        % 2a. 传感目标向量 a_phi_k^H  
        a_phi_k = squeeze(Txbeamsteering(:,k,:)); % Mtx x 1 row vector
        a_phi_k_H_list{k} = a_phi_k'; % 1 x Mtx column vector
        a_phi_k_fake = squeeze(Fake_Txbeamsteering(:,k,:));
        a_phi_k_H_list_fake{k} = a_phi_k_fake'; % 1 x Mtx column vector
        for i = 1:Nue
            % 2b. 通信信道 h_ik 和 Q_ik
            % 提取从 AP k 到 User i 的信道
            H_list_raw{i, k} = H_all(:, (1+(k-1)*Mtx):k*Mtx, i); % (MuexMtx)    
          
            
            
        end
    end
    
     
    % 2d. 传感目标函数权重 c_k
    c_k_vec = zeros(1, Ntx);
    for k = 1:Ntx
    % c_k_vec(k) = (1 / (Nrx * sigma_ap_linear)) * sum(sigmarcs^2 * betark(k, :));
        c_k_vec(k) = (1 / (Nrx )) * sum(sigmarcs^2 * betark(k, :));

    end
 

    % --- 3. 阶段一：寻找可行点 (P-Feas-CCCP) ---
    % fprintf('    [阶段一] 寻找可行点...\n');
    % p_init_total = P_T_lin * 0.001; 
    % p_init_per_beam = p_init_total / (Ntx * (Nue + 2));% 分配给波束的初始功率，包括 AN
    % scale = sqrt(p_init_per_beam); % 假设每个波束的功率
    % (Mtx x Nue x Ntx)
    % w_i_r = (randn(Mtx, Nue, Ntx) + 1j*randn(Mtx, Nue, Ntx)) * scale / sqrt(Mtx);
    % (Mtx x Ntx)
    % w_t_r = (randn(Mtx, Ntx) + 1j*randn(Mtx, Ntx)) * scale / sqrt(Mtx);
     p_wi_per_link=P_T_lin * 0.6/Nue*0.8;
     p_kickstart = P_T_lin * 0.3*0.8; 
     P_an=P_T_lin * 0.1*0.8; 
    for k = 1:Ntx
         for i = 1:Nue
            H_mat  = H_list_raw{i, k};       
            % 共轭转置得到波束方向 (Mtx x 1)
            % 注意：如果 h 是行向量，h' 就是列向量
            [~, ~, V] = svd(H_mat);         
            dominant_dir = V(:, 1); % 大小为 (Mtx x 1)
            % 归一化并赋值
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
    %         % 定义 CVX 变量
    %         variable w_i_vars(Mtx, Nue, Ntx) complex;
    %         variable w_t_vars(Mtx, Ntx) complex;
    %         variable s_var; % 松弛变量
    %         variable z_vars(Mtx, Ntx) complex;
    % 
    %         minimize(s_var);
    % 
    %         subject to
    %             for k = 1:Ntx
    %                 % (C2) 总功率约束 (per-AP)
    % 
    %                 power_expr_k = sum(sum_square_abs(w_i_vars(:, :, k))) + sum_square_abs(w_t_vars(:, k))+sum_square_abs(z_vars(:, k));
    %                 power_expr_k <= P_T_lin;
    %                 % 零空间约束 (C3)：人工噪声总和落在 U 的零空间
    %                 for i = 1:Nue
    %                     an_ik =h_ik_row_list{i, k} * z_vars(:, k);
    %                     abs(an_ik)<= epsilon_lin;
    %                 end
    %                % an_sense_leakage = a_phi_k_H_list{k} * z_vars(:, k);
    %                %  abs(an_sense_leakage) <= epsilon_lin; 
    %                 %(C4):波束图约束                      
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
    %                 % (C1'_feas) QoS 约束
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
    %         fprintf('    [阶段一] 迭代 %d: 求解失败。Status: %s\n', r, cvx_status);
    %         break;
    %     end
    % 
    %     fprintf('    [阶段一] 迭代 %d: s = %.6e\n', r, s_var);
    % 
    %     if s_var <= 0
    %         fprintf('    [阶段一] 找到可行解！\n');
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
    %     fprintf('    [阶段一] 未能找到可行解。算法终止。\n');
    %     gamma_t_opt = NaN;
    %     w_i_opt = {};
    %     w_t_opt = {};
    %     z_opt = {};
    %     solve_time = toc;
    %     success = false;
    %     return; 
    % end
    % --- 4. 阶段二：最大化传感 SINR (交替优化)(P-CCCP) ---


    
    MAX_ITER_MAIN = 30;
    last_obj_val = -Inf;
    fprintf('    [阶段二] 开始优化主问题...\n');

    for p = 1:MAX_ITER_MAIN
        u_opt = AN_update_receiver_mmse(Ntx, Nue, Mue, Mtx, H_list_raw, w_i_r, w_t_r, z_r, sigma_noise_lin);
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
            % --- 构建目标函数 ---
 
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
                    % (C2) 总功率约束 (per-AP)
                   
                    power_expr_k = sum(sum_square_abs(w_i_vars(:, :, k))) + sum_square_abs(w_t_vars(:, k))+ ...
                               sum_square_abs(z_vars(:, k));
                    power_expr_k <= P_T_lin;
                    %(C3)零空间约束
                    for i = 1:Nue
                        an_ik =h_ik_row_list{i, k} * z_vars(:, k);
                        abs(an_ik) <= epsilon_lin;

                    end
                    % an_sense_leakage = a_phi_k_H_list{k} * z_vars(:, k);
                    % abs(an_sense_leakage) <= epsilon_lin; 
                    %(C4):波束图约束                      
                    w_total_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k), z_vars(:,k)];
                    w_total_k_r = [squeeze(w_i_r(:,:,k)), w_t_r(:,k), z_r(:,k)];
                    W_ISAC_k_var = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k)];                 
                    expression grad_fake_term
                    grad_fake_term =2 * real( a_phi_k_H_list_fake{k} *  z_vars(:,k) * z_r(:,k)'*a_phi_k_H_list_fake{k}');
                    const_fake_term = real( a_phi_k_H_list_fake{k} * z_r(:,k) * z_r(:,k)'*a_phi_k_H_list_fake{k}' );
                    f=pow_pos( norm( a_phi_k_H_list{k} * W_ISAC_k_var, 'fro' ), 2 );
                    rou*f-grad_fake_term+const_fake_term<=0;           

            
                end
                for i = 1:Nue
                    % (C1'_feas) QoS 约束
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
            
                fprintf('    [阶段二] 迭代 %d: 求解失败。Status: %s\n', p, cvx_status);
                % 如果在优化过程中失败，我们回退到上一个可行的解
 
                break;
          
        end
   
        
        % --- 检查收敛 ---
        current_obj_val = 0;
        for k = 1:Ntx
            w_k_current = [squeeze(w_i_vars(:,:,k)), w_t_vars(:,k) ,z_vars(:,k)];
            current_obj_val = current_obj_val + c_k_vec(k) * real(a_phi_k_H_list{k} * w_k_current * w_k_current'*a_phi_k_H_list{k}');
        end
      


        
        fprintf('    [阶段二] 迭代 %d:  传感 SNR = %.8f dB\n', p, 10*log10(current_obj_val/ sigma_ap_linear));
        
        if abs(current_obj_val - last_obj_val) / (1e-9 + abs(last_obj_val)) < 4e-3
            fprintf('    [阶段二] 算法已收敛。\n');
            w_i_r = w_i_vars; 
            w_t_r = w_t_vars; 
            z_r = z_vars;
            last_obj_val = current_obj_val;
            break;
        end
        
        last_obj_val = current_obj_val;
        w_i_r = w_i_vars;
        w_t_r = w_t_vars;
        z_r = z_vars;
        
        if p == MAX_ITER_MAIN
            fprintf('    [阶段二] 达到最大迭代次数。\n');
        end
    end
 
        

  

    % --- 5. 最终结果 ---
    gamma_t_opt = last_obj_val/ sigma_ap_linear;
    w_i_opt = w_i_r;
    w_t_opt = w_t_r;
    z_opt = z_r;
    solve_time = toc; % 停止计时器
    success = true;
   
      % --- 6. *** 新增：最终约束验证 (详细报告版) *** ---
    fprintf('    [验证] 正在检查最终解的约束...\n');
    

    
    % 初始化“最差情况”的跟踪变量
    min_sinr_dB = Inf;
    max_leakage_lin = -Inf;
    max_power_lin = -Inf;
    min_beampattern_difference= Inf;
    interf_UE=0
    interf_sense=0;
    % 检查 (C1) 通信SINR 和 (C3) 传感泄漏
    for i = 1:Nue
        % 使用优化后的 U 验证 SINR
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


    for k = 1:Ntx
        for i = 1:Nue
          
            % % (C3) 零空间泄露
            h_ik_row_list{i, k} = u_opt{i}' * H_list_raw{i, k};

            leakage_ik_lin =h_ik_row_list{i, k} * z_opt(:, k);

            % 
            % % 跟踪最大泄漏
            max_leakage_lin = max(max_leakage_lin, leakage_ik_lin);
        end
        
       %(C4)波束图检测
       W_IASC_k_opt = [squeeze(w_i_opt(:,:,k)), w_t_opt(:,k)];
       W_total_k_opt = [W_IASC_k_opt,z_opt(:,k)];
       beampattern_true=real(a_phi_k_H_list{k} * W_total_k_opt * W_total_k_opt'*a_phi_k_H_list{k}');
       beampattern_fake=real(a_phi_k_H_list_fake{k} * W_total_k_opt * W_total_k_opt'*a_phi_k_H_list_fake{k}');
       beampattern_difference=beampattern_fake-beampattern_true;
       min_beampattern_difference= min(min_beampattern_difference,beampattern_difference);
       fprintf('真目标波束图: %.7e\n', beampattern_true);
       fprintf('假目标波束图: %.7e\n',beampattern_fake);

    end

    % 检查 (C2) 功率
    power_z=0;
    miu_z=0;
    power_total=0;
    power_1=0;
    power_2=0;
    for k = 1:Ntx
        power_k_i1 = sum(sum_square_abs(w_i_opt(:, 1, k)));
        power_k_i2 = sum(sum_square_abs(w_i_opt(:, 2, k)));
        % power_k_i3 = sum(sum_square_abs(w_i_opt(:, 3, k)));
        % power_k_i4 = sum(sum_square_abs(w_i_opt(:, 4, k)));
        power_k_t  = sum_square_abs(w_t_opt(:, k));
        power_k_z  = sum_square_abs(z_opt(:, k));
        fprintf('    [验证] --- AP%.1f功率 (约束 C2) ---\n',k);
        fprintf('    [验证] 用户1通信功率: %.7f \n,用户2通信功率: %.7f \n,感知功率: %.7f \n,AN功率 %.7f \n', power_k_i1 ,power_k_i2,power_k_t,power_k_z );
        fprintf('    [验证]            目标功率: %.3f \n', P_T_lin);
        power_total_k=power_k_i1+power_k_i2+power_k_t+power_k_z;
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
    % --- 统一输出验证结果 ---
    fprintf('    [验证] 最终传感 SNR: %.7f dB (%.7e 线性)\n', 10*log10(gamma_t_opt), gamma_t_opt);
    fprintf('    [验证]不包含z贡献的传感SNR:%.7f dB(%.7e 线性)\n',10*log10(val_noz), val_noz);
    fprintf('    [验证] --- 通信 SINR (约束 C1) ---\n');
    fprintf('    [验证] 最差 (最小) SINR: %.7f dB\n', min_sinr_dB);
    fprintf('    [验证]            目标 SINR: %.7f dB\n', Gamma_comm_dB);
   
    
 
    fprintf('    [验证] --- 零空间泄漏 (约束 C3) ---\n');

    fprintf('    [验证] 最差 (最大) 泄漏: %.7e (线性)\n', max_leakage_lin);

    fprintf('    [验证] (约束 C4)最小波束图差: %.7e\n',min_beampattern_difference);
    % 
    % 
    
end