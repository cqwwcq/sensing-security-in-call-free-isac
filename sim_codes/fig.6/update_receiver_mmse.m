function u_new = update_receiver_mmse(Ntx, Nue, Mue, Mtx, H_list, w_i, w_t, sigma_n)
    u_new = cell(Nue, 1);
    
    for i = 1:Nue
        % 1. 计算总接收信号协方差 R_yy (Mue x Mue)
        % R_yy = Sum_l P_l + P_sense + Noise
        R_yy = sigma_n * eye(Mue);
        
        % 所有用户的通信信号
        for l = 1:Nue
            if l~=i
            % h_l_agg = sum_k H_{i,k} * w_{l,k}
            h_w_l = zeros(Mue, 1);
            for k = 1:Ntx
                h_w_l = h_w_l + H_list{i, k} * w_i(:, l, k);
            end
            R_yy = R_yy + h_w_l * h_w_l';
            end
        end
        
        % 感知信号
        h_w_t = zeros(Mue, 1);
        for k = 1:Ntx
            h_w_t = h_w_t + H_list{i, k} * w_t(:, k);
        end
        R_yy = R_yy + h_w_t * h_w_t';
        
        % 2. 计算期望信号向量 v_i
        v_i = zeros(Mue, 1);
        for k = 1:Ntx
            v_i = v_i + H_list{i, k} * w_i(:, i, k);
        end
        
        % 3. MMSE 解: u = R_yy \ v_i
        % (注意：标准的 Max-SINR 接收机是 J^-1 * v，其中 J 是干扰+噪声。
        %  但在 MMSE 中，使用 R_yy^-1 * v 得到的方向是一样的，仅相差一个标量)
        u_vec = R_yy \ v_i;
        
        % 归一化
        if norm(u_vec) > 0
            u_new{i} = u_vec / norm(u_vec);
        else
            u_new{i} = ones(Mue, 1) / sqrt(Mue);
        end
    end
end
