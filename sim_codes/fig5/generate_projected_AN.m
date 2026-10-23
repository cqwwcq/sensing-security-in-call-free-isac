function Z_AN = generate_projected_AN(Nue,Mtx,Ntx,H, U_users, P_max_list, ratio_AN)
% GENERATE_PROJECTED_AN 生成基于零空间投影的随机人工噪声(AN)
% 
% 核心原理:
%   利用投影矩阵 P = I - G' * pinv(G') 将随机噪声投影到联合等效信道 G 的零空间。
%   确保 AN 对所有用户的接收波束合并后的信号正交 (Ui*Hi*z = 0)。
%
% 输入参数:
%   H           :  Cell数组, 物理信道矩阵. 
%   U_users     : [Mue x Nue] 矩阵, 所有用户的接收波束向量
%   P_max_list  : [Ntx x 1] 向量 (或标量), 每个 AP 的最大发射功率
%   ratio_AN    : 标量 (0~1), AN 功率占总功率的比例 (例如 0.25)
%   Mtx         : 标量, AP 的发射天线数量
%
% 输出参数:
%   Z_AN        : [Mtx x Ntx] 矩阵. 第 k 列是第 k 个 AP 的 AN 向量 z_k

    % --- 1. 参数提取与初始化 ---
    
        P_max_list = P_max_list * ones(Ntx, 1);
    
    
    % 初始化输出矩阵
    Z_AN = zeros(Mtx, Ntx);
    
    % --- 2. 逐个 AP 生成 AN ---
    for k = 1:Ntx
        
        % A. 构建联合等效信道矩阵 G_k (Nue x Mtx)
        % G_k 的每一行代表一个用户对该 AN 的"接收敏感方向"
        G_k = zeros(Nue, Mtx);
        
        for i = 1:Nue
            H_ik = H(:,(1+(k-1)*Mtx):k*Mtx,i);   % [Mue x Mtx]
            u_i  = U_users(:, i); % [Mue x 1]
           
            % 计算等效信道向量 g^H = u^H * H
            G_k(i, :) = u_i' * H_ik; 
        end
        
        % B. 自由度检查
        % 如果天线数 Mtx <= 用户数 Nue，零空间不存在或唯一的 0
        if Mtx <= Nue
            warning('AP %d: 天线数(%d) <= 用户数(%d)，无自由度生成正交AN，跳过。', k, Mtx, Nue);
            Z_AN(:, k) = zeros(Mtx, 1);
            continue;
        end
        
        % C. 计算投影矩阵 P_proj
        % 公式: P = I - G^dagger * G
        % 使用 pinv() 求伪逆以保证数值稳定性
        P_proj = eye(Mtx) - pinv(G_k) * G_k;
        
        % D. 生成全维随机噪声 r (Complex Gaussian)
        % r 是各向同性的，不包含任何指向性意图
        r_full = (randn(Mtx, 1) + 1j * randn(Mtx, 1));
        
        % E. 执行投影
        % 切除指向用户的分量，只保留零空间分量
        z_raw = P_proj * r_full;
        
        % F. 功率归一化
        P_target = P_max_list(k) * ratio_AN;
        current_energy = norm(z_raw)^2;
        
        if current_energy > 1e-12
            scale_factor = sqrt(P_target / current_energy);
            z_final = scale_factor * z_raw;
        else
            z_final = zeros(Mtx, 1);
        end
        
        % G. 存入结果
        Z_AN(:, k) = z_final;
        
        % (可选调试: 验证正交性)
        % residue = norm(G_k * z_final);
        % if residue > 1e-6
        %     fprintf('Warning: AP %d 正交性误差 %.2e\n', k, residue);
        % end
    end
end