clear;
clc;
close all;

%% ======================== 1. 仿真参数设置 ========================
% 说明：
% 1) 该脚本直接复用现有 Parameters / run_optimization / AN_run_optimization /
%    Em / AN_Em 等函数；
% 2) 互信息的随机变量定义为：
%    Theta   -> 某个 AP 上按网格穿越规则定义的真实方向类索引；
%    Theta^f -> 用于塑造假峰的假目标离散方向索引；
%    Theta^p -> 窃听者基于 EM 恢复信号检测到的方向类索引；
% 3) 统计时对全部 (Theta, Theta^f) 有效组合进行边缘化，再计算
%    每个 AP 上的 I(Theta; Theta^p)。

% ---------------- 系统维度参数 ----------------
Ntx = 8;                 % 发射 AP 数量
Mtx = 16;                % 每个 AP 的天线数
Nrx = 4;                 % 接收 AP 数量
Nue = 2;                 % 合法用户数
Mue = 16;                % 用户天线数，用户 1 视为内部窃听者

% ---------------- 优化与噪声参数 ----------------
P_T = 50;                % dBm
Gamma_comm = 10;         % dB
sigma_noise = -94;       % dBm
sigma_ap = -94;          % dBm
epsilon = 1e-10;
nvariance = 10^(sigma_noise / 10) / 1000;
rou_noAN=0;
rou_AN=1.3;
% ---------------- 目标状态与方向离散化参数 ----------------
% 默认真目标状态仍沿对角线取值，假目标状态默认与真目标共享同一候选集合，
% 但统计时不再做一一对应，而是遍历所有有效 (真目标, 假目标) 组合。
use_full_grid_as_true_state = false;
use_full_grid_as_fake_state = false;
grid_step = 20;
background_bin_deg = 15;

xcont = -195:grid_step:195;
ycont = 195:-grid_step:-195;
%
% xtarget_mat = -195:20:185;
% ytarget_mat = -xtarget_mat;


target_axis = -120:60:120;
[xtarget_mat, ytarget_mat] = meshgrid(target_axis, target_axis);


if use_full_grid_as_true_state
    [Xcont_true, Ycont_true] = meshgrid(xcont, ycont);
    PosTarget_all = [Xcont_true(:), Ycont_true(:)];
else
    PosTarget_all = [xtarget_mat(:), ytarget_mat(:)];
end
% PosTarget_all: 所有真目标离散状态的坐标集合，每一行是一个真目标位置 [x, y]。

N_target_state = size(PosTarget_all, 1);
if use_full_grid_as_fake_state
    [Xcont_fake, Ycont_fake] = meshgrid(xcont, ycont);
    PosFake_candidate_all = [Xcont_fake(:), Ycont_fake(:)];
else
    PosFake_candidate_all = PosTarget_all;
end
% PosFake_candidate_all: 假目标候选状态集合。
% 有 AN 情形下，对每个真目标状态都要遍历这里面的所有允许假目标状态。

N_fake_state = size(PosFake_candidate_all, 1);
valid_fake_idx_cell = cell(N_target_state, 1);
N_state_pair = 0;
for idx_target = 1:N_target_state
    PosTarget = PosTarget_all(idx_target, :);
    valid_fake_idx = find(~all(abs(PosFake_candidate_all - PosTarget) < 1e-12, 2));
    valid_fake_idx_cell{idx_target} = valid_fake_idx;
    N_state_pair = N_state_pair + numel(valid_fake_idx);
end
% valid_fake_idx_cell{i}: 给定第 i 个真目标状态时，所有合法假目标状态的索引。
% 这里排除了“假目标位置与真目标位置完全相同”的组合。
% N_state_pair: 全部合法 (真目标, 假目标) 组合的总数。

% ---------------- Monte Carlo 参数 ----------------
% 说明：该脚本需要重复调用优化器，运行时间较长。
% 若只做快速检查，可先把 num_MC 设小一些，例如 2 或 3。
% num_MC 的作用是对重复检测结果做平均，而不是把不同 Monte Carlo 的 Rall 先平均。
num_MC = 1;

% Em / AN_Em 内部固定恢复 128 个快拍；
% 对每个 Monte Carlo 样本，先把这些快拍对应的 Rall 在第三维求平均，
% 再基于平均后的协方差做一次方向检测。
% 若希望只用部分恢复时刻参与平均，可修改 snapshot_idx_vec。
snapshot_idx_vec = 1:128;
% snapshot_idx_vec: 参与统计的快拍编号集合。
% 若设成 1:128，表示把 EM 恢复出的所有快拍都先求均值后再做检测。

% 波束图扫描角度，与现有 LinesCont.m 保持一致。
theta_scan_num = 3604;
thetal_scan = linspace(0, 2 * pi, theta_scan_num);
asteer_scan = local_build_asteer_scan(thetal_scan, Mtx);

% 随机数种子，便于复现实验。
rng(20260424);

%% ======================== 2. 预分配统计量 ========================
PosAPs = [-200,-200;200,200;-200,200;200,-200;0,-200;0,200;-200,0;200,0];
% PosAPs: 8 个发射 AP 的坐标，和 Parameters.m 中的默认设置一致。


direction_state_cell = local_build_direction_state(PosAPs, PosTarget_all, grid_step, background_bin_deg);
% direction_state_cell{k}: 第 k 个 AP 对应的方向类描述结构体。
% 它定义了“哪个真目标网格属于哪个方向类”，以及每个方向类的角区间范围。

count_noAN_cell = cell(Ntx, 1);
count_AN_cell = cell(Ntx, 1);
for k = 1:Ntx
    N_theta_k = direction_state_cell{k}.num_state;
    count_noAN_cell{k} = zeros(N_theta_k, N_theta_k);
    count_AN_cell{k} = zeros(N_theta_k, N_theta_k);
end
% count_noAN_cell{k} / count_AN_cell{k}: 第 k 个 AP 的联合计数矩阵。
% 行索引表示真实方向类 Theta，列索引表示检测方向类 Theta^p。
% 两者分别用于无 AN 和有 AN 情形的统计。

fail_noAN = 0;
fail_AN = 0;
success_noAN = 0;
success_AN = 0;
% success_* / fail_*: 记录优化成功或失败的样本数。
% 失败样本不会进入后续互信息统计。

%% ======================== 3. 主循环：生成样本并统计方向转移 ========================
tic;
fprintf('开始互信息对比仿真：共 %d 个真目标状态，%d 个假目标候选状态，%d 个有效组合。\n', ...
    N_target_state, N_fake_state, N_state_pair);

for idx_target = 1:N_target_state
    PosTarget = PosTarget_all(idx_target, :);
    valid_fake_idx = valid_fake_idx_cell{idx_target};

    fprintf('\n[真目标状态 %d/%d] PosTarget = [%.1f, %.1f], 有效假目标组合数 = %d\n', ...
        idx_target, N_target_state, PosTarget(1), PosTarget(2), numel(valid_fake_idx));

    % 与 LinesCont.m 保持一致：生成一次场景，再分别进行无 AN 与有 AN 对比。
    sigma2 = sqrt(nvariance / 2) * (randn(Mue, 1) + 1j * randn(Mue, 1));
    % sigma2: 当前样本下，窃听者接收端使用的复高斯噪声向量。

    if isempty(valid_fake_idx)
        PosFake_reference = local_get_offset_fake_target(PosTarget);
    else
        PosFake_reference = PosFake_candidate_all(valid_fake_idx(1), :);
    end
    % PosFake_reference: 仅用于兼容 Parameters 接口的参考假目标坐标。
    % 无 AN 情形下它不参与误导设计，只用于生成与接口匹配的参数。

    [lambdac, betark, Dsn, PosAPs, Posue, PosRx, distUe, distTx, distRx, ...
        Txbeamsteering, Rxbeamsteering, Userbeamsteering, Fake_Txbeamsteering, ...
        alphark, sigmarcs, H, Z_shd] = ...
        Parameters(Ntx, Nrx, Mtx, Nue, Mue, PosTarget, PosFake_reference);

    %#ok<NASGU>
    % 变量 lambdac / Posue / PosRx / distUe / distTx / distRx / alphark / Z_shd
    % 在该脚本中不直接使用，但保留输出是为了和现有参数接口保持一致。

    % ---------------- 无 AN 情形 ----------------
    [~, w_i_noAN, w_t_noAN, ~, ~, status_noAN, ~, ~, ~, ~, ~, ~] = run_optimization(...
        Ntx, Nue, Mtx, Nrx, Mue, P_T, Gamma_comm, sigma_noise, sigma_ap, epsilon, ...
        H, Txbeamsteering, Userbeamsteering, sigmarcs, betark,rou_noAN);

    if status_noAN
        success_noAN = success_noAN + 1;
        W_noAN = local_pack_W_tensor(w_i_noAN, w_t_noAN, Mtx, Ntx);
        Rall_noAN = Em(Dsn, nvariance, Ntx, Nue, Mue, H, sigma2, W_noAN, Mtx);
        % W_noAN: 无 AN 情形下送入 Em 的发送波束张量，
        %         维度为 Mtx x (Nue+1) x Ntx，前 Nue 列是通信波束，最后一列是感知波束。
        % Rall_noAN: Em 输出的恢复协方差集合，
        %           维度为 Mtx x Mtx x 快拍 x AP。

        count_noAN_cell = local_update_direction_count(...
            count_noAN_cell, Rall_noAN, idx_target, direction_state_cell, ...
            snapshot_idx_vec, thetal_scan, asteer_scan, max(numel(valid_fake_idx), 1));
    else
        fail_noAN = fail_noAN + 1;
        fprintf('    无 AN 优化失败，本次样本跳过。\n');
    end

    for idx_fake_local = 1:numel(valid_fake_idx)
        idx_fake = valid_fake_idx(idx_fake_local);
        Fake_Target_Pos = PosFake_candidate_all(idx_fake, :);
        [angle_Fake_Tx, ~] = compute_angle_dist(PosAPs, Fake_Target_Pos);
        Fake_Txbeamsteering = beamsteering(angle_Fake_Tx.', Mtx);

        [~, w_i_AN, w_t_AN, ~, z_AN, ~, status_AN, ~, ~, ~, ~, ~, ~, ~, ~] = AN_run_optimization(...
            Ntx, Nue, Mtx, Nrx, Mue, P_T, Gamma_comm, sigma_noise, sigma_ap, epsilon, ...
            H, Txbeamsteering, Fake_Txbeamsteering, sigmarcs, betark,rou_AN);

        if status_AN
            success_AN = success_AN + 1;
            W_AN = local_pack_W_tensor(w_i_AN, w_t_AN, Mtx, Ntx);
            Rall_AN = AN_Em(Dsn, nvariance, Ntx, Nue, Mue, H, sigma2, W_AN, z_AN, Mtx);
            % W_AN: 有 AN 情形下送入 AN_Em 的通信+感知波束张量。
            % z_AN: 人工噪声波束矩阵，按 AP 分列存放。
            % Rall_AN: 有 AN 情形下，EM 恢复得到的发射协方差集合。

            count_AN_cell = local_update_direction_count(...
                count_AN_cell, Rall_AN, idx_target, direction_state_cell, ...
                snapshot_idx_vec, thetal_scan, asteer_scan, 1);
        else
            fail_AN = fail_AN + 1;
            fprintf('    有 AN 优化失败，真目标索引 %d，假目标索引 %d。\n', idx_target, idx_fake);
        end
    end
end

elapsed_time = toc;

%% ======================== 4. 由经验计数计算互信息 ========================
result_noAN = local_summarize_MI(count_noAN_cell);
result_AN = local_summarize_MI(count_AN_cell);
% result_noAN / result_AN: 两种情形下的统计结果结构体。
% 其中包含逐 AP 的互信息、经验条件分布、经验先验分布和命中率等量。

AP_index = (1:Ntx).';
MI_noAN_bit = result_noAN.MI_ap(:);
MI_AN_bit = result_AN.MI_ap(:);
MI_drop_bit = MI_noAN_bit - MI_AN_bit;
H_theta_noAN_bit = result_noAN.H_theta_ap(:);
H_theta_AN_bit = result_AN.H_theta_ap(:);
hit_noAN = result_noAN.hit_rate_ap(:);
hit_AN = result_AN.hit_rate_ap(:);
sample_num_noAN = result_noAN.sample_num_ap(:);
sample_num_AN = result_AN.sample_num_ap(:);
% MI_noAN_bit / MI_AN_bit: 每个 AP 上估计得到的互信息 I(Theta; Theta^p)。
% MI_drop_bit: AN 引入后互信息的下降量。
% H_theta_noAN_bit / H_theta_AN_bit: 真实方向类 Theta 的经验熵。
% hit_noAN / hit_AN: 真实方向类与检测方向类相同的比例，可理解为方向类判决命中率。
% sample_num_noAN / sample_num_AN: 每个 AP 实际参与统计的样本数。

results_MI = table(AP_index, MI_noAN_bit, MI_AN_bit, MI_drop_bit, ...
    H_theta_noAN_bit, H_theta_AN_bit, hit_noAN, hit_AN, ...
    sample_num_noAN, sample_num_AN);
% results_MI: 最终汇总表，每一行对应一个 AP，便于查看、导出和作图。

fprintf('\n======================== 互信息统计结果 ========================\n');
disp(results_MI);

fprintf('无 AN 平均互信息: %.6f bit\n', mean(MI_noAN_bit));
fprintf('有 AN 平均互信息: %.6f bit\n', mean(MI_AN_bit));
fprintf('平均互信息下降: %.6f bit\n', mean(MI_drop_bit));
fprintf('无 AN 成功样本数: %d, 失败样本数: %d\n', success_noAN, fail_noAN);
fprintf('有 AN 成功样本数: %d, 失败样本数: %d\n', success_AN, fail_AN);
fprintf('总耗时: %.2f s\n', elapsed_time);

%% ======================== 5. 结果可视化 ========================
figure('Name', 'Mutual Information', 'Color', 'w');
bar(AP_index, [MI_noAN_bit, MI_AN_bit], 'grouped');
grid on;
xlabel('number of AP');
ylabel('MI I(\Theta;\Theta^p) / bit');
title('Mutual Information');
legend('without AN', 'with AN', 'Location', 'best');

% figure('Name', '方向命中率对比 - 各 AP', 'Color', 'w');
% bar(AP_index, [hit_noAN, hit_AN], 'grouped');
% grid on;
% xlabel('AP 编号');
% ylabel('方向判决命中率');
% title('EM 恢复后主峰方向命中率对比');
% legend('无 AN', '有 AN', 'Location', 'best');

%% ======================== 6. 保存结果 ========================
save('MI_direction_detection_compare_results.mat', ...
    'results_MI', 'result_noAN', 'result_AN', 'direction_state_cell', ...
    'count_noAN_cell', 'count_AN_cell', 'PosTarget_all', 'PosFake_candidate_all', ...
    'valid_fake_idx_cell', 'N_target_state', 'N_fake_state', 'N_state_pair', ...
    'snapshot_idx_vec', 'grid_step', 'elapsed_time', ...
    'success_noAN', 'success_AN', 'fail_noAN', 'fail_AN');
% 保存结果时同时保留：
% 1) results_MI: 汇总表；
% 2) result_noAN / result_AN: 详细统计结构体；
% 3) direction_state_cell: 每个 AP 的方向类定义；
% 4) count_noAN_cell / count_AN_cell: 原始联合计数矩阵；
% 5) 目标状态集合和仿真参数，便于后续复现。


function direction_state_cell = local_build_direction_state(PosAPs, PosTarget_all, grid_step, background_bin_deg)
% 对每个 AP 构造方向类，规则尽量直接：
% 1) 目标分布范围内，默认每个目标方格各占一个方向状态；
% 2) 只有当两个目标方格从该 AP 看去几乎完全同角时，才将它们合并；
% 3) 为了容纳 EM 恢复误差，在目标方向之外再补一组较粗的背景方向分箱。
% 输出 direction_state_cell{k} 中各字段含义：
% theta_center_of_target: 每个真目标方格中心对应的角度。
% interval_list_target: 每个真目标方格对应的“射线穿过该网格”的角区间。
% state_index_of_target(i): 第 i 个真目标状态在该 AP 上属于哪个目标方向状态。
% target_center_state: 目标方向状态的中心角。
% target_interval_state: 目标方向状态对应的角区间。
% background_interval/background_center: 仅在目标方向并集之外构造的背景粗分箱角区间和中心角。
% background_state_offset: 背景状态编号相对于目标状态编号的偏移。
% num_target_state: 目标方向状态数量。
% num_state: 总方向状态数量 = 目标方向状态 + 背景方向状态。

Ntx_local = size(PosAPs, 1);
N_target_local = size(PosTarget_all, 1);
half_step = grid_step / 2;
background_bin_width = background_bin_deg * pi / 180;
duplicate_angle_tol = 1e-10;
direction_state_cell = cell(Ntx_local, 1);

for k = 1:Ntx_local
    [theta_center_of_target, ~] = compute_angle_dist(PosAPs(k, :), PosTarget_all);
    theta_center_of_target = local_wrapTo2Pi(theta_center_of_target(:));

    interval_list_target = zeros(N_target_local, 2);
    for idx_target_local = 1:N_target_local
        interval_list_target(idx_target_local, :) = local_compute_cell_angle_interval(...
            PosAPs(k, :), PosTarget_all(idx_target_local, :), half_step);
    end

    [theta_sorted, sort_idx] = sort(theta_center_of_target);
    target_member_cell = cell(N_target_local, 1);
    target_center_state = zeros(N_target_local, 1);
    target_state_count = 0;

    for idx_sorted = 1:N_target_local
        idx_target_local = sort_idx(idx_sorted);
        theta_now = theta_sorted(idx_sorted);

        if target_state_count == 0
            target_state_count = 1;
            target_member_cell{target_state_count} = idx_target_local;
            target_center_state(target_state_count) = theta_now;
        elseif abs(local_wrapToPi(theta_now - target_center_state(target_state_count))) <= duplicate_angle_tol
            target_member_cell{target_state_count}(end + 1) = idx_target_local;
        else
            target_state_count = target_state_count + 1;
            target_member_cell{target_state_count} = idx_target_local;
            target_center_state(target_state_count) = theta_now;
        end
    end

    if target_state_count > 1 && ...
            abs(local_wrapToPi(target_center_state(1) - target_center_state(target_state_count))) <= duplicate_angle_tol
        target_member_cell{1} = [target_member_cell{target_state_count}, target_member_cell{1}];
        target_member_cell(target_state_count) = [];
        target_center_state(target_state_count) = [];
        target_state_count = target_state_count - 1;
    else
        target_member_cell = target_member_cell(1:target_state_count);
        target_center_state = target_center_state(1:target_state_count);
    end

    state_index_of_target = zeros(N_target_local, 1);
    target_interval_state = zeros(target_state_count, 2);
    for idx_state = 1:target_state_count
        member_idx = target_member_cell{idx_state};
        state_index_of_target(member_idx) = idx_state;
        target_interval_state(idx_state, :) = local_merge_circular_intervals(interval_list_target(member_idx, :));
        target_center_state(idx_state) = local_circular_mean(theta_center_of_target(member_idx));
    end

    background_interval = local_build_background_intervals(target_interval_state, background_bin_width);
    background_center = mean(background_interval, 2);
    num_background_state = size(background_interval, 1);

    direction_state_cell{k}.theta_center_of_target = theta_center_of_target;
    direction_state_cell{k}.interval_list_target = interval_list_target;
    direction_state_cell{k}.state_index_of_target = state_index_of_target;
    direction_state_cell{k}.target_member_cell = target_member_cell;
    direction_state_cell{k}.target_center_state = target_center_state(:);
    direction_state_cell{k}.target_interval_state = target_interval_state;
    direction_state_cell{k}.num_target_state = target_state_count;
    direction_state_cell{k}.background_interval = background_interval;
    direction_state_cell{k}.background_center = background_center;
    direction_state_cell{k}.background_state_offset = target_state_count;
    direction_state_cell{k}.state_center = [target_center_state(:); background_center(:)];
    direction_state_cell{k}.num_state = target_state_count + num_background_state;
end
end


function asteer_scan = local_build_asteer_scan(thetal_scan, Mtx)
% 构造与 LinesCont.m 一致的波束扫描导向矩阵。
% 输出维度为 Mtx x Ntheta。
phase_mat = 1i * pi * ((1:(Mtx - 1)).') * sin(thetal_scan);
asteer_scan = [ones(1, numel(thetal_scan)); exp(phase_mat)];
end


function W_tensor = local_pack_W_tensor(w_i_opt, w_t_opt, Mtx, Ntx)
% 将 run_optimization / AN_run_optimization 的输出整理成 Em / AN_Em 需要的格式。
% 输入：
% w_i_opt: 通信波束，维度为 Mtx x Nue x Ntx。
% w_t_opt: 感知波束，维度为 Mtx x Ntx。
% 输出：
% W_tensor: 维度为 Mtx x (Nue+1) x Ntx，最后一列是感知波束。
W_tensor = cat(2, w_i_opt, reshape(w_t_opt, Mtx, 1, Ntx));
end


function count_cell = local_update_direction_count(count_cell, Rall, idx_target, ...
    direction_state_cell, snapshot_idx_vec, thetal_scan, asteer_scan, sample_weight)
% 基于 EM 恢复结果 Rall，逐 AP 对选定快拍先求平均，再检测主峰方向并更新经验联合计数。
% 关键输入含义：
% idx_target: 当前样本对应的真目标状态索引。
% Rall: Em 或 AN_Em 输出的恢复协方差集合。
% snapshot_idx_vec: 参与平均的快拍编号。
% sample_weight: 当前样本计入联合计数矩阵时的权重。
% 更新后的 count_cell{k}(i,j) 含义：
% 第 k 个 AP 上，真实方向类 = i 且检测方向类 = j 的累计权重。

if nargin < 9
    sample_weight = 1;
end

for k = 1:numel(direction_state_cell)
    theta_true_idx = direction_state_cell{k}.state_index_of_target(idx_target);
    R_temp = sum(Rall(:, :, snapshot_idx_vec, k), 3);
    theta_detect = local_detect_theta_from_R(R_temp, thetal_scan, asteer_scan);
    theta_detect_idx = local_find_direction_state_idx(theta_detect, direction_state_cell{k});
    count_cell{k}(theta_true_idx, theta_detect_idx) = count_cell{k}(theta_true_idx, theta_detect_idx) + sample_weight;
end
end


function theta_detect = local_detect_theta_from_R(R_temp, thetal_scan, asteer_scan)
% 对某个 AP 的平均恢复协方差做主峰方向检测。
% 输出 theta_detect 是连续角度值，不是方向类编号。
% 后续还需要再映射到离散方向类索引。

R_temp = (R_temp + R_temp') / 2;
temp_mat = asteer_scan' * R_temp;
py_vec = real(sum(temp_mat .* asteer_scan.', 2)).';
[~, ind_max] = max(py_vec);
theta_detect = thetal_scan(ind_max);
end


function theta_idx = local_find_direction_state_idx(theta_val, direction_state_struct)
% 将检测角度量化到方向状态：
% 1) 若角度落入某个目标方向状态区间，则归入最近的目标方向状态；
% 2) 若角度不在任何目标方向区间内，则归入较粗的背景方向状态；
% 3) 若仍存在边界误差，则回退到最近的状态中心。
% 输出 theta_idx 是方向类编号，后续作为联合计数矩阵的列索引。

theta_idx = [];
theta_val = local_wrapTo2Pi(theta_val);
num_target_state = direction_state_struct.num_target_state;
target_candidate_flag = false(num_target_state, 1);

for idx_state = 1:num_target_state
    interval_now = direction_state_struct.target_interval_state(idx_state, :);
    if local_is_angle_in_interval(theta_val, interval_now)
        target_candidate_flag(idx_state) = true;
    end
end

target_candidate = find(target_candidate_flag);

if ~isempty(target_candidate)
    delta_target = abs(local_wrapToPi(theta_val - direction_state_struct.target_center_state(target_candidate).'));
    [~, idx_min] = min(delta_target);
    theta_idx = target_candidate(idx_min);
    return;
end

for idx_bg = 1:size(direction_state_struct.background_interval, 1)
    interval_now = direction_state_struct.background_interval(idx_bg, :);
    if local_is_angle_in_interval(theta_val, interval_now)
        theta_idx = direction_state_struct.background_state_offset + idx_bg;
        return;
    end
end

delta_theta = abs(local_wrapToPi(theta_val - direction_state_struct.state_center.'));
[~, theta_idx] = min(delta_theta);
end


function result_struct = local_summarize_MI(count_cell)
% 由经验联合计数计算：
% 1) 联合概率 p(theta, theta_p)
% 2) 条件概率 q(theta_p | theta)
% 3) 边缘概率 p(theta), p(theta_p)
% 4) 互信息 I(theta; theta_p)
% 输出结构体 result_struct 的主要字段：
% count_cell: 原始联合计数矩阵。
% q_theta_cell{k}: 第 k 个 AP 的经验条件分布 q(theta^p | theta)。
% p_theta_cell{k}: 第 k 个 AP 的真实方向类经验先验分布。
% r_theta_p_cell{k}: 第 k 个 AP 的检测方向类经验边缘分布。
% MI_ap(k): 第 k 个 AP 的互信息估计值。
% H_theta_ap(k): 第 k 个 AP 的真实方向类经验熵。
% hit_rate_ap(k): 第 k 个 AP 的方向类命中率。
% sample_num_ap(k): 第 k 个 AP 实际参与统计的样本数量。

Ntx = numel(count_cell);

result_struct.count_cell = count_cell;
result_struct.q_theta_cell = cell(Ntx, 1);
result_struct.p_theta_cell = cell(Ntx, 1);
result_struct.r_theta_p_cell = cell(Ntx, 1);
result_struct.MI_ap = zeros(Ntx, 1);
result_struct.H_theta_ap = zeros(Ntx, 1);
result_struct.hit_rate_ap = zeros(Ntx, 1);
result_struct.sample_num_ap = zeros(Ntx, 1);

for k = 1:Ntx
    count_mat = count_cell{k};
    total_count = sum(count_mat(:));
    result_struct.sample_num_ap(k) = total_count;

    if total_count == 0
        result_struct.q_theta_cell{k} = zeros(size(count_mat));
        result_struct.p_theta_cell{k} = zeros(size(count_mat, 1), 1);
        result_struct.r_theta_p_cell{k} = zeros(1, size(count_mat, 2));
        continue;
    end

    joint_prob = count_mat / total_count;
    p_theta = sum(joint_prob, 2);
    r_theta_p = sum(joint_prob, 1);

    q_theta = zeros(size(count_mat));
    row_sum = sum(count_mat, 2);
    valid_row = row_sum > 0;
    q_theta(valid_row, :) = count_mat(valid_row, :) ./ row_sum(valid_row);

    MI_val = 0;
    [row_idx, col_idx] = find(joint_prob > 0);
    for idx = 1:numel(row_idx)
        i = row_idx(idx);
        j = col_idx(idx);
        MI_val = MI_val + joint_prob(i, j) * log2(joint_prob(i, j) / (p_theta(i) * r_theta_p(j)));
    end

    H_theta = -sum(p_theta(p_theta > 0) .* log2(p_theta(p_theta > 0)));
    hit_rate = trace(count_mat) / total_count;

    result_struct.q_theta_cell{k} = q_theta;
    result_struct.p_theta_cell{k} = p_theta;
    result_struct.r_theta_p_cell{k} = r_theta_p;
    result_struct.MI_ap(k) = max(real(MI_val), 0);
    result_struct.H_theta_ap(k) = real(H_theta);
    result_struct.hit_rate_ap(k) = real(hit_rate);
end
end


function theta_out = local_wrapTo2Pi(theta_in)
% 本地实现，避免依赖额外工具箱。
theta_out = mod(theta_in, 2 * pi);
theta_out(theta_out < 0) = theta_out(theta_out < 0) + 2 * pi;
end


function theta_out = local_wrapToPi(theta_in)
% 将角度差限制在 [-pi, pi)。
theta_out = mod(theta_in + pi, 2 * pi) - pi;
end


function interval_out = local_compute_cell_angle_interval(pos_ap, pos_target, half_step)
% 计算“从 AP 出发的射线穿过目标网格”对应的角度区间。
% 这里用目标网格四个顶点相对于 AP 的最小包络角区间来近似。
% 输入：
% pos_ap: 单个 AP 坐标 [x, y]。
% pos_target: 单个目标网格中心坐标 [x, y]。
% half_step: 网格边长的一半。
% 输出 interval_out = [theta_start, theta_end]。

corner_pos = [
    pos_target(1) - half_step, pos_target(2) - half_step;
    pos_target(1) - half_step, pos_target(2) + half_step;
    pos_target(1) + half_step, pos_target(2) - half_step;
    pos_target(1) + half_step, pos_target(2) + half_step];

[corner_angle, ~] = compute_angle_dist(pos_ap, corner_pos);
corner_angle = local_wrapTo2Pi(corner_angle(:).');
interval_out = local_min_cover_interval(corner_angle);
end


function interval_out = local_min_cover_interval(angle_vec)
% 在 [0, 2pi) 上为一组角度寻找最短包络区间。
% 若区间跨越 2pi，则展开到 [a, b] 且允许 b > 2pi，
% 后续统一使用 theta、theta+2pi 两份复制来判断重叠和包含关系。
% 这个函数的作用是把“网格四个顶点张成的角度范围”压缩成一个最短包络区间。

angle_sorted = sort(local_wrapTo2Pi(angle_vec));
angle_ext = [angle_sorted, angle_sorted(1) + 2 * pi];
angle_gap = diff(angle_ext);
[~, idx_gap_max] = max(angle_gap);
interval_start = angle_ext(idx_gap_max + 1);
interval_end = angle_ext(idx_gap_max) + 2 * pi;
interval_out = local_canonicalize_interval([interval_start, interval_end]);
end


function background_interval = local_build_background_intervals(target_interval_state, background_bin_width)
% 仅在目标方向状态并集之外构造背景粗分箱。
% 也就是说，背景方向只覆盖“目标方向区间没有覆盖到的角区域”。

two_pi = 2 * pi;
linear_interval = zeros(0, 2);

for idx_interval = 1:size(target_interval_state, 1)
    interval_now = target_interval_state(idx_interval, :);
    linear_interval = [linear_interval; local_interval_to_linear_segments(interval_now)]; %#ok<AGROW>
end

merged_interval = local_merge_linear_intervals(linear_interval);
gap_interval = local_build_gap_intervals(merged_interval, two_pi);
background_interval = local_split_intervals_by_width(gap_interval, background_bin_width);
end


function linear_interval = local_interval_to_linear_segments(interval_in)
% 将圆周角区间拆成 [0, 2pi] 上的线性区间段。

two_pi = 2 * pi;
interval_in = local_canonicalize_interval(interval_in);

if interval_in(2) <= two_pi
    linear_interval = interval_in;
else
    linear_interval = [interval_in(1), two_pi; 0, interval_in(2) - two_pi];
end
end


function merged_interval = local_merge_linear_intervals(linear_interval)
% 合并 [0, 2pi] 上互相重叠的线性区间。

if isempty(linear_interval)
    merged_interval = zeros(0, 2);
    return;
end

linear_interval = sortrows(linear_interval, 1);
merged_interval = zeros(size(linear_interval));
merged_count = 1;
merged_interval(merged_count, :) = linear_interval(1, :);

for idx_interval = 2:size(linear_interval, 1)
    interval_now = linear_interval(idx_interval, :);
    if interval_now(1) <= merged_interval(merged_count, 2)
        merged_interval(merged_count, 2) = max(merged_interval(merged_count, 2), interval_now(2));
    else
        merged_count = merged_count + 1;
        merged_interval(merged_count, :) = interval_now;
    end
end

merged_interval = merged_interval(1:merged_count, :);
end


function gap_interval = local_build_gap_intervals(merged_interval, two_pi)
% 构造目标方向并集之外的补集区间。

if isempty(merged_interval)
    gap_interval = [0, two_pi];
    return;
end

gap_interval = zeros(0, 2);

if merged_interval(1, 1) > 0
    gap_interval = [gap_interval; 0, merged_interval(1, 1)]; %#ok<AGROW>
end

for idx_interval = 1:(size(merged_interval, 1) - 1)
    gap_left = merged_interval(idx_interval, 2);
    gap_right = merged_interval(idx_interval + 1, 1);
    if gap_right > gap_left
        gap_interval = [gap_interval; gap_left, gap_right]; %#ok<AGROW>
    end
end

if merged_interval(end, 2) < two_pi
    gap_interval = [gap_interval; merged_interval(end, 2), two_pi]; %#ok<AGROW>
end
end


function split_interval = local_split_intervals_by_width(interval_mat, max_width)
% 将补集区间按给定的最大角宽做粗分箱。

split_interval = zeros(0, 2);

for idx_interval = 1:size(interval_mat, 1)
    interval_now = interval_mat(idx_interval, :);
    interval_width = interval_now(2) - interval_now(1);
    if interval_width <= 0
        continue;
    end

    num_bin = max(1, ceil(interval_width / max_width));
    edge_vec = linspace(interval_now(1), interval_now(2), num_bin + 1);
    split_interval = [split_interval; [edge_vec(1:end-1).', edge_vec(2:end).']]; %#ok<AGROW>
end
end


function interval_merge = local_merge_circular_intervals(interval_mat)
% 将一组互相连通的圆周角区间对齐到同一坐标系后再合并。
% interval_mat 的每一行对应一个待合并的原始角区间。
% 输出 interval_merge 是这些区间合并后的总角区间。

interval_aligned = zeros(size(interval_mat));
interval_aligned(1, :) = interval_mat(1, :);
merged_center = mean(interval_aligned(1, :));

for idx_interval = 2:size(interval_mat, 1)
    interval_now = interval_mat(idx_interval, :);
    center_now = mean(interval_now);
    shift_id = round((merged_center - center_now) / (2 * pi));
    interval_now = interval_now + shift_id * 2 * pi;

    shift_candidate = [-1, 0, 1] * 2 * pi;
    best_shift = 0;
    best_dist = inf;
    for idx_shift = 1:numel(shift_candidate)
        interval_try = interval_now + shift_candidate(idx_shift);
        dist_now = abs(mean(interval_try) - merged_center);
        if dist_now < best_dist
            best_dist = dist_now;
            best_shift = shift_candidate(idx_shift);
        end
    end

    interval_aligned(idx_interval, :) = interval_now + best_shift;
    merged_center = mean([min(interval_aligned(1:idx_interval, 1)), max(interval_aligned(1:idx_interval, 2))]);
end

interval_merge = local_canonicalize_interval([min(interval_aligned(:, 1)), max(interval_aligned(:, 2))]);
end


function interval_out = local_canonicalize_interval(interval_in)
% 将角区间平移到更易读的标准表示：
% 1) 如果区间本来不跨越 2pi，则尽量让它落在 [0, 2pi) 内；
% 2) 如果区间确实跨越 2pi，则保留 [a, b] 且 b > 2pi 的展开表示。

two_pi = 2 * pi;
interval_out = interval_in;

shift_num = floor(interval_out(1) / two_pi);
interval_out = interval_out - shift_num * two_pi;

while interval_out(1) < 0
    interval_out = interval_out + two_pi;
end

while interval_out(1) >= two_pi
    interval_out = interval_out - two_pi;
end
end


function theta_mean = local_circular_mean(theta_vec)
% 计算一组圆周角的均值。

theta_mean = angle(mean(exp(1i * theta_vec)));
theta_mean = local_wrapTo2Pi(theta_mean);
end


function is_in = local_is_angle_in_interval(theta_val, interval_now)
% 判断一个角度是否落入可能跨越 2pi 的区间。
% 之所以要检查 theta、theta-2pi、theta+2pi，
% 是为了兼容跨越 0/2pi 边界的圆周角区间。

theta_candidate = [theta_val, theta_val - 2 * pi, theta_val + 2 * pi];
is_in = false;
for idx_theta = 1:numel(theta_candidate)
    if theta_candidate(idx_theta) >= interval_now(1) && theta_candidate(idx_theta) <= interval_now(2)
        is_in = true;
        return;
    end
end
end


function Fake_Target_Pos = local_get_offset_fake_target(PosTarget)
% 当候选集合无法提供有效假目标时，使用一个固定偏移位置作为接口回退。

if PosTarget(1) < 0
    Fake_Target_Pos = PosTarget + [80, -100];
else
    Fake_Target_Pos = PosTarget - [80, -100];
end
end