
clc; clear; close all;

%% 1. 系统参数设置
xtarget = -75;
ytarget = -xtarget;
PosTarget = [xtarget,ytarget];
Fake_Target_Pos=PosTarget + [80, -100];
% Ntx = 8; %number of transmit APs
% Mtx = 64; %number of antennas at all the APs
Ntx = 8; %number of transmit APs
Mtx = 16; %number of antennas at all the APs
Nrx = 4; %number of receive APs
Nue = 2; %number of users
Mue = [16]; %number of antennas at the users
theta_scan = -90:0.5:90; % 扫描角度范围
k_sel=7;
d_lambda=1/2;
%% 2. 构建分布式阵列流形 (Steering Vector)
steering_vec = @(theta) exp(1j * 2 * pi * d_lambda * (0:Mtx-1)' * sind(theta));

[lambdac,betark,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,Userbeamsteering,Fake_Txbeamsteering,alphark,sigmarcs,H,Z_shd]=Parameters(Ntx,Nrx,Mtx,Nue,Mue,PosTarget,Fake_Target_Pos);
% 获取特定方向的导向矢量
for k=1:Ntx
    a_target(:,k) = squeeze(Txbeamsteering(:,k,:));% Mtx x 1 row vector
    a_fake(:,k)   = squeeze(Fake_Txbeamsteering(:,k,:));
end

%% 3. 波束赋形权重计算 (模拟优化结果)

% --- 方案 A: 基准方案 (无防御 / No Defense) ---


W1 =load('D:\研究生毕业论文\童栎丞论文\童栎丞小论文英文版\sim_codes\fig4\WMC0.mat');
W1=W1.W;

% --- 方案 B: 本文提出来的方案 (Proposed with AN) ---
W2  =load('D:\研究生毕业论文\童栎丞论文\童栎丞小论文英文版\sim_codes\fig4\WMC0.mat');
Z =load('D:\研究生毕业论文\童栎丞论文\童栎丞小论文英文版\sim_codes\fig4\ZMC0.mat');
W2=W2.W;
Z=Z.Z;
%% 4. 计算波束图 (Beampattern Calculation)
P_benchmark = zeros(size(theta_scan));
P_proposed  = zeros(size(theta_scan));
itery = 1;
MC=5;
w_benchmark=W1(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC);
w_proposed=W2(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC);
Z_AN=Z(:,itery,:,MC);
w_proposed=cat(2,w_proposed,Z_AN);
for i = 1:length(theta_scan)
    a_tmp = steering_vec(theta_scan(i));
    
    % 计算阵列增益 
    P_benchmark(i) = a_tmp'*w_benchmark(:,:,k_sel)*w_benchmark(:,:,k_sel)'*a_tmp;
    P_proposed(i)  = a_tmp'*w_proposed(:,:,k_sel)*w_proposed(:,:,k_sel)'*a_tmp;
end

% 归一化为 dB
P_benchmark_dB = 10 * log10(P_benchmark / max(P_benchmark));
P_proposed_dB  = 10 * log10(P_proposed / max(P_proposed));

%% 6. 绘图
% 创建一个大的 Figure 窗口
% figure('Color', 'w', 'Position', [100, 100, 1000, 600]);

% --- 子图 1: 基准方案 (Benchmark) ---
% subplot(1, 2, 1); 
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
plot(theta_scan, P_benchmark_dB, 'Color', [0.2, 0.4, 0.6], 'LineWidth', 1.5);
grid on; box on;
xlabel('Angle (degree)', 'FontSize', 25, 'FontWeight', 'bold');
ylabel('Normalized Gain (dB)', 'FontSize', 25, 'FontWeight', 'bold');
xlim([-90, 90]); ylim([-30, 0]);
% title('(a) Conventional Design', 'FontSize', 18, 'FontWeight', 'bold');
% 可选：标注真实目标区域
xline(0, ':k'); % 示意线，实际需根据真实角度画

% --- 子图 2: 提议方案 (Proposed) ---
% subplot(1, 2, 2); 
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
plot(theta_scan, P_proposed_dB, 'Color', [0.85, 0.32, 0.1], 'LineWidth', 1.5);
grid on; box on;
xlabel('Angle (degree)', 'FontSize', 25, 'FontWeight', 'bold');
ylabel('Normalized Gain (dB)', 'FontSize', 25, 'FontWeight', 'bold');
xlim([-90,90]); ylim([-30, 0]);
% title('(b) Proposed Design', 'FontSize', 18, 'FontWeight', 'bold');

% 添加文字标注，突出"虚假波束"
% 自动寻找峰值位置进行标注 (简单演示)
% [pks, locs] = findpeaks(P_proposed_dB, theta_scan, 'MinPeakHeight', -5, 'MinPeakDistance', 20);
% for k = 1:length(pks)
%     text(locs(k), pks(k)+2, sprintf('%.1f^o', locs(k)), ...
%         'HorizontalAlignment', 'center', 'FontSize', 9);
% end