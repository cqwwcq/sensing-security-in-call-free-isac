clear; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%坐标位置分布图绘制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. 数据定义
PosAPs = [-200,-200;200,200;-200,200;200,-200;0,-200;0,200;-200,0;200,0];


Posue = [50,50;-50 -50];

PosRx = [150,150;-150,150;150,-150;-150,-150]; %location of Receive APs for Nrx= 4

xtarget = -75;
ytarget = -xtarget; % 即 75
PosTarget = [xtarget, ytarget];
Fake_Target_Pos=PosTarget + [80, -100];
%% 2. 绘图设置
figure('Color', 'w', 'Name', 'System Topology'); 
hold on; 
grid on; 
box on;
axis equal; % 保持横纵坐标比例一致，防止图形变形

%% 3. 绘制各个元素

% --- 绘制发射 AP (Tx APs) ---
% 使用蓝色正方形表示
h_ap = scatter(PosAPs(:,1), PosAPs(:,2), 80, 's', ...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'c', ...
    'LineWidth', 1.5, 'DisplayName', 'Tx APs');

% *** 关键步骤：为 AP 添加编号标签 ***
for i = 1:size(PosAPs, 1)
    % 在点旁边添加文本，(i) 就是输入的顺序索引
    text(PosAPs(i,1) + 5, PosAPs(i,2) + 5, sprintf('AP-%d', i), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
end

% --- 绘制接收 AP (Rx APs) ---
% 使用绿色菱形表示
h_rx = scatter(PosRx(:,1), PosRx(:,2), 80, 'd', ...
    'MarkerEdgeColor', 'g', 'MarkerFaceColor', [0.8 1 0.8], ...
    'LineWidth', 1.5, 'DisplayName', 'Rx APs');

% (可选) 为 Rx AP 添加标签
for i = 1:size(PosRx, 1)
    text(PosRx(i,1) + 3, PosRx(i,2) - 5, sprintf('Rx-%d', i), ...
        'FontSize', 8, 'Color', [0 0.5 0]);
end

% --- 绘制用户 (Users) ---
% 使用红色圆圈表示
h_ue = scatter(Posue(:,1), Posue(:,2), 80, 'o', ...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', [1 0.8 0.8], ...
    'LineWidth', 1.5, 'DisplayName', 'Users (UE)');

% (可选) 为用户添加标签
for i = 1:size(Posue, 1)
    text(Posue(i,1) + 3, Posue(i,2), sprintf('UE-%d', i), ...
        'FontSize', 9, 'Color', 'r');
end

% --- 绘制目标 (Target) ---
% 使用黑色五角星表示
h_tg = scatter(PosTarget(1), PosTarget(2), 150, 'p', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', ...
    'LineWidth', 1.5, 'DisplayName', 'Target');

text(PosTarget(1) + 3, PosTarget(2), 'Target', ...
    'FontSize', 10, 'FontWeight', 'bold');

% --- 绘制假目标 (Fake Target) ---
% 使用紫色五角星 / 或叉号区分
h_fake = scatter(Fake_Target_Pos(1), Fake_Target_Pos(2), 150, '^', ...
    'MarkerEdgeColor', [0.5 0 0.5], ...
    'MarkerFaceColor', [0.9 0.7 0.9], ...
    'LineWidth', 1.5, 'DisplayName', 'Fake Target');

text(Fake_Target_Pos(1) + 3, Fake_Target_Pos(2), 'Fake Target', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.5 0 0.5]);
%% 4. 图形美化
xlabel('X Position (m)', 'FontSize', 12);
ylabel('Y Position (m)', 'FontSize', 12);
title('CF-mMIMO-ISAC System Layout', 'FontSize', 14);

% 设置图例
legend([h_ap, h_rx, h_ue, h_tg], 'Location', 'bestoutside', 'FontSize', 10);

% 自动调整坐标轴范围，留出一点边距
margin = 20;
xlim([min([PosAPs(:,1); Posue(:,1); PosRx(:,1); xtarget]) - margin, ...
      max([PosAPs(:,1); Posue(:,1); PosRx(:,1); xtarget]) + margin]);
ylim([min([PosAPs(:,2); Posue(:,2); PosRx(:,2); ytarget]) - margin, ...
      max([PosAPs(:,2); Posue(:,2); PosRx(:,2); ytarget]) + margin]);

