clear; clc; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%gaussANVSfakeAN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% --- 伪峰 AN，MTX=16 ---
power_pseudo_16 = [0.6308, 0.6144, 0.5984, 0.5827, 0.5582, 0.5351, 0.5072, 0.4795, ...
                0.4475, 0.4065, 0.3665, 0.3163, 0.2565, 0.1879, 0.1039, 0];

det_pseudo_16 = [0.32, 0.47, 0.67, 0.95, 1.57, 2.37, 3.52, 4.64, 5.62, ...
              6.32, 6.62, 6.75, 6.79, 6.80, 6.80, 6.80] / 10;

% --- 高斯 AN，MTX=16 ---
%
power_gauss_16 = [0.88, 0.85, 0.80, 0.75, 0.70, 0.63, ...
                 0.5984, 0.5827, 0.5582, 0.5351, 0.5072, ...
                 0.4795, 0.4475, 0.4065, 0.3665, 0.3163,0.2565,0.1879,0.1039,0];

det_gauss_16 = [0.82, 1.19, 2.06, 3.20, 4.41, 5.71, ...
             6.09, 6.24, 6.42, 6.55, 6.66, ...
             6.74, 6.80, 6.84, 6.87, 6.88, 6.89, 6.90, 6.90, 6.90
] / 10;
% --- 伪峰 AN，MTX=8 ---
power_pseudo_8 = [0.6748,0.6644,0.6444,0.6252,0.6064,0.5811,0.5588,0.5322,...
0.4976,0.4612,0.4129,0.3667,0.3049,0.2234,0.1235,0.0019];

det_pseudo_8 = [0.60, 0.67, 0.83, 1.02, 1.22, 1.54, 1.87, ...
                2.28, 2.85, 3.44, 4.09, 4.54, 4.92, 5.16, 5.26, 5.29] / 10;

% --- 高斯 AN，MTX=8 ---
%
power_gauss_8 = [0.9,0.85,0.8,0.75,0.7,0.6748,0.6644,0.6444,...
                0.6252,0.6064,0.5811,0.5588,0.5322,0.4976,0.4612,0.4129,0.3667,0.3049,0.2234,0.1235,0.0019
];

det_gauss_8 = [0.10, 0.21, 0.47, 0.97, 1.80, 2.33, 2.55, 2.99, ...
              3.40, 3.76, 4.18, 4.47, 4.74, 4.97, 5.12, 5.22, 5.26, 5.29, 5.30, 5.30, 5.30
] / 10;


% --- 翻转使横坐标从小到大 ---
power_pseudo_16 = fliplr(power_pseudo_16);
det_pseudo_16   = fliplr(det_pseudo_16);
power_gauss_16  = fliplr(power_gauss_16);
det_gauss_16    = fliplr(det_gauss_16);

power_pseudo_8 = fliplr(power_pseudo_8);
det_pseudo_8   = fliplr(det_pseudo_8);
power_gauss_8  = fliplr(power_gauss_8);
det_gauss_8    = fliplr(det_gauss_8);

% --- 平滑 ---
det_pseudo_smooth_16 = smoothdata(det_pseudo_16, 'gaussian', 4);
det_gauss_smooth_16  = smoothdata(det_gauss_16,  'gaussian', 4);
det_pseudo_smooth_8 = smoothdata(det_pseudo_8, 'gaussian', 4);
det_gauss_smooth_8  = smoothdata(det_gauss_8,  'gaussian', 4);
% --- 绘图 (IEEE-style) ---
figure('Color', 'w');

hold on; box on; grid on;

% 颜色按方案区分：红 = Fake-peak (proposed)；蓝 = Gaussian (baseline)
% 线型按阵列大小区分：实线 = Mtx=16；虚线 = Mtx=8
clr_fake  = [0.84, 0.15, 0.16];
clr_gauss = [0.12, 0.27, 0.60];

% Fake-peak AN, Mtx=16: red solid, filled circle
plot(power_pseudo_16 * 100, det_pseudo_smooth_16 * 100, '-o', ...
    'Color', clr_fake, 'LineWidth', 1.6, ...
    'MarkerSize', 6, 'MarkerFaceColor', clr_fake, ...
    'DisplayName', 'Fake-peak AN (M_{tx}=16)');

% Gaussian AN, Mtx=16: blue solid, filled square
plot(power_gauss_16 * 100, det_gauss_smooth_16 * 100, '-s', ...
    'Color', clr_gauss, 'LineWidth', 1.6, ...
    'MarkerSize', 6, 'MarkerFaceColor', clr_gauss, ...
    'DisplayName', 'Gaussian AN (M_{tx}=16)');

% Fake-peak AN, Mtx=8: red dashed, hollow triangle
plot(power_pseudo_8 * 100, det_pseudo_smooth_8 * 100, '--^', ...
    'Color', clr_fake, 'LineWidth', 1.6, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
    'DisplayName', 'Fake-peak AN (M_{tx}=8)');

% Gaussian AN, Mtx=8: blue dashed, hollow diamond
plot(power_gauss_8 * 100, det_gauss_smooth_8 * 100, '--d', ...
    'Color', clr_gauss, 'LineWidth', 1.6, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'w', ...
    'DisplayName', 'Gaussian AN (M_{tx}=8)');

% --- 坐标轴 ---
xlabel('AN Power Ratio \eta (%)', 'FontSize', 11);
ylabel('Detection Rate (%)', 'FontSize', 11);

xlim([0, 90]);
ylim([0, 75]);

% --- 图例 ---
legend('Location', 'southwest', 'FontSize', 9, 'Box', 'on', ...
       'EdgeColor', [0.6 0.6 0.6]);

set(gca, 'FontSize', 10, 'LineWidth', 0.8, ...
    'GridAlpha', 0.25, 'GridLineStyle', '-', ...
    'TickDir', 'in', 'TickLength', [0.012 0.025], ...
    'XMinorTick', 'on', 'YMinorTick', 'on');