clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Per-AP directional mutual information I(Theta ; Theta^p)
%  Bar chart: without AN vs. with the proposed (fake-peak) AN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- AP index and MI values ---
ap_idx     = 1:8;
mi_without = [3.39, 2.96, 3.40, 3.36, 3.54, 3.74, 3.69, 3.73];  % no AN
mi_with    = [0.74, 0.60, 0.83, 0.62, 1.16, 0.70, 0.82, 0.67];  % with proposed AN

% --- Plot ---
figure('Color', 'w');

b = bar(ap_idx, [mi_without; mi_with]', 'grouped', 'LineWidth', 0.6);
b(1).FaceColor = [0.12, 0.27, 0.60];   % blue  : without AN
b(2).FaceColor = [0.84, 0.39, 0.17];   % orange: with AN

hold on; box on; grid on;

% --- Axes ---
xlabel('Number of AP', 'FontSize', 11);
ylabel('Mutual Information', 'FontSize', 11);
xticks(ap_idx);
ylim([0, 4]);
yticks(0:0.5:4);

% --- Legend ---
legend({'Without AN', 'With AN'}, ...
       'Location', 'northwest', 'FontSize', 10, ...
       'Box', 'on', 'EdgeColor', [0.6 0.6 0.6]);

set(gca, 'FontSize', 10, 'LineWidth', 0.8, ...
    'GridAlpha', 0.25, 'GridLineStyle', '-', ...
    'YMinorTick', 'on', 'TickDir', 'in');
