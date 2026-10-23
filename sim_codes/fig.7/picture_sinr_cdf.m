clear; clc; close all;

sinr_noAN = load('sinr_noAN.mat');
sinr_gaussAN = load('sinr_gaussAN.mat');
sinr_peakAN = load('sinr_peakAN.mat');

sinr_noAN = sinr_noAN.min_sinr;
sinr_gaussAN = sinr_gaussAN.min_sinr;
sinr_peakAN = sinr_peakAN.min_sinr;

% CDF plot in the same IEEE-style format as Fig. 5.
figure('Color', 'w');
hold on; box on; grid on;

data_list = {sinr_noAN, sinr_gaussAN, sinr_peakAN};
clr_noAN = [0.15, 0.15, 0.15];
clr_gauss = [0.12, 0.27, 0.60];
clr_fake = [0.84, 0.15, 0.16];
color_list = {clr_noAN, clr_gauss, clr_fake};
style_list = {'-', '--', '-'};
label_list = {'Without AN', 'Gaussian AN', 'Fake-peak AN'};

for k = 1:3
    raw = data_list{k}(:);
    raw = raw(isfinite(raw));

    x_sorted = sort(raw);
    f_sorted = (1:length(x_sorted))' / length(x_sorted);

    [x_unique, ia, ~] = unique(x_sorted, 'last');
    f_unique = f_sorted(ia);

    xx = linspace(min(x_unique), max(x_unique), 500);
    ff = interp1(x_unique, f_unique, xx, 'pchip');
    ff = max(0, min(1, ff));

    plot(xx, ff, style_list{k}, ...
        'Color', color_list{k}, ...
        'LineWidth', 1.6, ...
        'DisplayName', label_list{k});
end

legend('Location', 'southeast', 'FontSize', 9, 'Box', 'on', ...
       'EdgeColor', [0.6 0.6 0.6]);

xlabel('Communication SINR (dB)', 'FontSize', 11);
ylabel('CDF', 'FontSize', 11);

ylim([0 1]);
set(gca, 'FontSize', 10, 'LineWidth', 0.8, ...
    'GridAlpha', 0.25, 'GridLineStyle', '-', ...
    'TickDir', 'in', 'TickLength', [0.012 0.025], ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
