clear; clc; close all;

data = load('final_Poscont.mat');
final_Poscont = data.final_Poscont;

PosTarget = [-75, 75];
MC = 1000;
Nseed = 10;

err_noAN = zeros(MC * Nseed, 1);
err_gaussAN = zeros(MC * Nseed, 1);
err_AN = zeros(MC * Nseed, 1);
err_AN_210 = zeros(MC * Nseed, 1);
err_AN_269 = zeros(MC * Nseed, 1);

idx = 1;
for s = 1:Nseed
    for mc = 1:MC
        err_gaussAN(idx) = norm(final_Poscont(mc, :, 1, s) - PosTarget);
        err_AN(idx) = norm(final_Poscont(mc, :, 2, s) - PosTarget);
        err_noAN(idx) = norm(final_Poscont(mc, :, 3, s) - PosTarget);
        err_AN_210(idx) = norm(final_Poscont(mc, :, 4, s) - PosTarget);
        err_AN_269(idx) = norm(final_Poscont(mc, :, 5, s) - PosTarget);
        idx = idx + 1;
    end
end

% CDF plot in the same IEEE-style format as Fig. 5.
figure('Color', 'w');
hold on; box on; grid on;

data_list = {err_noAN, err_gaussAN, err_AN, err_AN_210, err_AN_269};
clr_noAN = [0.15, 0.15, 0.15];
clr_gauss = [0.12, 0.27, 0.60];
clr_fake1 = [0.84, 0.15, 0.16];
clr_fake2 = [0.93, 0.55, 0.13];
clr_fake3 = [0.45, 0.28, 0.70];
color_list = {clr_noAN, clr_gauss, clr_fake1, clr_fake2, clr_fake3};
style_list = {'-', '--', '-', '-.', ':'};
label_list = {'Without AN', 'Gaussian AN', 'Fake-peak AN 1', ...
              'Fake-peak AN 2', 'Fake-peak AN 3'};

for k = 1:numel(data_list)
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

xlabel('Localization Error (m)', 'FontSize', 11);
ylabel('CDF', 'FontSize', 11);

ylim([0 1]);
set(gca, 'FontSize', 10, 'LineWidth', 0.8, ...
    'GridAlpha', 0.25, 'GridLineStyle', '-', ...
    'TickDir', 'in', 'TickLength', [0.012 0.025], ...
    'XMinorTick', 'on', 'YMinorTick', 'on');
