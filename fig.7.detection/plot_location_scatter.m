function plot_location_scatter(sweep_file, out_pdf)
%PLOT_LOCATION_SCATTER Response-letter figure for Reviewer 1, Comment 3.
%
%   Locations estimated by the dominant-peak attack over all channel
%   realizations, without AN (left) and with the proposed design at
%   eta ~ 0.5 (right).  Estimates are cell centers, so each marker's area is
%   proportional to the number of realizations at that cell.  Also shown:
%   the target, the virtual location, and the transmit and receive APs.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    sweep_file = fullfile(here, 'detection_sweep.mat');
end
if nargin < 2
    out_pdf = fullfile(here, 'location_scatter_response.pdf');
end
S = load(sweep_file);
cfg = S.cfg;
sty = paper_fig_style();
[~, j] = min(abs(mean(S.eta_prop, 2) - 0.5));

est = {reshape(S.est_noan, 2, []).', squeeze(S.est_prop(j, :, :)).'};
titles = {'No AN', sprintf('Proposed, \\eta = %.2f', mean(S.eta_prop(j, :)))};
colors = {sty.blue, sty.orange};

fig = sty.new();
set(fig, 'Position', [2 2 2*sty.size_cm(1) sty.size_cm(2)]);
delete(gca);
for s = 1:2
    ax = subplot(1, 2, s);
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    set(ax, 'FontName', sty.font, 'FontSize', 9, 'GridLineStyle', ':');
    [loc, ~, idx] = unique(est{s}, 'rows');
    cnt = accumarray(idx, 1);
    scatter(ax, loc(:, 1), loc(:, 2), 8 + 400*cnt/max(cnt), colors{s}, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'DisplayName', 'Estimates');
    plot(ax, cfg.PosTarget(1), cfg.PosTarget(2), 'kx', 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'DisplayName', 'Target');
    plot(ax, cfg.PosFake(1), cfg.PosFake(2), 'k+', 'MarkerSize', 10, ...
        'LineWidth', 1.5, 'DisplayName', 'Virtual location');
    plot(ax, cfg.PosTx(:, 1), cfg.PosTx(:, 2), 'k^', 'MarkerSize', 6, ...
        'DisplayName', 'Transmit APs');
    plot(ax, cfg.PosRx(:, 1), cfg.PosRx(:, 2), 'ks', 'MarkerSize', 6, ...
        'DisplayName', 'Receive APs');
    axis(ax, 'equal');
    xlim(ax, cfg.area(1:2) + [-10 10]);
    ylim(ax, cfg.area(3:4) + [-10 10]);
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    title(ax, titles{s}, 'FontWeight', 'normal');
end
sty.export(fig, out_pdf);
end
