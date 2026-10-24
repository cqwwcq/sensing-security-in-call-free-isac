function plot_localization_success(res_file, out_pdf)
%PLOT_LOCALIZATION_SUCCESS Fig. 2: target-localization success rate vs eta.
%
%   As in the published figure, the null-space-only baseline is plotted at
%   the proposed design's eta (paired operating points), the MI-minimizing
%   baseline at its own realized eta, and the no-AN baseline as a line.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    res_file = fullfile(here, 'detection_results.mat');
end
if nargin < 2
    out_pdf = fullfile(here, '..', ...
        'Fake-Peak Artificial Noise for Sensing-Location Privacy in Cell-Free Massive MIMO ISAC', ...
        'figure', 'localization_success_tvt.pdf');
end
O = load(res_file);
sty = paper_fig_style();

fig = sty.new();
ax = gca;
xmax = 0.9;
plot(ax, [0 xmax], 100*O.p_noan*[1 1], '--', 'Color', sty.blue, 'LineWidth', sty.lw, ...
    'DisplayName', 'No AN');
plot(ax, O.eta, 100*O.p_null, '-s', 'Color', sty.green, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Null-space-only Gaussian AN');
if isfield(O, 'p_zou')
    plot(ax, O.eta_zou, 100*O.p_zou, '-d', 'Color', sty.purple, 'LineWidth', sty.lw, ...
        'MarkerSize', sty.ms, 'DisplayName', 'MI-minimizing AN (Zou et al.)');
end
plot(ax, O.eta, 100*O.p_prop, '-o', 'Color', sty.orange, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Proposed fake-peak AN');
xlim(ax, [0 xmax]);
ylim(ax, [0 max(45, 100*max([O.p_noan, O.p_null, O.p_prop]) + 5)]);
xlabel(ax, 'Realized AN-Power Ratio \eta');
ylabel(ax, 'Target-Localization Success Rate (%)');
sty.export(fig, out_pdf);
end
