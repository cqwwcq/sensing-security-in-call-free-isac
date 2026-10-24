function plot_sensing_snr(res_file, out_pdf)
%PLOT_SENSING_SNR Fig. 3: target sensing SNR vs eta at the paired points.
%
%   Adds a dashed curve for the proposed design with the known AN echo
%   exploited at the receivers, gamma_t + sum_k c_k T_k (Reviewer 2,
%   Comment 1).  The solid curve is (2), with the AN echo discarded.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    res_file = fullfile(here, 'detection_results.mat');
end
if nargin < 2
    out_pdf = fullfile(here, '..', ...
        'Fake-Peak Artificial Noise for Sensing-Location Privacy in Cell-Free Massive MIMO ISAC', ...
        'figure', 'sensing_snr_final.pdf');
end
O = load(res_file);
sty = paper_fig_style();

fig = sty.new();
ax = gca;
xmax = 0.9;
plot(ax, [0 xmax], O.snr_noan_mean*[1 1], '--', 'Color', sty.blue, 'LineWidth', sty.lw, ...
    'DisplayName', 'No AN (reference)');
plot(ax, O.eta, O.snr_null_mean, '-s', 'Color', sty.green, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Null-space-only Gaussian AN');
if isfield(O, 'snr_zou_mean')
    plot(ax, O.eta_zou, O.snr_zou_mean, '-d', 'Color', sty.purple, 'LineWidth', sty.lw, ...
        'MarkerSize', sty.ms, 'DisplayName', 'MI-minimizing AN (Zou et al.)');
end
plot(ax, O.eta, O.snr_prop_mean, '-o', 'Color', sty.orange, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Proposed fake-peak AN');
plot(ax, O.eta, O.snrz_prop_mean, '--o', 'Color', sty.orange, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'MarkerFaceColor', 'w', ...
    'DisplayName', 'Proposed, AN echo exploited');
xlim(ax, [0 xmax]);
xlabel(ax, 'Realized AN-Power Ratio \eta');
ylabel(ax, 'Target Sensing SNR (dB)');
sty.export(fig, out_pdf);
end
