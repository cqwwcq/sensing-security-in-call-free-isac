function plot_detection_pd(res_file, out_pdf)
%PLOT_DETECTION_PD Response-letter figure for Reviewer 3, Comment 2.
%
%   Detection probability of the GLRT at Pfa = cfg.Pfa versus an RCS
%   reduction relative to Table I, for the no-AN baseline and the proposed
%   design at eta ~ 0.5.  Markers: Monte Carlo over Swerling-I echo draws on
%   the designed beams (DETECTION_PD).  Lines: the closed-form model of the
%   response letter, with exponentially distributed per-receive-AP SNRs of
%   mean gamma_t.  Agreement validates the numbers quoted in the letter.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    res_file = fullfile(here, 'detection_results.mat');
end
if nargin < 2
    out_pdf = fullfile(here, 'detection_pd_response.pdf');
end
O = load(res_file);
cfg = O.cfg;
sty = paper_fig_style();
[~, j] = min(abs(O.eta - 0.5));
off = cfg.pd_offset_dB;

Ns = cfg.Nue + 1;
df = 2 * cfg.Nrx * Ns;
tau = chi2inv(1 - cfg.Pfa, df);
G = sum(-log(rand(2e5, cfg.Nrx)), 2);        % sum of Nrx unit exponentials
pd_model = @(g_dB, N) arrayfun(@(o) ...
    mean(ncx2cdf(tau, df, 2*N*10^((g_dB + o)/10)*G, 'upper')), off);

fig = sty.new();
ax = gca;
mk = {'o', 's'};
for b = 1:numel(cfg.pd_N_list)
    N = cfg.pd_N_list(b);
    plot(ax, off, O.pd_noan_mean(:, b), mk{b}, 'Color', sty.blue, ...
        'MarkerSize', sty.ms, 'DisplayName', sprintf('No AN, N = %d', N));
    plot(ax, off, pd_model(O.snr_noan_mean, N), '-', 'Color', sty.blue, ...
        'LineWidth', sty.lw, 'HandleVisibility', 'off');
    plot(ax, off, squeeze(O.pd_prop_mean(j, :, b)), mk{b}, 'Color', sty.orange, ...
        'MarkerSize', sty.ms, 'DisplayName', ...
        sprintf('Proposed (\\eta = %.2f), N = %d', O.eta(j), N));
    plot(ax, off, pd_model(O.snr_prop_mean(j), N), '-', 'Color', sty.orange, ...
        'LineWidth', sty.lw, 'HandleVisibility', 'off');
end
ylim(ax, [0 1]);
xlabel(ax, 'RCS Change Relative to Table I (dB)');
ylabel(ax, sprintf('Detection Probability (P_{fa} = 10^{%d})', round(log10(cfg.Pfa))));
sty.export(fig, out_pdf);
end
