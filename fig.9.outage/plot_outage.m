function plot_outage(res_file, out_pdf)
%PLOT_OUTAGE Fig. 4: communication outage probability vs CSI NMSE.
%
%   Plots the four published variants and, if present in the results file,
%   the robust no-AN baseline (Reviewer 1, Comment 5).  Also prints the
%   values quoted in the paper.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    res_file = fullfile(here, 'outage_results.mat');
end
if nargin < 2
    out_pdf = fullfile(here, '..', ...
        'Fake-Peak Artificial Noise for Sensing-Location Privacy in Cell-Free Massive MIMO ISAC', ...
        'figure', 'outage_vs_nmse.pdf');
end
O = load(res_file);
sty = paper_fig_style();

fig = sty.new();
ax = gca;
x = O.nmse;
plot(ax, x, 100*O.m_noan, '-+', 'Color', sty.green, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'No AN');
plot(ax, x, 100*O.m_null, '-d', 'Color', sty.purple, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Null-space-only Gaussian AN');
plot(ax, x, 100*O.m_nom, '-^', 'Color', sty.orange, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Proposed nominal');
plot(ax, x, 100*O.m_rob, '-s', 'Color', sty.blue, 'LineWidth', sty.lw, ...
    'MarkerSize', sty.ms, 'DisplayName', 'Proposed robust');
if isfield(O, 'm_rnoan')
    plot(ax, x, 100*O.m_rnoan, '--x', 'Color', sty.gray, 'LineWidth', sty.lw, ...
        'MarkerSize', sty.ms, 'DisplayName', 'Robust no AN');
end
xlim(ax, [min(x) max(x)]);
ylim(ax, [0 100]);
xticks(ax, x(1):4:x(end));
xlabel(ax, 'CSI NMSE (dB)');
ylabel(ax, 'Communication Outage Probability (%)');
sty.export(fig, out_pdf);

for v = [-10 -4]
    j = find(x == v, 1);
    if isempty(j)
        continue;
    end
    fprintf('NMSE %d dB: outage no AN %.2f %%, null-only %.2f %%, nominal %.2f %%, robust %.2f %%', ...
        v, 100*O.m_noan(j), 100*O.m_null(j), 100*O.m_nom(j), 100*O.m_rob(j));
    if isfield(O, 'm_rnoan')
        fprintf(', robust no AN %.2f %%', 100*O.m_rnoan(j));
    end
    fprintf('\n  sensing SNR no AN %.2f dB, nominal %.2f dB, robust %.2f dB', ...
        O.m_sna(j), O.m_sn(j), O.m_sr(j));
    if isfield(O, 'm_srn')
        fprintf(', robust no AN %.2f dB', O.m_srn(j));
    end
    fprintf('\n');
end
end
