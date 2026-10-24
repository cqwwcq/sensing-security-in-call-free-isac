function plot_multi_hypothesis(res_file, out_pdf)
%PLOT_MULTI_HYPOTHESIS Response-letter figure for Reviewer 2, Comment 5.
%
%   Localization success of the proposed design versus eta under the
%   dominant-peak attack of the paper and under the multi-hypothesis attack
%   of LOCALIZE_MULTI with its three decision rules and with two guesses.
%   Dashed lines: the same attacks without AN.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
if nargin < 1
    res_file = fullfile(here, 'detection_results.mat');
end
if nargin < 2
    out_pdf = fullfile(here, 'multi_hypothesis_response.pdf');
end
O = load(res_file);
sty = paper_fig_style();

rules = {'p',     'p_noan',     sty.orange, '-o', 'Dominant peak (paper)'
         'p_mhv', 'p_mhv_noan', sty.green,  '-s', 'Two peaks, more votes'
         'p_mhr', 'p_mhr_noan', sty.purple, '-d', 'Two peaks, random choice'
         'p_mha', 'p_mha_noan', sty.blue,   '-^', 'Two peaks, mechanism-aware'
         'p_mhs', 'p_mhs_noan', sty.gray,   '-x', 'Two guesses'};

fig = sty.new();
ax = gca;
xmax = 0.9;
for r = 1:size(rules, 1)
    y = O.([rules{r, 1} '_prop']);
    plot(ax, O.eta, 100*y, rules{r, 4}, 'Color', rules{r, 3}, 'LineWidth', sty.lw, ...
        'MarkerSize', sty.ms, 'DisplayName', rules{r, 5});
    plot(ax, [0 xmax], 100*O.(rules{r, 2})*[1 1], '--', 'Color', rules{r, 3}, ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
end
xlim(ax, [0 xmax]);
ylim(ax, [0 100]);
xlabel(ax, 'Realized AN-Power Ratio \eta');
ylabel(ax, 'Target-Localization Success Rate (%)');
sty.export(fig, out_pdf);
end
