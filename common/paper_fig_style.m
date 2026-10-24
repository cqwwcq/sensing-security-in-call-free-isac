function sty = paper_fig_style()
%PAPER_FIG_STYLE Colors, markers and export helper matching the TVT figures.
sty.blue   = [0.00 0.45 0.74];
sty.green  = [0.47 0.67 0.19];
sty.purple = [0.49 0.18 0.56];
sty.orange = [0.85 0.33 0.10];
sty.gray   = [0.40 0.40 0.40];
sty.lw = 1.2;
sty.ms = 5;
sty.font = 'Times New Roman';
sty.size_cm = [12.8 10.2];      % about the size of the published PDFs
sty.new = @new_figure;
sty.export = @export_figure;
end

function fig = new_figure()
sty = paper_fig_style();
fig = figure('Units', 'centimeters', 'Position', [2 2 sty.size_cm], 'Color', 'w');
ax = axes(fig);
hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
set(ax, 'FontName', sty.font, 'FontSize', 9, 'GridLineStyle', ':', 'GridAlpha', 0.4);
end

function export_figure(fig, file)
ax = findobj(fig, 'Type', 'axes');
set(findall(fig, '-property', 'FontName'), 'FontName', 'Times New Roman');
legend(ax(1), 'Location', 'best', 'FontSize', 8);
exportgraphics(fig, file, 'ContentType', 'vector');
fprintf('Wrote %s\n', file);
end
