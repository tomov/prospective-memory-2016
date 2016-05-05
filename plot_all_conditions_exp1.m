function plot_all_conditions_exp1( stats, ymin, ymax, slope, intercept, show_legend)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

focal_titles = {'Nonfocal', 'Focal'};
markers = {'b-*', 'r-*'};
emphasis_titles = {'Lo', 'Hi'};
og_only_titles = {'PM', 'NPM'};
handles = [];

legend_titles = {};
for FOCAL = 1:-1:0
    xticklabels = {};
    values = [];
    errors = [];
    xes = [];
    x = 1;
    for OG_ONLY = 1:-1:0
        for EMPHASIS = 0:1
            stat = stats(stats(:, 1) == OG_ONLY & stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, :);
            M = stat(4) * slope + intercept;
            SD = stat(5) * slope;
            values = [values, M];
            errors = [errors, SD];
            xes = [xes, x];
            x = x + 1;
            xlabel = sprintf('%s,%s', emphasis_titles{EMPHASIS+1}, og_only_titles{OG_ONLY+1});
            xticklabels = [xticklabels,  {xlabel}];
        end
    end
    hold on;
    errorbar(xes, values, errors);
    handle = plot(xes, values, markers{FOCAL+1}, 'LineWidth', 2, 'MarkerSize', 6);
    handles = [handles, handle];
    legend_titles = [legend_titles; focal_titles(FOCAL+1)];
end

hold off;
axis([0.5 x-0.5 ymin ymax]);
set(gca, 'XTickLabel', [xticklabels, {''}]);
if show_legend
    h = legend(handles, legend_titles);
    set(h, 'FontSize', 13);
end

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);




