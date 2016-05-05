function plot_all_conditions_exp3( stats, ymin, ymax, slope, intercept, show_legend)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

target_titles = {'1 Target', '6 Targets'};
markers = {'b-o', 'r-*'};
emphasis_titles = {'Low', 'High'};
og_only_titles = {'PM', 'No PM'};
handles = [];

legend_titles = {};
target_idx = 0;
for TARGETS = [1,6]
    target_idx = target_idx + 1;
    xticklabels = {};
    values = [];
    errors = [];
    xes = [];
    x = 1;
    for OG_ONLY = 1:-1:0
        for EMPHASIS = 0
            stat = stats(stats(:, 1) == OG_ONLY & stats(:, 2) == TARGETS & stats(:, 3) == EMPHASIS, :);
            M = stat(4) * slope + intercept;
            SD = stat(5) * slope;
            values = [values, M];
            errors = [errors, SD];
            xes = [xes, x];
            x = x + 1;
            xlabel = sprintf('%s', og_only_titles{OG_ONLY+1});
            xticklabels = [xticklabels,  {xlabel}, {''}];
        end
    end
    hold on;
    errorbar(xes, values, errors);
    handle = plot(xes, values, markers{target_idx}, 'LineWidth', 2, 'MarkerSize', 6);
    handles = [handles, handle];
    legend_titles = [legend_titles; target_titles(target_idx)];
end

hold off;
axis([0.5 x-0.5 ymin ymax]);
set(gca, 'XTickLabel', [{''}, xticklabels, {''}]);
if show_legend
    h = legend(handles, legend_titles);
    set(h, 'FontSize', 15);
end


h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);




