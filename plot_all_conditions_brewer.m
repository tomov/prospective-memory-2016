function plot_all_conditions_exp1( stats, ymin, ymax, slope, intercept, show_legend)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

xticklabels = {'Baseline', 'Focal', 'Nonfocal'};
capacity_titles = {'Low WM', 'High WM'};
markers = {'r-*', 'b-*'};
handles = [];

legend_titles = {};
for CAPACITY = 0:1
    values = [];
    errors = [];
    xes = [];
    x = 1;
    for OG_ONLY = 1:-1:0
        for FOCAL = 1:-1:0
            if OG_ONLY && ~FOCAL
                continue % no such data point
            end            
            stat = stats(stats(:, 1) == OG_ONLY & stats(:, 2) == FOCAL & stats(:, 3) == CAPACITY, :);
            M = stat(4) * slope + intercept;
            SD = stat(5) * slope;
            values = [values, M];
            errors = [errors, SD];
            xes = [xes, x];
            x = x + 1;
        end
    end
    hold on;
    errorbar(xes, values, errors);
    handle = plot(xes, values, markers{CAPACITY+1}, 'LineWidth', 2, 'MarkerSize', 6);
    handles = [handles, handle];
    legend_titles = [legend_titles; capacity_titles(CAPACITY+1)];
end

hold off;
axis([0.5 x-0.5 ymin ymax]);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', xticklabels);
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




