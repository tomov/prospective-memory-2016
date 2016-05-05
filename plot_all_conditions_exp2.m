function plot_all_conditions_exp2( stats, ymin, ymax, slope, intercept, show_legend, og_only_range)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

focal_titles = {'Nonfoc', 'Focal'};
markers = {'b-*', 'r-*', 'g-*', 'm-*'};
emphasis_titles = {'Low', 'High'};
block_titles = {'Block #1', 'Block #2', 'Block #3', 'Block #4'};
og_only_titles = {'PM', 'NPM'};
handles = [];

legend_titles = {};
plot_id = 1;
for OG_ONLY = og_only_range
    for FOCAL = 0:1
        xticklabels = {};
        values = [];
        errors = [];
        xes = [];
        x = 1;
        for BLOCK = 1:4
            stat = stats(stats(:, 1) == OG_ONLY & stats(:, 2) == FOCAL & stats(:, 6) == BLOCK, :);
            M = stat(4) * slope + intercept;
            SD = stat(5) * slope;
            values = [values, M];
            errors = [errors, SD];
            xes = [xes, x];
            x = x + 1;
            xticklabel = sprintf('%s', block_titles{BLOCK});
            xticklabels = [xticklabels,  {xticklabel}];
        end
        hold on;
        errorbar(xes, values, errors);
        handle = plot(xes, values, markers{plot_id}, 'LineWidth', 2, 'MarkerSize', 6);
        plot_id = plot_id + 1;
        handles = [handles, handle];
        legend_title = sprintf('%s, %s', focal_titles{FOCAL+1}, og_only_titles{OG_ONLY+1});
        legend_titles = [legend_titles; legend_title];
    end
end

hold off;
axis([0 x ymin ymax]);
set(gca, 'XTickLabel', [{''}, xticklabels, {''}]);
if show_legend
    h = legend(handles, legend_titles);
    set(h, 'FontSize', 9);
end


h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);



