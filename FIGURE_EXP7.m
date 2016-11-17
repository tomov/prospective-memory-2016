subjects_per_condition = 100;
blocks = 20; % a "block" is a few OG + PM trials, followed by a few third task trials (with one target item), followed by a few OG + PM trials

% from experiment 1
RT_slope = 8.9868;
RT_intercept = 360;

IT_target_RTs = extra(:, 17:17 + blocks - 1) * RT_slope + RT_intercept;
IT_nontarget_RTs = extra(:, 17 + blocks:17 + 2 * blocks - 1) * RT_slope + RT_intercept;
IT_tar_hits = 100 * extra(:, 17 + 2 * blocks:17 + 3 * blocks - 1);
IT_tar_pm_hits = 100 * extra(:, 17 + 3 * blocks:17 + 4 * blocks - 1); % commission errors
IT_nontar_hits = 100 * extra(:, 17 + 4 * blocks:17 + 5 * blocks - 1);


figure;

%
% plot RTs
%

subplot(2, 1, 1);

title('Aftereffects of Intention Decay');
xlabel('Third Task Block #');
ylabel('RT (msec)');

handles = [];

hold on;

% plot nontargets
Ms = mean(IT_nontarget_RTs);
SEMs = std(IT_nontarget_RTs) / sqrt(subjects_per_condition);
errorbar(1:blocks, Ms, SEMs);
handle = plot(1:blocks, Ms, 'b-*', 'LineWidth', 2, 'MarkerSize', 6);
handles = [handles, handle];
axis([0 21 1000 1600]);

% plot targets
Ms = mean(IT_target_RTs);
SEMs = std(IT_target_RTs) / sqrt(subjects_per_condition);
errorbar(1:blocks, Ms, SEMs);
handle = plot(1:blocks, Ms, 'g-*', 'LineWidth', 2, 'MarkerSize', 6);
handles = [handles, handle];
axis([0 21 1000 1600]);

hold off;
h = legend(handles, {'Non-targets', 'Targets'});
set(h, 'FontSize', 13);


h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);


%
% plot accuracies
%

subplot(2, 1, 2);

xlabel('Third Task Block #');
ylabel('Accuracy (%)');

handles = [];

hold on;

% plot nontargets
Ms = mean(IT_nontar_hits);
SEMs = std(IT_nontar_hits) / sqrt(subjects_per_condition);
errorbar(1:blocks, Ms, SEMs);
handle = plot(1:blocks, Ms, 'b-*', 'LineWidth', 2, 'MarkerSize', 6);
handles = [handles, handle];

% plot targets
Ms = mean(IT_tar_hits);
SEMs = std(IT_tar_hits) / sqrt(subjects_per_condition);
errorbar(1:blocks, Ms, SEMs);
handle = plot(1:blocks, Ms, 'g-*', 'LineWidth', 2, 'MarkerSize', 6);
handles = [handles, handle];
axis([0 21 0 100]);

% plot target misses but PM hits
Ms = nanmean(IT_tar_pm_hits);
SEMs = nanstd(IT_tar_pm_hits) / sqrt(subjects_per_condition);
errorbar(1:blocks, Ms, SEMs);
handle = plot(1:blocks, Ms, 'r-*', 'LineWidth', 2, 'MarkerSize', 6);
handles = [handles, handle];
axis([0 21 0 100]);

hold off;
h = legend(handles, {'Non-targets', 'Targets', 'Target Miss, PM Hit'});
set(h, 'FontSize', 13);


h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);

