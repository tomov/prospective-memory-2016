%load('exp5-data.mat')

EM2005_with_stats_exp5

figure;

% plot aftereffects of intention
%

subplot(2, 2, 1);

empirical_Ms = empirical_stats(2:3, 4)';
empirical_SEMs = empirical_stats(2:3, 5)';
barweb(empirical_Ms, empirical_SEMs, 1, {}, ...
    'Empirical Data', 'Third Task Trial Type', 'Third Task RT (msec)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([500 700]);



subplot(2, 2, 2);

simulation_Ms = simulation_stats(2:3, 4)' * RT_slope + RT_intercept;
simulation_SEMs = simulation_stats(2:3, 5)' * RT_slope;
%
%simulation_Ms = [mean(subjects(:, 14))    mean(subjects(:, 12)) ] * RT_slope + RT_intercept;
%simulation_SEMs = [std(subjects(:, 14))    std(subjects(:, 12)) ] * RT_slope / sqrt(simulation_subjects_per_condition);


barweb(simulation_Ms, simulation_SEMs, 1, {}, ...
    'Simulation Data', 'Third Task Trial Type', 'Third Task RT (msec)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([1070 1170]);



% Plot the accuracies
%


subplot(2, 2, 3);

empirical_Ms = [96 96]; % don't think E&M report a difference here
empirical_SEMs = [3 3] / sqrt(empirical_subjects_per_condition);
barweb(empirical_Ms, empirical_SEMs, 1, {}, ...
    '', 'Third Task Trial Type', 'Third Task Accuracy (%)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([0 100]);



subplot(2, 2, 4);

simulation_Ms = [mean(subjects(:, 15)) mean(subjects(:, 13))];
simulation_SEMs = [std(subjects(:, 15)) std(subjects(:, 13))] / sqrt(simulation_subjects_per_condition);
barweb(simulation_Ms, simulation_SEMs, 1, {}, ...
    '', 'Third Task Trial Type', 'Third Task Accuracy (%)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([0 100]);
