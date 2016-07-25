%load('exp1-data.mat')

B2010_with_stats

figure;

subplot(4, 2, 1);
title('Empirical Data');
ylabel('PM Hit Rate (%)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 10 11]), 00, 100, 1, 0, true);

subplot(4, 2, 2);
title('Simulation Data');
plot_all_conditions_brewer(simulation_stats(:, [1:3 10 11]), 00, 100, 1, 0, false);

text(-10,10.2,'Test title spanning two subplots -- Some fine tuning will be required')

subplot(4, 2, 3);
ylabel('OG RT (msec)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 4 5]), 650, 1850, 1, 0, false);

subplot(4, 2, 4);
plot_all_conditions_brewer(simulation_stats(:, [1:3 4 5]), 650, 1850, RT_slope, RT_intercept, false);


subplot(4, 2, 5);
ylabel('1st PM RT (msec)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 13 14]), 600, 2500, 1, 0, false);

subplot(4, 2, 6);
plot_all_conditions_brewer(simulation_stats(:, [1:3 8 9]), 600, 2500, RT_slope, RT_intercept, false);

%
% For sanity check
%

subplot(4, 2, 7);
ylabel('OG Accuracy (%)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 6 7]), 00, 100, 1, 0, false);

subplot(4, 2, 8);
plot_all_conditions_brewer(simulation_stats(:, [1:3 6 7]), 00, 100, 1, 0, false);


