%load('exp1-data-newww.mat')

EM2005_with_stats_exp1

figure;

subplot(4, 2, 1);
title('Empirical Data');
ylabel('OG RT (msec)');
plot_all_conditions_exp1(empirical_stats(:, [1:3 4 5]), 1000, 1700, 1, 0, true);

subplot(4, 2, 2);
title('Simulation Data');
%ylabel(OG_RT_label_cycles_to_msec);
plot_all_conditions_exp1(simulation_stats(:, [1:3 4 5]), 1000, 1700, RT_slope, RT_intercept, false);

subplot(4, 2, 3);
ylabel('OG Accuracy (%)');
plot_all_conditions_exp1(empirical_stats(:, [1:3 6 7]), 0, 100, 1, 0, false);

subplot(4, 2, 4);
plot_all_conditions_exp1(simulation_stats(:, [1:3 6 7]), 0, 100, 1, 0, false);

subplot(4, 2, 5);
ylabel('PM Hit Rate (%)');
plot_all_conditions_exp1(empirical_stats(:, [1:3 10 11]), 0, 100, 1, 0, false);

subplot(4, 2, 6);
plot_all_conditions_exp1(simulation_stats(:, [1:3 10 11]), 0, 100, 1, 0, false);

%
% For sanity
%

subplot(4, 2, 7);
ylabel('PM RT (msec)');
plot_all_conditions_exp1(empirical_stats(:, [1:3 8 9]), 1000, 2300, 1, 0, false);

subplot(4, 2, 8);
plot_all_conditions_exp1(simulation_stats(:, [1:3 8 9]), 1000, 2300, RT_slope, RT_intercept, false);

