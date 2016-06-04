%load('exp3-data-newww.mat')

EM2005_with_stats_exp3


figure;

subplot(3, 2, 1);
title('Empirical Data');
ylabel('OG RT (msec)');
plot_all_conditions_exp3(empirical_stats(:, [1 12 3 4 5]), 4500, 5500, 1, 0, false);

subplot(3, 2, 2);
title('Simulation Data');
%ylabel(OG_RT_label_cycles_to_msec);
plot_all_conditions_exp3(simulation_stats(:, [1 12 3 4 5]), 4500, 5500, RT_slope, RT_intercept, true);

subplot(3, 2, 3);
ylabel('OG Accuracy (%)');
plot_all_conditions_exp3(empirical_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

subplot(3, 2, 4);
plot_all_conditions_exp3(simulation_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

subplot(3, 2, 5);
ylabel('PM Hit Rate (%)');
plot_all_conditions_exp3(empirical_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);

subplot(3, 2, 6);
plot_all_conditions_exp3(simulation_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);
