%load('exp2-data.mat')

EM2005_with_stats_exp2


fit_as_experiment_1 = false;
if fit_as_experiment_1
    RT_slope = 10;
    RT_intercept = 205;
end

figure;

subplot(3, 2, 1);
title('Empirical Data');
ylabel('OG RT (msec)');
plot_all_conditions_exp2(empirical_stats(:, [1:3 4 5 12]), 900, 1300, 1, 0, false, [0 1]);

subplot(3, 2, 2);
title('Simulation Data');
plot_all_conditions_exp2(simulation_stats(:, [1:3 4 5 12]), 900, 1300, RT_slope, RT_intercept, true, [0 1]);

subplot(3, 2, 3);
ylabel('OG Accuracy (%)');
plot_all_conditions_exp2(empirical_stats(:, [1:3 6 7 12]), 20, 100, 1, 0, false, [0 1]);

subplot(3, 2, 4);
plot_all_conditions_exp2(simulation_stats(:, [1:3 6 7 12]), 20, 100, 1, 0, false, [0 1]);

subplot(3, 2, 5);
ylabel('PM Hit Rate (%)');
plot_all_conditions_exp2(empirical_stats(:, [1:3 10 11 12]), 0, 100, 1, 0, false, [0]);

subplot(3, 2, 6);
plot_all_conditions_exp2(simulation_stats(:, [1:3 10 11 12]), 0, 100, 1, 0, false, [0]);
