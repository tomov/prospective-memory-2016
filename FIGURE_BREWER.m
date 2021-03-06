%load('exp1-data.mat')

B2010_with_stats

figure;

fit_with_same_RTs = false;
fit_as_experiment_1 = false;
if fit_as_experiment_1
    OG_RT_slope = 10;
    OG_RT_intercept = 205;
    PM_RT_slope = OG_RT_slope;
    PM_RT_intercept = OG_RT_intercept;
elseif fit_with_same_RTs
    OG_RT_slope = RT_slope;
    OG_RT_intercept = RT_intercept;
    PM_RT_slope = RT_slope;
    PM_RT_intercept = RT_intercept;
end

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
plot_all_conditions_brewer(empirical_stats(:, [1:3 4 5]), 650, 850, 1, 0, false);

subplot(4, 2, 4);
plot_all_conditions_brewer(simulation_stats(:, [1:3 4 5]), 650 + fit_as_experiment_1 * 200, 850 + fit_as_experiment_1 * 700, OG_RT_slope, OG_RT_intercept, false);


subplot(4, 2, 5);
ylabel('1st PM RT (msec)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 13 14]), 600, 2500, 1, 0, false);

subplot(4, 2, 6);
plot_all_conditions_brewer(simulation_stats(:, [1:3 8 9]), 600, 2500, PM_RT_slope, PM_RT_intercept, false);

%
% For sanity check
%

subplot(4, 2, 7);
ylabel('OG Accuracy (%)');
plot_all_conditions_brewer(empirical_stats(:, [1:3 6 7]), 00, 100, 1, 0, false);

subplot(4, 2, 8);
plot_all_conditions_brewer(simulation_stats(:, [1:3 6 7]), 00, 100, 1, 0, false);


