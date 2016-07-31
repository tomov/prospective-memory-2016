%load('exp5-data.mat')

EM2005_with_stats_exp5

%{
figure;
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'Third Task Trial Type');
h = legend({'Target', 'Nontarget'});
set(h, 'FontSize', 15);
ylabel('Third Task RT (ms)');
ylim([1100, 1200]);
 %}       

figure;

subplot(1, 2, 1);

empirical_Ms = empirical_stats(2:3, 4)';
empirical_SEMs = empirical_stats(2:3, 5)';
barweb(empirical_Ms, empirical_SEMs, 1, {}, ...
    'Empirical Data', 'Third Task Trial Type', 'Third Task RT (msec)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([500 700]);

subplot(1, 2, 2);

simulation_Ms = simulation_stats(2:3, 4)' * RT_slope + RT_intercept;
simulation_SEMs = simulation_stats(2:3, 5)' * RT_slope;
barweb(simulation_Ms, simulation_SEMs, 1, {}, ...
    'Simulation Data', 'Third Task Trial Type', 'Third Task RT (msec)');
h = legend({'Non-target', 'Target'});
set(h, 'FontSize', 10);
ylim([1050 1160]);


