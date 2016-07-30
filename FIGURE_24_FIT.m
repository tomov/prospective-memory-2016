%  load('exp3-data.mat')

EM2005_with_stats_exp4



figure;
scatter(simulation_cycles, empirical_RTs, 'fill');
clear xlabel ylabel;
xlabel('Simulation RTs (cycles)');
ylabel('Empirical RTs (msec)');
text(83, 5150, OG_RT_label_cycles_to_msec, 'fontsize', 12);
lsline
title(sprintf('R^2 = %.4f', rsq));

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);