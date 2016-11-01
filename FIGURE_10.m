% load('exp1-data.mat')

EM2005_with_stats_exp1

figure;

scatter(simulation_cycles, empirical_RTs, 'fill');
clear xlabel ylabel;
xlabel('Simulation RTs (cycles)');
ylabel('Empirical RTs (msec)');
lsline
text(86, 1550, OG_RT_label_cycles_to_msec, 'fontsize', 12);
title(sprintf('R^2 = %.4gf', rsq));

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);
