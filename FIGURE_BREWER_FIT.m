B2010_with_stats

figure;

subplot(2, 1, 1);

scatter(simulation_cycles, empirical_RTs, 'fill');
clear xlabel ylabel;
xlabel('Simulation OG RTs (cycles)');
ylabel('Empirical OG RTs (msec)');
lsline
text(86, 1550, RT_label_cycles_to_msec, 'fontsize', 12);
title(sprintf('R^2 = %.4f', rsq));

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);

% OG RT's

subplot(2, 2, 3);

scatter(simulation_OG_cycles, empirical_OG_RTs, 'fill');
clear xlabel ylabel;
xlabel('Simulation OG RTs (cycles)');
ylabel('Empirical OG RTs (msec)');
lsline
text(86, 1550, OG_RT_label_cycles_to_msec, 'fontsize', 12);
title(sprintf('R^2 = %.4f', OG_rsq));

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);

% PM RT's

subplot(2, 2, 4);

scatter(simulation_PM_cycles, empirical_PM_RTs, 'fill');
clear xlabel ylabel;
xlabel('Simulation PM RTs (cycles)');
ylabel('Empirical PM RTs (msec)');
lsline
text(86, 1550, PM_RT_label_cycles_to_msec, 'fontsize', 12);
title(sprintf('R^2 = %.4f', PM_rsq));

h = get(gca, 'xlabel');
set(h, 'FontSize', 15);
h = get(gca, 'ylabel');
set(h, 'FontSize', 15);
h = get(gca, 'title');
set(h, 'FontSize', 15);
