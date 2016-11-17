EM2005_with_stats_exp1

fprintf('\n\n');

OG_only = logical(subjects(:, 1));
focality = logical(subjects(:, 2));
emphasis = logical(subjects(:, 3));

sem = @(x) std(x) / sqrt(length(x));

% OG accuracy PM vs. No PM
OG_acc = subjects(:, 5);
fprintf(['Accuracy on the ongoing task was near ceiling in all simulations (M = %.2f +- %.2f %%), ' ...
    'both in the presence of a PM instruction (M = %.2f +- %.2f %%) as well as in its absence (%.2f +- %.2f %%)\n'], ...
    mean(OG_acc(:)), sem(OG_acc(:)), ...
    mean(OG_acc(~OG_only)), sem(OG_acc(~OG_only)), ...
    mean(OG_acc(OG_only)), sem(OG_acc(OG_only)));

fprintf('\n');

% OG RT PM vs No PM
OG_RTs = subjects(:, 4) * RT_slope + RT_intercept;
[p, table] = anovan(OG_RTs, {OG_only}, 'model','full', 'display', 'off');
fprintf(['Ongoing RTs were slower in the presence of a PM instruction (M = %.2f +- %.2f ms), ', ...
    'compared to performing the OG task alone (M = %.2f +- %.2f ms; F = %f, p = %e)\n'], ...
    mean(OG_RTs(~OG_only)), sem(OG_RTs(~OG_only)), ...
    mean(OG_RTs(OG_only)), sem(OG_RTs(OG_only)), ...
    table{2,6}, p(1));

% OG RT focality
[p, table] = anovan(OG_RTs(~OG_only), {focality(~OG_only)}, 'model','full', 'display', 'off');
fprintf(['In the presence of a PM instruction, ', ...
    'ongoing RTs were significantly faster when the PM target was focal (M = %.2f +- %.2f ms) ', ...
    'compared to when it was nonfocal (M = %.2f +- %.2f ms; F = %f, p = %e).\n'], ... 
    mean(OG_RTs(~OG_only & focality)), sem(OG_RTs(~OG_only & focality)), ...
    mean(OG_RTs(~OG_only & ~focality)), sem(OG_RTs(~OG_only & ~focality)), ...
    table{2,6}, p(1));

% OG RT emphasis
[p, table] = anovan(OG_RTs(~OG_only), {emphasis(~OG_only)}, 'model','full', 'display', 'off');
fprintf(['Also in the presence of a PM instruction, ', ...
    'ongoing RTs were significantly faster when the PM task had a lower priority (or low emphasis; M = %.2f +- %.2f ms) ', ...
    'compared to when it had a higher priority (or high emphasis; M = %.2f +- %.2f ms; F = %f, p = %e).\n'], ... 
    mean(OG_RTs(~OG_only & ~emphasis)), sem(OG_RTs(~OG_only & ~emphasis)), ...
    mean(OG_RTs(~OG_only & emphasis)), sem(OG_RTs(~OG_only & emphasis)), ...
    table{2,6}, p(1));

% OG RT focality X emphasis
[p, table] = anovan(OG_RTs(~OG_only), {focality(~OG_only) emphasis(~OG_only)}, 'model','full', 'display', 'off');
fprintf(['There was a nominal although insignificant interaction between these two effects (F = %f, p = %f) ', ...
    'with focal, low emphasis runs having the fastest OG RTs (or lowest OG cost; M = %.2f +- %.2f ms), ', ...
    'followed by focal, high emphasis runs (M = %.2f +- %.2f ms), ', ...
    'followed by nonfocal, low emphasis runs (M = %.2f +- %.2f ms), ', ...
    'followed by nonfocal, high emphasis runs which had the slowest OG RTs (or highest OG cost; M = %.2f +- %.2f ms)\n'], ...
    table{4,6}, p(3), ...
    mean(OG_RTs(~OG_only & focality & ~emphasis)), sem(OG_RTs(~OG_only & focality & ~emphasis)), ...
    mean(OG_RTs(~OG_only & focality & emphasis)), sem(OG_RTs(~OG_only & focality & emphasis)), ...
    mean(OG_RTs(~OG_only & ~focality & ~emphasis)), sem(OG_RTs(~OG_only & ~focality & ~emphasis)), ...
    mean(OG_RTs(~OG_only & ~focality & emphasis)), sem(OG_RTs(~OG_only & ~focality & emphasis)));

fprintf('\n');

% PM hit rate focality
PM_hit = subjects(:, 7);
[p, table] = anovan(PM_hit(~OG_only), {focality(~OG_only)}, 'model','full', 'display', 'off');
fprintf(['PM hit rates were significantly higher when the PM target was focal (M = %.2f +- %.2f %%) ', ...
    'compared to when it was nonfocal (M = %.2f +- %.2f %%; F = %f, p = %e).\n'], ... 
    mean(PM_hit(~OG_only & focality)), sem(PM_hit(~OG_only & focality)), ...
    mean(PM_hit(~OG_only & ~focality)), sem(PM_hit(~OG_only & ~focality)), ...
    table{2,6}, p(1));

% PM hit rate emphasis
[p, table] = anovan(PM_hit(~OG_only), {emphasis(~OG_only)}, 'model','full', 'display', 'off');
fprintf(['PM hit rates were significantly higher when the PM task had a higher emphasis (M = %.2f +- %.2f %%) ', ...
    'compared to when it had a lower emphasis (M = %.2f +- %.2f %%; F = %f, p = %f).\n'], ... 
    mean(PM_hit(~OG_only & emphasis)), sem(PM_hit(~OG_only & emphasis)), ...
    mean(PM_hit(~OG_only & ~emphasis)), sem(PM_hit(~OG_only & ~emphasis)), ...
    table{2,6}, p(1));

% PM hit rate focality X emphasis
[p, table] = anovan(PM_hit(~OG_only), {focality(~OG_only) emphasis(~OG_only)}, 'model','full', 'display', 'off');
fprintf('There was a significant interaction between these two effects (F = %f, p = %f).\n', ...
    table{4,6}, p(3));
[p, table] = anovan(PM_hit(~OG_only & focality), {emphasis(~OG_only & focality)}, 'model','full', 'display', 'off');
fprintf(['Specifically, for focal targets, the priority of the PM task was irrelevant and PM hit rates were near-ceiling ', ...
    'both when the emphasis was high (M = %.2f +- %.2f %%) ', ...
    'or low (M = %.2f +- %.2f %%; F = %f, p = %f).\n'], ...
    mean(PM_hit(~OG_only & focality & emphasis)), sem(PM_hit(~OG_only & focality & emphasis)), ...
    mean(PM_hit(~OG_only & focality & ~emphasis)), sem(PM_hit(~OG_only & focality & ~emphasis)), ...
    table{2,6}, p(1));
[p, table] = anovan(PM_hit(~OG_only & emphasis), {focality(~OG_only & emphasis)}, 'model','full', 'display', 'off');
fprintf(['For nonfocal targets, however, PM hit rates were significantly higher ', ...
    'when the emphasis was high (M = %.2f +- %.2f %%), ', ...
    'compared to when it was low (M = %.2f +- %.2f %%; F = %f, p = %f).\n'], ...
    mean(PM_hit(~OG_only & ~focality & emphasis)), sem(PM_hit(~OG_only & ~focality & emphasis)), ...
    mean(PM_hit(~OG_only & ~focality & ~emphasis)), sem(PM_hit(~OG_only & ~focality & ~emphasis)), ...
    table{2,6}, p(1));

% Intention superiority effect
PM_miss_OG_RT = subjects(~isnan(subjects(:, 11)), 11) * RT_slope + RT_intercept;
[p, table] = anovan([PM_miss_OG_RT; OG_RTs], {[ones(length(PM_miss_OG_RT), 1); zeros(length(OG_RTs), 1)]}, 'model','full', 'display', 'off');
fprintf(['Our model showed a robust intention superiority effect: ongoing RTs (M = %.2f +- %.2f ms) ', ...
    'were slower than RTs on PM miss trials that resulted in correct OG responses ', ...
    '(i.e. trials on which the target item was presented but the model performed the OG task instead; ', ...
    'M = %.2f +- %.2f ms; F = %f, p = %f).\n'], ...
    mean(OG_RTs), sem(OG_RTs), ...
    mean(PM_miss_OG_RT), sem(PM_miss_OG_RT), ...
    table{2,6}, p(1));