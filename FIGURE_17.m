load('exp2-data.mat')

EM2005_with_stats_exp2




Ms = zeros(1, 2);
SEMs = zeros(1, 2);
for FOCAL = 1:-1:0
    samples = blocks(blocks(:, 1) == 0 & blocks(:, 2) == FOCAL, 7);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    Ms(2 - FOCAL) = M;
    SEMs(2 - FOCAL) = SEM;
end

figure;

subplot(1, 2, 1);
barweb([93 61], [16 32]/sqrt(subjects_per_condition), 1, {}, ...
    'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
legend({'Focal', 'Nonfocal'});
ylim([30 100]);

subplot(1, 2, 2);
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'PM Condition');
ylim([30 100]);