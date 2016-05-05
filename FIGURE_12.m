%load('exp1-data.mat')

EM2005_with_stats_exp1







Ms = zeros(1, 2);
SEMs = zeros(1, 2);
for FOCAL = 1:-1:0
    samples = subjects(subjects(:, 1) == 0 & subjects(:, 2) == FOCAL, 7);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    Ms(2 - FOCAL) = M;
    SEMs(2 - FOCAL) = SEM;
    fprintf('focal = %d, %.3f +- %.3f\n', FOCAL, M, SEM);
end

figure;

subplot(3, 2, 1);
barweb([90 67], [16 33]/sqrt(subjects_per_condition), 1, {}, ...
    'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
h = legend({'Foc', 'Nonfoc'});
set(h, 'FontSize', 10);
ylim([30 100]);

subplot(3, 2, 2);
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'PM Condition');
ylim([30 100]);








Ms = zeros(1, 2);
SEMs = zeros(1, 2);
for EMPHASIS = 0:1
    samples = subjects(subjects(:, 1) == 0 & subjects(:, 3) == EMPHASIS, 7);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    Ms(EMPHASIS + 1) = M;
    SEMs(EMPHASIS + 1) = SEM;
    fprintf('emphasis = %d, %.3f +- %.3f\n', EMPHASIS, M, SEM);
end

subplot(3, 2, 3);
barweb([70 87], [32 22]/sqrt(subjects_per_condition), 1, {}, ...
    '', 'PM Condition', 'PM Hit rate (%)');
h =legend({'Low', 'High'});
set(h, 'FontSize', 10);
ylim([30 100]);

subplot(3, 2, 4);
barweb(Ms, SEMs, 1, {}, ...
    '', 'PM Condition');
ylim([30 100]);








titles = {'', ''};
sources = {empirical_stats, simulation_stats};
for s_id = 1:2
    subplot(3, 2, s_id + 4);

    stats = sources{s_id};
    Ms = zeros(2);
    SEMs = zeros(2);
    for FOCAL = 1:-1:0
        for EMPHASIS = 0:1
            M = stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 10);
            SEM = stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 11);
            Ms(2 - FOCAL, EMPHASIS + 1) = M;
            SEMs(2 - FOCAL, EMPHASIS + 1) = SEM;
            fprintf('focal = %d, emphasis = %d, %.3f +- %.3f\n', FOCAL, EMPHASIS, M, SEM);
        end
    end

    barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 1
        h = legend({'Low', 'High'});
        set(h, 'FontSize', 10);
        ylabel('PM Hit rate (%)');
    end
    ylim([30 100]);
end


