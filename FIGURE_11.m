% load('exp1-data.mat')

EM2005_with_stats_exp1







Ms = zeros(1, 2);
SEMs = zeros(1, 2);
for FOCAL = 1:-1:0
    samples = subjects(subjects(:, 2) == FOCAL, 4);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    M = M * RT_slope + RT_intercept;
    SEM = SEM * RT_slope;
    Ms(2 - FOCAL) = M;
    SEMs(2 - FOCAL) = SEM;
    fprintf('focal = %d, %.3f +- %.3f\n', FOCAL, M, SEM);
end

figure;

subplot(3, 2, 1);
barweb([1145.63 1335.73], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
    'Empirical Data', 'PM Condition', 'Ongoing RT (msec)');
h = legend({'Focal', 'Nonfoc'});
set(h, 'FontSize', 10);
ylim([1000 1400]);

subplot(3, 2, 2);
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'PM Condition');
ylim([1000 1400]);








Ms = zeros(1, 2);
SEMs = zeros(1, 2);
for EMPHASIS = 0:1
    samples = subjects(subjects(:, 3) == EMPHASIS, 4);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    M = M * RT_slope + RT_intercept;
    SEM = SEM * RT_slope;
    Ms(EMPHASIS + 1) = M;
    SEMs(EMPHASIS + 1) = SEM;
    fprintf('emphasis = %d, %.3f +- %.3f\n', EMPHASIS, M, SEM);
end

subplot(3, 2, 3);
barweb([1190.11 1291.26], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
    '', 'PM Condition', 'Ongoing RT (msec)');
h = legend({'Low', 'High'});
set(h, 'FontSize', 10);
ylim([1000 1400]);

subplot(3, 2, 4);
barweb(Ms, SEMs, 1, {}, ...
    '', 'PM Condition');
ylim([1000 1400]);









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
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 4);
            SEM = stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 5);
            if s_id == 2
                M = M * RT_slope + RT_intercept;
                SEM = SEM * RT_slope;
            end
            Ms(2 - FOCAL, EMPHASIS + 1) = M;
            SEMs(2 - FOCAL, EMPHASIS + 1) = SEM;
            M_ogonly = stats(stats(:,1) == 1 & stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 4);
            M_ogonly = M_ogonly * RT_slope + RT_intercept; 
            fprintf('focal = %d, emphasis = %d, %.3f +- %.3f, cost = %.3f\n', FOCAL, EMPHASIS, M, SEM, M - M_ogonly);
        end
    end

    barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 1
        h = legend({'Low', 'High'});
        set(h, 'FontSize', 10);
        ylabel('Ongoing RT (msec)');
    end
    ylim([1000 1700]);
end

