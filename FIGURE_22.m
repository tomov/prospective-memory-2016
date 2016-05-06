%load('exp3-data.mat')

EM2005_with_stats_exp3




Ms = zeros(1, 2);
SEMs = zeros(1, 2);
targets_idx = 0;
for TARGETS = [1,6]
    targets_idx = targets_idx + 1;
    samples = subjects(subjects(:, 1) == 0 & subjects(:, 9) == TARGETS, 7);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples) / 2); % TODO notice / 2 => b/c we technically are merging subjects in high/low emphasis conditions
    Ms(targets_idx) = M;
    SEMs(targets_idx) = SEM;
end

figure;

subplot(3, 2, 1);
barweb([80 72], [28 25]/sqrt(subjects_per_condition), 1, {}, ...
    'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
ylim([40 100]);

subplot(3, 2, 2);
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'PM Condition');
legend({'1 Target', '6 Targets'});
ylim([40 100]);








Ms = zeros(1, 2);
SEMs = zeros(1, 2);
targets_idx = 0;
for TARGETS = [1,6]
    targets_idx = targets_idx + 1;
    samples = subjects(subjects(:, 1) == 0 & subjects(:, 9) == TARGETS, 5);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples) / 2); % TODO notice / 2 => b/c we technically are merging subjects in high/low emphasis conditions
    Ms(targets_idx) = M;
    SEMs(targets_idx) = SEM;
end

subplot(3, 2, 3);
barweb([69 70], [10 11]/sqrt(subjects_per_condition), 1, {}, ...
    '', 'PM Condition', 'OG Accuracy (%)');
ylim([40 100]);

subplot(3, 2, 4);
barweb(Ms, SEMs, 1, {}, ...
    '', 'PM Condition');
legend({'1 Target', '6 Targets'});
ylim([40 100]);









titles = {'', ''};
sources = {empirical_stats, simulation_stats};
for s_id = 1:2
    subplot(3, 2, s_id + 4);

    stats = sources{s_id};
    Ms = zeros(2);
    SEMs = zeros(2);
    targets_idx = 0;
    for TARGETS = [1,6]
        targets_idx = targets_idx + 1;
        for OG_ONLY = 1:-1:0
            M = stats(stats(:,1) == OG_ONLY & ...
                stats(:, 12) == TARGETS, 4);
            SEM = stats(stats(:,1) == OG_ONLY & ...
                stats(:, 12) == TARGETS, 5);
            if s_id == 2
                M = M * RT_slope + RT_intercept;
                SEM = SEM * RT_slope;
            end
            Ms(2 - OG_ONLY, targets_idx) = M;
            SEMs(2 - OG_ONLY, targets_idx) = SEM;

            if s_id == 2
                fprintf('targets = %d, og only = %d: $%.2f \\pm %.2f$\n', TARGETS, OG_ONLY,  M, SEM);
            end
        end
    end

    barweb(Ms, SEMs, 1, {'1 Target', '6 Targets'}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 1
        ylabel('Ongoing RT (msec)');
    else
        legend({'No-PM', 'PM'});
    end
    ylim([4500 5500]);
end
