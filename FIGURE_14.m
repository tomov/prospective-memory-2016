% load('exp1-data.mat')

EM2005_with_stats_exp1



figure;



titles = {'', ''};
OG_ONLY = 0; % only consider OG with the PM task

Ms = [];
SEMs = [];
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        OG_RTs = data(data(:, 1) == 0 & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS, 4) * RT_slope + RT_intercept;
        OG_RT_M = nanmean(OG_RTs);
        OG_RT_SEM = nanstd(OG_RTs) / sqrt(subjects_per_condition); % SD -> SEM
        PM_miss_RTs = data(data(:, 1) == 0 & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS, 11) * RT_slope + RT_intercept;
        PM_miss_RT_M = nanmean(PM_miss_RTs);
        PM_miss_RT_SEM = nanmean(PM_miss_RTs) / sqrt(subjects_per_condition); % SD -> SEM

        Ms = [Ms; OG_RT_M PM_miss_RT_M];
        SEMs = [SEMs; OG_RT_SEM PM_miss_RT_SEM];
        fprintf('focal = %d, emphasis = %d, OG RT %.3f +- %.3f, PM miss RT %.3f +- %.3f\n', FOCAL, EMPHASIS, OG_RT_M, OG_RT_SEM, PM_miss_RT_M, PM_miss_RT_SEM);
    end
end

barweb(Ms, SEMs, 1, {'Foc,Lo', 'Foc,Hi', 'Nonfoc,Lo', 'Nonfoc,Hi'}, ...
    'Intention Superiority Effect', 'PM Condition');
h = legend({'OG', 'PM miss'});
set(h, 'FontSize', 10);
ylabel('RT (msec)');
ylim([500 2700]);

