load('exp2-data.mat')

EM2005_with_stats_exp2




PM_RTs = data(:, 6);

all = [];
groups = [];

titles = {'', 'Simulation Data'};
sources = {empirical_stats, simulation_stats};
for s_id = 2:2

    figure;

    stats = sources{s_id};
    Ms = zeros(4,2);
    SEMs = zeros(4,2);
    for BLOCK = 1:4
        for FOCAL = 1:-1:0
            samples = PM_RTs(data(:, 2) == FOCAL & data(:, 10) == BLOCK);
            samples = samples(~isnan(samples));

            if FOCAL && (BLOCK == 1 || BLOCK == 4)
                all = [all; samples];
                groups = [groups; repmat(BLOCK, length(samples), 1)];
            end

            M = mean(samples);
            M = M * RT_slope + RT_intercept;

            SEM = std(samples);
            SEM = SEM / sqrt(length(samples)) * RT_slope;

            Ms(BLOCK, 2 - FOCAL) = M;
            SEMs(BLOCK, 2 - FOCAL) = SEM;
            fprintf('focal = %d, block = %d, $%.2f \\pm %.2f$\n', FOCAL, BLOCK, M, SEM);
        end
    end

    barweb(Ms, SEMs, 1, {'Block #1', 'Block #2', 'Block #3', 'Block #4'}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 2
        h = legend({'Focal', 'Nonfocal'});
        set(h, 'FontSize', 15);
        ylabel('PM hit RT (ms)');
    end
    %ylim([30 100]);
end



[p, table] = anova1(all, groups, 'off')