load('exp1-data.mat')

EM2005_with_stats_exp1

titles = {'', 'Simulation Data'};
sources = {empirical_stats, simulation_stats};
for s_id = 2:2

    figure;

    stats = sources{s_id};
    Ms = zeros(2);
    SEMs = zeros(2);
    for FOCAL = 1:-1:0
        for EMPHASIS = 0:1
            M = stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 8);
            SEM = stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 9);
            M = M * RT_slope + RT_intercept;
            SEM = SEM * RT_slope;
            Ms(2 - FOCAL, EMPHASIS + 1) = M;
            SEMs(2 - FOCAL, EMPHASIS + 1) = SEM;
            fprintf('focal = %d, emphasis = %d, %.3f +- %.3f\n', FOCAL, EMPHASIS, M, SEM);
        end
    end

    barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 2
        h = legend({'Low', 'High'});
        set(h, 'FontSize', 10);
        ylabel('PM hit RT (ms)');
    end
    %ylim([30 100]);
end
