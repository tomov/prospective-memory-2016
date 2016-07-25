%load('exp3-data.mat')

EM2005_with_stats_exp3




titles = {'', 'Simulation Data'};
sources = {empirical_stats, simulation_stats};
for s_id = 2:2

    figure;

    stats = sources{s_id};
    Ms = zeros(1,2);
    SEMs = zeros(1,2);
        t_id = 0;
        for TARGETS = [1,6]
            t_id = t_id + 1;
            M = stats(stats(:,1) == 0 & ...
                stats(:, 2) == 1 & stats(:, 12) == TARGETS, 8);
            SEM = stats(stats(:,1) == 0 & ...
                stats(:, 2) == 1 & stats(:, 12) == TARGETS, 9);
            M = M * RT_slope + RT_intercept;
            SEM = SEM * RT_slope;
            Ms(1, t_id) = M;
            SEMs(1, t_id) = SEM;
            fprintf('focal = %d, targets = %d, %.3f +- %.3f\n', FOCAL, TARGETS, M, SEM);
        end

    barweb(Ms, SEMs, 1, {}, ...
        titles{s_id}, 'PM Condition');
    if s_id == 2
        h = legend({'1 Target', '6 Targets'});
        set(h, 'FontSize', 10);
        ylabel('PM hit RT (ms)');
    end
    ylim([4000 5000]);
end
