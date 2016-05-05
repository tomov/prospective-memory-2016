load('exp4-data.mat')

EM2005_with_stats_exp4


figure;


titles = {'', 'Simulation Data'};
    sources = {empirical_stats, simulation_stats};
    for s_id = 2:2

        figure;

        stats = sources{s_id};
        Ms = zeros(1,2);
        SEMs = zeros(1,2);
            t_id = 0;
            for EMPHASIS = 0:1
                t_id = t_id + 1;
                M = stats(stats(:,1) == 0 & ...
                    stats(:, 3) == EMPHASIS, 8);
                SEM = stats(stats(:,1) == 0 & ...
                    stats(:, 3) == EMPHASIS, 9);
                M = M * RT_slope + RT_intercept;
                SEM = SEM * RT_slope * 4;
                Ms(1, t_id) = M;
                SEMs(1, t_id) = SEM;
                fprintf('emphasis = %d, %.3f +- %.3f\n', EMPHASIS, M, SEM);
            end

        barweb(Ms, SEMs, 1, {}, ...
            titles{s_id}, 'Group');
        if s_id == 2
            h = legend({'Low Cost', 'High Cost'});
            set(h, 'FontSize', 15);
            ylabel('PM hit RT (ms)');
        end
        ylim([6000 11000]);
    end
    h = get(gca, 'xlabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'ylabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'title');
    set(h, 'FontSize', 15);    
