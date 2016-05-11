


Block1_RT = data(:, 6);
Block1_RT = Block1_RT(data(:,1) == 0 & data(:, 9) == 1 & data(:, 2) == 1);
Block1_RT = Block1_RT(~isnan(Block1_RT));

mean(Block1_RT) * RT_slope + RT_intercept
std(Block1_RT) * RT_slope / sqrt(length(Block1_RT))

Block4_RT = data(:, 6);
Block4_RT = Block4_RT(data(:,1) == 0 & data(:, 9) == 6 & data(:, 2) == 1);
Block4_RT = Block4_RT(~isnan(Block4_RT));

mean(Block4_RT) * RT_slope + RT_intercept
std(Block4_RT) * RT_slope / sqrt(length(Block4_RT))


all = [Block1_RT; Block4_RT];
groups = [repmat({'1'}, length(Block1_RT), 1); repmat({'6'}, length(Block4_RT), 1)];

[p, table] = anova1(all, groups, 'off')







%{
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
        ylim([5000 6000]);
    end

%}