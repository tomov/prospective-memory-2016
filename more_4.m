

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






GroupA = data(:, 6);
GroupA = GroupA(data(:,1) == 0 & data(:, 3) == 0);

mean(GroupA) * RT_slope + RT_intercept
std(GroupA) * RT_slope / sqrt(length(GroupA))

GroupB = data(:, 6);
GroupB = GroupB(data(:,1) == 0 & data(:, 3) == 1);

mean(GroupB) * RT_slope + RT_intercept
std(GroupB) * RT_slope / sqrt(length(GroupB))


all = [GroupA; GroupB];
groups = [repmat({'1'}, length(GroupA), 1); repmat({'4'}, length(GroupB), 1)];

[p, table] = anova1(all, groups, 'off')




%{
Block1_RT = data(:, 4);
Block1_RT = Block1_RT(data(:,1) == 1 & data(:, 3) == 0);

mean(Block1_RT) * RT_slope + RT_intercept
std(Block1_RT) * RT_slope / sqrt(length(Block1_RT))

Block4_RT = data(:, 4);
Block4_RT = Block4_RT(data(:,1) == 1 & data(:, 3) == 1);

mean(Block4_RT) * RT_slope + RT_intercept
std(Block4_RT) * RT_slope / sqrt(length(Block4_RT))


all = [Block1_RT; Block4_RT];
groups = [repmat({'1'}, length(Block1_RT), 1); repmat({'4'}, length(Block4_RT), 1)];

[p, table] = anova1(all, groups, 'off')

%}