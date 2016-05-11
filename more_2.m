
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
    

%{
        stats = simulation_stats;
        Ms = zeros(2);
        SEMs = zeros(2);
        stats = simulation_stats;
        for FOCAL = 1:-1:0
            for BLOCK = 1:4
                M = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 10);
                SEM = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 11);
                M_ogonly = stats(stats(:,1) == 1 & stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 4);
                M_ogonly = M_ogonly * RT_slope + RT_intercept; 
                fprintf('$%.1f\\%% \\pm %.1f\\%%$, ', M, SEM);
            end
            fprintf('\n');
        end
%}
        

%{

Block1_RT = data(:, 7);
Block1_RT = Block1_RT(data(:,1) == 0 & data(:, 10) == 1 & data(:, 2) == 1);

mean(Block1_RT)
std(Block1_RT) / sqrt(length(Block1_RT))

Block4_RT = data(:, 7);
Block4_RT = Block4_RT(data(:,1) == 0 & data(:, 10) == 4 & data(:, 2) == 1);

mean(Block4_RT)
std(Block4_RT) / sqrt(length(Block4_RT))


all = [Block1_RT; Block4_RT];
groups = [repmat({'1'}, length(Block1_RT), 1); repmat({'4'}, length(Block4_RT), 1)];

[p, table] = anova1(all, groups, 'off')
%}



%{
Block1_RT = data(:, 4);
Block1_RT = Block1_RT(data(:,1) == 0 & data(:, 10) == 1 & data(:, 2) == 0);

mean(Block1_RT) * RT_slope + RT_intercept
std(Block1_RT) * RT_slope / sqrt(length(Block1_RT))

Block4_RT = data(:, 4);
Block4_RT = Block4_RT(data(:,1) == 0 & data(:, 10) == 4 & data(:, 2) == 0);

mean(Block4_RT) * RT_slope + RT_intercept
std(Block4_RT) * RT_slope / sqrt(length(Block4_RT))


all = [Block1_RT; Block4_RT];
groups = [repmat({'1'}, length(Block1_RT), 1); repmat({'4'}, length(Block4_RT), 1)];

[p, table] = anova1(all, groups, 'off')
%}


%{
        stats = simulation_stats;
        Ms = zeros(2);
        SEMs = zeros(2);
        stats = simulation_stats;
        for FOCAL = 1:-1:0
            for BLOCK = 1:4
                M = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 4);
                SEM = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 5);
                M = M * RT_slope + RT_intercept;
                SEM = SEM * RT_slope;
                M_ogonly = stats(stats(:,1) == 1 & stats(:, 2) == FOCAL & stats(:, 12) == BLOCK, 4);
                M_ogonly = M_ogonly * RT_slope + RT_intercept; 
                fprintf('focal = %d, block = %d, %.3f +- %.3f, cost = %.3f\n', FOCAL, BLOCK, M, SEM, M - M_ogonly);
            end
            
            
            M = mean(stats(stats(:,1) == 0 & ...
                stats(:, 2) == FOCAL, 4));
            M = M * RT_slope + RT_intercept;
            M_ogonly = stats(stats(:,1) == 1 & stats(:, 2) == FOCAL, 4);
            M_ogonly = mean(M_ogonly);
            M_ogonly = M_ogonly * RT_slope + RT_intercept; 
            fprintf('focal = %d, %.3f, cost = %.3f\n', FOCAL, M, M - M_ogonly);
        end


%}

%{
for BLOCK = 1:4
    for FOCAL = 1:-1:0
        for OG_ONLY = 1:-1:0
            for col = 4
                samples = blocks(blocks(:, 1) == OG_ONLY & blocks(:, 2) == FOCAL & blocks(:, 10) == BLOCK, col);
                M = mean(samples);
                SD = std(samples);
                assert(length(samples) == subjects_per_condition);
                SEM = SD / sqrt(subjects_per_condition);
                M = M * RT_slope + RT_intercept;
                SEM = SEM * RT_slope;
                fprintf('& %.2f $\\pm$ %.2f ', M, SEM);
            end
        end
    end
    fprintf('\n');
end

%}