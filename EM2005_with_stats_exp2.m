%function stat = EM2005_analyze_stats( subjects, subjects_extra)
% get stats from subjects and analyze so you can fit with fitparam()
%{
 subjects[subj_id, :] ==
 samples variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
    9 - subject id
    10 - block id
 statistics order = 
    4 - OG_RT, 
    5 - OG_Hit, 
    6 - PM_RT, 
    7 - PM_Hit,
    8 - PM_miss_OG_hit
    11 - first_PM_RT
    12 - PM_miss_OG_RT 
 (see EM2005 exp 2)
%}

DO_PLOT = false;
blocks = data;


% -------------- define the empirical stats (Table 2 from E&M 2005)

%{
 stats variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS  (N/A in this case, just to make compatible with exp 1
    12 - block id
 statistics order = 
    4 - OG_RT_M,
    5 - OG_RT_SEM,
    6 - OG_Hit_M,
    7 - OG_Hit_SEM,
    8 - PM_RT_M,
    9 - PM_RT_SEM
    10 - PM_Hit_M,
    11 - PM_Hit_SEM
%}
SD_cols = [5 7 9 11];

subjects_per_condition = 24;

empirical_stats = [
    % "trial" (= block) 1
    1 1 0, 1108.76, 203.35,   96, 3, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    low emph
    0 1 0, 1160.54, 195.68,   96, 3, NaN, NaN, 92, 28,   1;  % PM,    focal,    low emph
    1 0 0,  969.02, 116.40,   96, 3, NaN, NaN, NaN, NaN, 1;  % no-PM, nonfocal, low emph
    0 0 0, 1240.36, 215.61,   96, 3, NaN, NaN, 71, 46,   1;  % PM,    nonfocal, low emph

    % "trial" (= block) 2
    1 1 0, 1122.50, 171.40,   96, 3, NaN, NaN, NaN, NaN, 2;  % no-PM, focal,    low emph
    0 1 0, 1096.72, 136.01,   96, 3, NaN, NaN, 96, 20,   2;  % PM,    focal,    low emph
    1 0 0,  971.36, 124.89,   96, 3, NaN, NaN, NaN, NaN, 2;  % no-PM, nonfocal, low emph
    0 0 0, 1190.24, 205.00,   96, 3, NaN, NaN, 71,  46,  2;  % PM,    nonfocal, low emph

    % "trial" (= block) 3
    1 1 0, 1096.18, 155.22,   96, 3, NaN, NaN, NaN, NaN, 3;  % no-PM, focal,    low emph
    0 1 0, 1175.57, 248.35,   96, 3, NaN, NaN, 96, 20,   3;  % PM,    focal,    low emph
    1 0 0,  957.19, 119.83,   96, 3, NaN, NaN, NaN, NaN, 3;  % no-PM, nonfocal, low emph
    0 0 0, 1122.55, 172.26,   96, 3, NaN, NaN, 63, 50,   3;  % PM,    nonfocal, low emph

    % "trial" (= block) 4
    1 1 0, 1108.99, 203.15,   96, 3, NaN, NaN, NaN, NaN, 4;  % no-PM, focal,    low emph
    0 1 0, 1120.33, 170.92,   96, 3, NaN, NaN, 88, 34,   4;  % PM,    focal,    low emph
    1 0 0,  952.68, 112.84,   96, 3, NaN, NaN, NaN, NaN, 4;  % no-PM, nonfocal, low emph
    0 0 0, 1090.58, 185.73,   96, 3, NaN, NaN, 42, 50,   4;  % PM,    nonfocal, low emph
];



% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);

% ------------- calculate simulation stats (Table 2 from E&M 2005)

simulation_stats = [];
EMPHASIS = 0;
% order here matters -- must be same as empirical_data above for line
% regression
for BLOCK = 1:4
    for FOCAL = 1:-1:0
        for OG_ONLY = 1:-1:0
            stat = [OG_ONLY, FOCAL, EMPHASIS];
            for col = 4:7
                samples = blocks(blocks(:, 1) == OG_ONLY & blocks(:, 2) == FOCAL & blocks(:, 10) == BLOCK, col);
                %assert(sum(isnan(samples)) == 0);
                %samples = samples(~isnan(samples));
                M = mean(samples);
                SD = std(samples);
               % assert(length(samples) == subjects_per_condition);
                stat = [stat, M, SD];
            end
            stat = [stat, BLOCK];
            simulation_stats = [simulation_stats; stat];
        end
    end
end

% convert SD's to SEM's in simulation data
simulation_stats(:, SD_cols) = simulation_stats(:, SD_cols) / sqrt(subjects_per_condition);



% -------------- run linear regression to find slope and intercept for RT's

simulation_cycles = simulation_stats(:, 4);
empirical_RTs = empirical_stats(:, 4);

% TODO HACK FIXME -- to "fix" E&M's wacky results, simply make the
% Focal, PM and the Focal, No-PM condition same as the Nonfocal, No-PM
% condition
% a little fraud b/c the Focal, PM is slightly above the No-PM conditions generally 
RTs_to_ignore_when_fitting = empirical_stats(:, 1) == 1 & empirical_stats(:, 2) == 1;
%RTs_to_ignore_when_fitting = RTs_to_ignore_when_fitting | (empirical_stats(:, 1) == 0 & empirical_stats(:, 2) == 1);
empirical_RTs = empirical_RTs(~RTs_to_ignore_when_fitting);
simulation_cycles = simulation_cycles(~RTs_to_ignore_when_fitting);

p = polyfit(simulation_cycles, empirical_RTs, 1);
RT_slope = p(1);
RT_intercept = p(2);
yfit = polyval(p, simulation_cycles);

yresid = empirical_RTs - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(empirical_RTs)-1) * var(empirical_RTs);
rsq = 1 - SSresid/SStotal;

%TODO FIXME use slope and intercept from experiment 1 UGH or not...
%RT_slope = 8;
%RT_intercept = 263;

OG_RT_label_cycles_to_msec = sprintf('OG RT (msec = cycles * %.1f + %.1f)', RT_slope, RT_intercept);

if DO_PLOT
    figure;
    scatter(simulation_cycles, empirical_RTs, 'fill');
    clear xlabel ylabel;
    xlabel('Simulation RTs (cycles)');
    ylabel('Empirical RTs (msec)');
    lsline
    text(min(xlim) + 1, 1220, OG_RT_label_cycles_to_msec, 'fontsize', 14);
    title(sprintf('R^2 = %.4f', rsq));
    h = get(gca, 'xlabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'ylabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'title');
    set(h, 'FontSize', 15);

end    





% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% ----------------------------- COMPARE EMPIRICAL AND SIMULATION DATA ----
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------





% ---------------------------------------------------
% ----------------------------- PM PERFORMANCE ------
% ---------------------------------------------------




PM_hit = blocks(:, 7);

% ----------------------- PM hit rate in focal vs. nonfocal ----

PM_hit_focal = PM_hit(blocks(:, 2) == 1);
PM_hit_nonfocal = PM_hit(blocks(:, 2) == 0);
[p, table] = anova1([PM_hit_focal PM_hit_nonfocal], {'Focal', 'Nonfocal'}, 'off');

fprintf('\n\n----- PM Performance: Focal vs. Nonfocal ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 18.38\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,5}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for FOCAL = 1:-1:0
        samples = blocks(blocks(:, 1) == 0 & blocks(:, 2) == FOCAL, 7);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        Ms(2 - FOCAL) = M;
        SEMs(2 - FOCAL) = SEM;
    end
    
    figure;
    
    subplot(1, 2, 1);
    barweb([93 61], [16 32]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
    legend({'Focal', 'Nonfocal'});
    ylim([30 100]);

    subplot(1, 2, 2);
    barweb(Ms, SEMs, 1, {}, ...
        'Simulation Data', 'PM Condition');
    ylim([30 100]);
end



% --------------------- PM hit rate in diff. blocks ----


[p, table] = anovan(PM_hit, {blocks(:, 10)}, 'model','full', 'display', 'off');

fprintf('\n\n----- PM Performance: Blocks ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 2.99\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


% --------------------- PM hit rate in diff blocks, split by focal vs. nonfocal

empirical_Fs = [
    8.76;  % nonfocal
    0.5;   % focal (F < 1)
    ];
focal_titles = {'Nonfocal', 'Focal'};
for FOCAL = 1:-1:0
    samples = blocks(blocks(:, 1) == 0 & blocks(:, 2) == FOCAL, :);
    [p, table] = anovan(samples(:, 7), {samples(:, 10)}, 'model','full', 'display', 'off');
    fprintf('\n\n----- PM hit rate: %s ------\n', focal_titles{FOCAL+1});
    fprintf('\n  Empirical Data -------\n');
    fprintf('                 F = %.2f\n', empirical_Fs(FOCAL+1, EMPHASIS+1));
    fprintf('\n  Simulation Data -------\n');
    fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));
end






% ---------------------------------------------------
% ----------------------------- OG PERFORMANCE ------
% ---------------------------------------------------




% ------------------ OG accuracy

OG_hit = blocks(:, 5);

[p, table] = anovan(OG_hit, {blocks(:, 1) blocks(:, 2) blocks(:, 10)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG accuracy: 2x2x4 ANOVA ------\n');
table(1:8,6)
p
fprintf('none of these should be significant (i.e. you SHOULD have F < 1 everywhere)\n');




% ----------------- OG RT: PM vs. No PM


OG_RTs = blocks(:, 4);

[p, table] = anovan(OG_RTs, {blocks(:, 1)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: PM vs. No PM ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 33.15\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


% ----------------- OG RT: 2x2x2 ANOVA


[p, table] = anovan(OG_RTs, {blocks(:, 1) blocks(:, 2) blocks(:, 10)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: 2x2x4 ANOVA ------\n');
table(1:8,6)
p
fprintf(' F(3, 138) = the 3-way interaction => (empirical data) F = 3.71');

% ----------------- OG RT cost: focal vs. nonfocal


M_OG_only = mean(OG_RTs(blocks(:, 1) == 1));

empirical_Fs = [
    50.38; % nonfocal
    1.08   % focal
    ];
empirical_Fs_blocks = [
    7.08; % nonfocal
    0.5   % focal ( F < 1)
    ];

focal_titles = {'Nonfocal', 'Focal'};
for FOCAL = 1:-1:0
    samples = blocks(blocks(:, 2) == FOCAL, :);
    [p, table] = anovan(samples(:,4), {samples(:, 1)}, 'model','full', 'display', 'off');
    fprintf('\n\n----- OG RT Cost: %s ------\n', focal_titles{FOCAL+1});
    fprintf('\n  Empirical Data -------\n');
    fprintf('                 F = %.2f\n', empirical_Fs(FOCAL+1));
    fprintf('\n  Simulation Data -------\n');
    fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

    % remove the OG_ONLY runs and only compare first and last block
    samples = blocks(blocks(:, 1) == 0 & blocks(:, 2) == FOCAL & (blocks(:, 10) == 1 | blocks(:, 10) == 4), :);
    [p, table] = anovan(samples(:,4), {samples(:, 10)}, 'model','full', 'display', 'off');
    fprintf('\n\n----- OG RT Cost steady decrease: %s ------\n', focal_titles{FOCAL+1});
    fprintf('\n  Empirical Data -------\n');
    fprintf('                 F = %.2f\n', empirical_Fs_blocks(FOCAL+1));
    fprintf('\n  Simulation Data -------\n');
    fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));
end


%{

[p, table] = anovan(OG_RTs, {subjects(:, 2)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: Focal vs. Nonfocal ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 22.87\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for FOCAL = 1:-1:0
        samples = subjects(subjects(:, 2) == FOCAL, 4);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        Ms(2 - FOCAL) = M * RT_slope + RT_intercept;
        SEMs(2 - FOCAL) = SEM * RT_slope;
    end
    
    figure;
    
    subplot(3, 2, 1);
    barweb([1145.63 1335.73], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'Ongoing RT (msec)');
    legend({'Focal', 'Nonfocal'});
    ylim([1000 1400]);

    subplot(3, 2, 2);
    barweb(Ms, SEMs, 1, {}, ...
        'Simulation Data', 'PM Condition');
    ylim([1000 1400]);
end


% ----------------- OG RT: high vs. low emphasis

[p, table] = anovan(OG_RTs, {subjects(:, 3)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: High vs. Low emphasis ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 6.47\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for EMPHASIS = 0:1
        samples = subjects(subjects(:, 3) == EMPHASIS, 4);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        Ms(EMPHASIS + 1) = M * RT_slope + RT_intercept;
        SEMs(EMPHASIS + 1) = SEM * RT_slope;
    end
    
    subplot(3, 2, 3);
    barweb([1190.11 1291.26], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'Ongoing RT (msec)');
    legend({'Low Emphasis', 'High Emphasis'});
    ylim([1000 1400]);

    subplot(3, 2, 4);
    barweb(Ms, SEMs, 1, {}, ...
        'Simulation Data', 'PM Condition');
    ylim([1000 1400]);
end







% ----------------- OG RT: 2x2x2 ANOVA



% ----------------- OG RT: cost qualified

% cost is implied in figures above, NEXT...



% ----------------- OG RT: cost, interaction focality & emphasis


M_OG_only = mean(OG_RTs(subjects(:, 1) == 1));

empirical_Fs = [
    61.52, 127.96; % nonfocal: low, high
    1.73, 6.15     % focal: low, high
    ];
focal_titles = {'Nonfocal', 'Focal'};
emphasis_titles = {'Low Emphasis', 'High Emphasis'};
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        samples = subjects(subjects(:, 2) == FOCAL & subjects(:, 3) == EMPHASIS, :);
        [p, table] = anovan(samples(:,4), {samples(:, 1)}, 'model','full', 'display', 'off');
        fprintf('\n\n----- OG RT Cost: %s, %s ------\n', focal_titles{FOCAL+1}, emphasis_titles{EMPHASIS+1});
        fprintf('\n  Empirical Data -------\n');
        fprintf('                 F = %.2f\n', empirical_Fs(FOCAL+1, EMPHASIS+1));
        fprintf('\n  Simulation Data -------\n');
        fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));
    end
end


if DO_PLOT
    titles = {'Empirical Data', 'Simulation Data'};
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
            end
        end

        barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
            titles{s_id}, 'PM Condition');
        if s_id == 1
            legend({'Low Emphasis', 'High Emphasis'});
            ylabel('Ongoing RT (msec)');
        end
        ylim([1000 1700]);
    end
end




%}


% ---------------------------------------------------
% ----------------------------- TABLE 1 -------------
% ---------------------------------------------------


if DO_PLOT
    figure;

    subplot(3, 2, 1);
    title('Empirical Data');
    ylabel('OG RT (msec)');
    plot_all_conditions_exp2(empirical_stats(:, [1:3 4 5 12]), 900, 1300, 1, 0, false, [0 1]);

    subplot(3, 2, 2);
    title('Simulation Data');
    plot_all_conditions_exp2(simulation_stats(:, [1:3 4 5 12]), 900, 1300, RT_slope, RT_intercept, true, [0 1]);

    subplot(3, 2, 3);
    ylabel('OG Accuracy (%)');
    plot_all_conditions_exp2(empirical_stats(:, [1:3 6 7 12]), 20, 100, 1, 0, false, [0 1]);

    subplot(3, 2, 4);
    plot_all_conditions_exp2(simulation_stats(:, [1:3 6 7 12]), 20, 100, 1, 0, false, [0 1]);

    subplot(3, 2, 5);
    ylabel('PM Hit Rate (%)');
    plot_all_conditions_exp2(empirical_stats(:, [1:3 10 11 12]), 20, 100, 1, 0, false, [0]);

    subplot(3, 2, 6);
    plot_all_conditions_exp2(simulation_stats(:, [1:3 10 11 12]), 20, 100, 1, 0, false, [0]);
end
