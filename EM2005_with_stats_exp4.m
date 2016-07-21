%function stat = EM2005_analyze_stats( subjects, subjects_extra)
% get stats from subjects and analyze so you can fit with fitparam()
%{
 subjects[subj_id, :] ==
 samples variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
    9 - TARGETS
 statistics order = 
    4 - OG_RT, 
    5 - OG_Hit, 
    6 - PM_RT, 
    7 - PM_Hit,
    8 - PM_miss_OG_hit
    10 - first_PM_RT
    11 - PM_miss_OG_RT
 (see EM2005 exp 2)
%}

DO_PLOT = false;
subjects = data;


% -------------- define the empirical stats (Table 2 from E&M 2005)

%{
 stats variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - STRATEGY (!!!!!!!!!!!!!! not emphasis!!)
    12 - # of targets
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

subjects_per_condition = 104 / 2; % divide by 2! b/c "EMPHASIS" == "high cost / low cost" is not a real condition => samples will be in half

empirical_stats = [
    % strategy 0 = no cost participants (no monitoring)
    1 1 0, 7201.28,  1520.06,       74, 10, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    low emph, 1 targets
    0 1 0, 6706.30,  1337.78,       73, 10, NaN, NaN, 94,   14, 1;  % PM, focal,    low emph, 1 targets
    
    % strategy 1 = cost participants (monitoring)
    1 1 1, 6797.26, 1484.68,        74, 10, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    high emph, 1 targets
    0 1 1, 7533.95, 1769.18,        73, 10, NaN, NaN, 95,   17, 1;  % PM, focal,    high emph, 1 targets
];


% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);


% -------------- change meaning of third column (formerly EMPHASIS) !!!
% now it means whether it's "low cost" (0) or "high cost" (1) group

% first step -- clear the nan's
assert(sum(isnan(subjects(:, 4))) < 6); % don't want too many nan's
subjects = subjects(~isnan(subjects(:, 4)), :);

for OG_ONLY = 1:-1:0
    which = subjects(:, 1) == OG_ONLY & ~isnan(subjects(:, 4));
    m = median(subjects(which, 4));
    subjects(which, 3) = (subjects(which, 4) > m);
end

% ------------- calculate simulation stats

simulation_stats = [];
STRATEGY = 1;
FOCAL = 1;
TARGETS = 1;
% order here matters -- must be same as empirical_data above for line
% regression
for STRATEGY = 0:1
    for OG_ONLY = 1:-1:0
        stat = [OG_ONLY, FOCAL, STRATEGY];
        for col = 4:7
            samples = subjects(subjects(:, 1) == OG_ONLY & subjects(:, 3) == STRATEGY & subjects(:, 9) == TARGETS, col);
            samples = samples(~isnan(samples));
            M = mean(samples);
            SD = std(samples);
            SEM = SD / sqrt(length(samples));
            % this assert blows b/c we get one or two NaN's... ndb
            % instead, assert we're within a reasonable range
            %assert(length(samples) == subjects_per_condition);
            %assert(length(samples) >= subjects_per_condition - 4)
            
            %if col == 4
            %    fprintf('OG RT: emphasis = %d, og only = %d: $%.2f \\pm %.2f$\n', EMPHASIS, OG_ONLY,  M * RT_slope + RT_intercept, SEM * RT_slope);
            %elseif col == 7
            %    fprintf('PM hit rate: emphasis = %d, og only = %d: $%.2f \\pm %.2f$\n', EMPHASIS, OG_ONLY,  M, SEM);
            %end
            stat = [stat, M, SD];
        end
        stat = [stat, TARGETS];
        simulation_stats = [simulation_stats; stat];
    end
end

% convert SD's to SEM's in simulation data
simulation_stats(:, SD_cols) = simulation_stats(:, SD_cols) / sqrt(subjects_per_condition);



% -------------- run linear regression to find slope and intercept for RT's

empirical_RTs = empirical_stats(:, 4);
simulation_cycles = simulation_stats(:, 4);

p = polyfit(simulation_cycles, empirical_RTs, 1);
RT_slope = p(1);
RT_intercept = p(2);
yfit = polyval(p, simulation_cycles);

yresid = empirical_RTs - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(empirical_RTs)-1) * var(empirical_RTs);
rsq = 1 - SSresid/SStotal;


OG_RT_label_cycles_to_msec = sprintf('OG RT (msec = cycles * %.1f + %.1f)', RT_slope, RT_intercept);

if DO_PLOT
    figure;
    scatter(simulation_cycles, empirical_RTs, 'fill');
    clear xlabel ylabel;
    xlabel('Simulation RTs (cycles)');
    ylabel('Empirical RTs (msec)');
    text(73.4, 7400, OG_RT_label_cycles_to_msec, 'fontsize', 14);
    lsline
    title(sprintf('R^2 = %.4f', rsq));
    
    h = get(gca, 'xlabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'ylabel');
    set(h, 'FontSize', 15);
    h = get(gca, 'title');
    set(h, 'FontSize', 15);    
end


%{

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
    
    subplot(3, 2, 1);
    barweb([93 61], [16 32]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
    legend({'Focal', 'Nonfocal'});
    ylim([30 100]);

    subplot(3, 2, 2);
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
    samples = blocks(blocks(:, 1) == 0 & blocks(:, 2) == FOCAL, :)
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

    [p, table] = anovan(samples(:,4), {samples(:, 1), samples(:, 10)}, 'model','full', 'display', 'off');
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
% ----------------------------- TABLE -------------
% ---------------------------------------------------

%}

if DO_PLOT
    figure;

    subplot(3, 2, 1);
    title('Empirical Data');
    ylabel('OG RT (msec)');
    plot_all_conditions_exp4(empirical_stats(:, [1 12 3 4 5]), 6000, 8000, 1, 0, false);

    subplot(3, 2, 2);
    title('Simulation Data');
    %ylabel(OG_RT_label_cycles_to_msec);
    plot_all_conditions_exp4(simulation_stats(:, [1 12 3 4 5]), 6000, 8000, RT_slope, RT_intercept, true);

    subplot(3, 2, 3);
    ylabel('OG Accuracy (%)');
    plot_all_conditions_exp4(empirical_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 4);
    plot_all_conditions_exp4(simulation_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 5);
    ylabel('PM Hit Rate (%)');
    plot_all_conditions_exp4(empirical_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);

    subplot(3, 2, 6);
    plot_all_conditions_exp4(simulation_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);
end
