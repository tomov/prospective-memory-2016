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
 (see EM2005)
%}

DO_PLOT = false;
subjects = data;

% note that you can run this with the data returned from experiment 2
% it will just ignore the blocks and the subject id's


% -------------- define the empirical stats (Table 1 from E&M 2005)

%{
 stats variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
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

SD_cols = [5,7,9,11]; % SDs. we convert those to SEM by dividing by subjects_per_condition
subjects_per_condition = 24;

empirical_stats = [
    1 1 0, 1073.25, 112.04, 97, 2, NaN, NaN, NaN, NaN;  % no-PM, focal,    low emph
    0 1 0, 1120.87, 116.48, 97, 2, NaN, NaN, 88, 16;    % PM,    focal,    low emph
    1 1 1, 1149.25, 137.58, 97, 2, NaN, NaN, NaN, NaN;  % no-PM, focal,    high emph
    0 1 1, 1239.17, 175.42, 97, 2, NaN, NaN, 92, 16;    % PM,    focal,    high emph
    1 0 0, 1140.92, 172.87, 97, 2, NaN, NaN, NaN, NaN;  % no-PM, nonfocal, low emph
    0 0 0, 1425.39, 378.52, 97, 2, NaN, NaN, 53, 34;    % PM,    nonfocal, low emph
    1 0 1, 1183.17, 164.43, 97, 2, NaN, NaN, NaN, NaN;  % no-PM, nonfocal, high emph
    0 0 1, 1593.43, 300.86, 97, 2, NaN, NaN, 81, 27;    % PM,    nonfocal, high emph 
];


% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);


% ------------- calculate simulation stats (Table 1 from E&M 2005)

simulation_stats = [];
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        for OG_ONLY = 1:-1:0            
            stat = [OG_ONLY, FOCAL, EMPHASIS];
            for col = 4:7
                samples = subjects(subjects(:, 1) == OG_ONLY & subjects(:, 2) == FOCAL & subjects(:, 3) == EMPHASIS, col);
                samples = samples(~isnan(samples));
                M = mean(samples);
                SD = std(samples);
                %assert(length(samples) == subjects_per_condition);
                stat = [stat, M, SD];
            end
            simulation_stats = [simulation_stats; stat];
        end
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

OG_RT_label_cycles_to_msec = sprintf('OG RT (msec) = cycles * %.1f + %.1f', RT_slope, RT_intercept);

% from simulation 4
% TODO is this okay? assert(size(unique(subjects(:, [2 3 9]), 'rows'), 1) == 1); % this approach only works if we have a single condition; for reference, see sanity check at the end of EM2005.m
RT_costs = subjects(subjects(:, 1) == 0, 4) - subjects(subjects(:, 1) == 1, 4);



if DO_PLOT
    figure;
    scatter(simulation_cycles, empirical_RTs, 'fill');
    clear xlabel ylabel;
    xlabel('Simulation RTs (cycles)');
    ylabel('Empirical RTs (msec)');
    lsline
    text(86, 1550, OG_RT_label_cycles_to_msec, 'fontsize', 12);
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


% ---------- OG accuracy -- PM vs No-PM task

[p, table] = anovan(subjects(:, 5), {subjects(:,1)}, 'model','full', 'display', 'off');
fprintf('\n\n----- OG accuracy -- PM vs. No-PM; non-significant in E&M 2005 -----\n');
fprintf('                 F = %.4f, p = %f ; %.4f vs. %.4f\n', table{2,6}, p(1), mean(subjects(subjects(:, 1) == 1, 5)), mean(subjects(subjects(:, 1) == 0, 5)));

% ---------- OG accuracy -- Focality

fprintf('\n\n----- OG accuracy -- Focality; non-significant in E&M 2005 -----\n');
[p, table] = anovan(subjects(:, 5), {subjects(:,2)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f ; %.4f vs. %.4f\n', table{2,6}, p(1), mean(subjects(subjects(:, 2) == 1, 5)), mean(subjects(subjects(:, 2) == 0, 5)));


% ---------- OG accuracy -- Emphasis

fprintf('\n\n----- OG accuracy -- Emphasis; non-significant in E&M 2005? -----\n');
[p, table] = anovan(subjects(:, 5), {subjects(:,3)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f; %.4f vs. %.4f\n', table{2,6}, p(1), mean(subjects(subjects(:, 3) == 0, 5)), mean(subjects(subjects(:, 3) == 1, 5)));

% ---------- OG accuracy -- PM vs No-PM task x Focality interaction

% TODO am I doing this right? why are these different than the ones above?
fprintf('\n\n----- OG accuracy -- PM vs. No-PM x Focality interaction; non-significant in E&M 2005 -----\n');
[p, table] = anovan(subjects(:, 5), {subjects(:,1) subjects(:,2)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f\n', table{4,6}, p(3));

fprintf('\n\n----- OG accuracy -- Focality insignificant in OG-only -----\n');
[p, table] = anovan(subjects(subjects(:, 1) == 1, 5), {subjects(subjects(:, 1) == 1, 2)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


% ---------- OG accuracy -- PM vs No-PM task x Emphasis interaction

fprintf('\n\n----- OG accuracy -- PM vs. No-PM x Emphasis interaction; non-significant in E&M 2005 -----\n');
[p, table] = anovan(subjects(:, 5), {subjects(:,1) subjects(:,3)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f\n', table{4,6}, p(3));

fprintf('\n\n----- OG accuracy -- Emphasis insignificant in OG-only -----\n');
[p, table] = anovan(subjects(subjects(:, 1) == 1, 5), {subjects(subjects(:, 1) == 1, 3)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


% ---------- OG accuracy -- PM vs No-PM task x Focality x Emphasis interaction

fprintf('\n\n----- OG accuracy -- PM vs. No-PM x Focality x Emphasis interaction; non-significant in E&M 2005 -----\n');
[p, table] = anovan(subjects(:, 5), {subjects(:,1) subjects(:,2) subjects(:,3)}, 'model','full', 'display', 'off');
fprintf('                 F = %.4f, p = %f\n', table{8,6}, p(7));
fprintf('nonfocal, high emph %.4f vs. focal, low emph %.4f vs. no-PM %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 1, 5)), ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 0, 5)), ...
    mean(subjects(subjects(:, 1) == 1, 5)));





% ---------------------------------------------------
% ----------------------------- PM PERFORMANCE ------
% ---------------------------------------------------




PM_hit = subjects(:, 7);

% ----------------------- PM hit rate in focal vs. nonfocal ----

PM_hit_focal = PM_hit(subjects(:, 2) == 1);
PM_hit_nonfocal = PM_hit(subjects(:, 2) == 0);
[p, table] = anova1([PM_hit_focal PM_hit_nonfocal], {'Focal', 'Nonfocal'}, 'off');

fprintf('\n\n----- PM Performance: Focal vs. Nonfocal ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 20.03\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,5}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for FOCAL = 1:-1:0
        samples = subjects(subjects(:, 1) == 0 & subjects(:, 2) == FOCAL, 7);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        Ms(2 - FOCAL) = M;
        SEMs(2 - FOCAL) = SEM;
        fprintf('focal = %d, %.3f +- %.3f\n', FOCAL, M, SEM);
    end
    
    figure;
    
    subplot(3, 2, 1);
    barweb([90 67], [16 33]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
    h = legend({'Foc', 'Nonfoc'});
    set(h, 'FontSize', 10);
    ylim([30 100]);

    subplot(3, 2, 2);
    barweb(Ms, SEMs, 1, {}, ...
        'Simulation Data', 'PM Condition');
    ylim([30 100]);
end



% --------------------- PM hit rate in high emphasis vs. low emphasis ----

PM_hit_low = PM_hit(subjects(:, 3) == 0);
PM_hit_high = PM_hit(subjects(:, 3) == 1);
[p, table] = anova1([PM_hit_high PM_hit_low], {'High', 'Low'}, 'off');

fprintf('\n\n----- PM Performance: High vs. Low Emphasis ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 10.41\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,5}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for EMPHASIS = 0:1
        samples = subjects(subjects(:, 1) == 0 & subjects(:, 3) == EMPHASIS, 7);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        Ms(EMPHASIS + 1) = M;
        SEMs(EMPHASIS + 1) = SEM;
        fprintf('emphasis = %d, %.3f +- %.3f\n', EMPHASIS, M, SEM);
    end
    
    subplot(3, 2, 3);
    barweb([70 87], [32 22]/sqrt(subjects_per_condition), 1, {}, ...
        '', 'PM Condition', 'PM Hit rate (%)');
    h =legend({'Low', 'High'});
    set(h, 'FontSize', 10);
    ylim([30 100]);

    subplot(3, 2, 4);
    barweb(Ms, SEMs, 1, {}, ...
        '', 'PM Condition');
    ylim([30 100]);
end


% -------------- PM hit rate in high emph vs. low emph. for different focalities 

[p, table] = anovan(PM_hit, {subjects(:, 2) subjects(:, 3)}, 'model','full', 'display', 'off');

fprintf('\n\n----- PM Performance: Interaction between Focality and Emphasis ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('    Focality:    F = 20.03\n');
fprintf('    Emphasis:    F = 10.41\n');
fprintf('    Interaction: F = 5.73\n');
fprintf('\n  Simulation Data -------\n');
fprintf('    Focality:    F = %.4f, p = %f\n', table{2,6}, p(1));
fprintf('    Emphasis:    F = %.4f, p = %f\n', table{3,6}, p(2));
fprintf('    Interaction: F = %.4f, p = %f\n', table{4,6}, p(3));

if DO_PLOT
    titles = {'', ''};
    sources = {empirical_stats, simulation_stats};
    for s_id = 1:2
        subplot(3, 2, s_id + 4);

        stats = sources{s_id};
        Ms = zeros(2);
        SEMs = zeros(2);
        for FOCAL = 1:-1:0
            for EMPHASIS = 0:1
                M = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 10);
                SEM = stats(stats(:,1) == 0 & ...
                    stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 11);
                Ms(2 - FOCAL, EMPHASIS + 1) = M;
                SEMs(2 - FOCAL, EMPHASIS + 1) = SEM;
                fprintf('focal = %d, emphasis = %d, %.3f +- %.3f\n', FOCAL, EMPHASIS, M, SEM);
            end
        end

        barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
            titles{s_id}, 'PM Condition');
        if s_id == 1
            h = legend({'Low', 'High'});
            set(h, 'FontSize', 10);
            ylabel('PM Hit rate (%)');
        end
        ylim([30 100]);
    end
end










%%% ----------------- PM RT's ----------------------------------

PM_RTs = subjects(:, 6);

[p, table] = anovan(PM_RTs, {subjects(:, 2) subjects(:, 3)}, 'model','full', 'display', 'off');

fprintf('\n\n----- PM RTs NEW: Interaction between Focality and Emphasis ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('    Focality:    F = %.4f, p = %f\n', table{2,6}, p(1));
fprintf('    Emphasis:    F = %.4f, p = %f\n', table{3,6}, p(2));
fprintf('    Interaction: F = %.4f, p = %f\n', table{4,6}, p(3));

if DO_PLOT
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
end








% ---------------------------------------------------
% ----------------------------- OG PERFORMANCE ------
% ---------------------------------------------------


% ------------------ OG accuracy

OG_hit = subjects(:, 5);

[p, table] = anovan(OG_hit, {subjects(:, 1) subjects(:, 2) subjects(:, 3)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG accuracy: 2x2x2 ANOVA ------\n');
table(1:8,6)
p
fprintf('===> shit, some of these are significant... oh well....\n');



% ----------------- OG RT: focal vs. nonfocal

OG_RTs = subjects(:, 4);

[p, table] = anovan(OG_RTs, {subjects(:, 2)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: Focal vs. Nonfocal ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 22.87\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

fprintf('  ====>  %.4f +- %.4f vs. %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 2) == 1, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 2) == 1, 4)) * RT_slope / sqrt(subjects_per_condition), ...
    mean(subjects(subjects(:, 2) == 0, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 2) == 0, 4)) * RT_slope / sqrt(subjects_per_condition) ...
    );

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for FOCAL = 1:-1:0
        samples = subjects(subjects(:, 2) == FOCAL, 4);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        M = M * RT_slope + RT_intercept;
        SEM = SEM * RT_slope;
        Ms(2 - FOCAL) = M;
        SEMs(2 - FOCAL) = SEM;
        fprintf('focal = %d, %.3f +- %.3f\n', FOCAL, M, SEM);
    end
    
    figure;
    
    subplot(3, 2, 1);
    barweb([1145.63 1335.73], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'Ongoing RT (msec)');
    h = legend({'Focal', 'Nonfoc'});
    set(h, 'FontSize', 10);
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

fprintf('  ====>  %.4f +- %.4f vs. %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 3) == 1, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 3) == 1, 4)) * RT_slope  / sqrt(subjects_per_condition), ...
    mean(subjects(subjects(:, 3) == 0, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 3) == 0, 4)) * RT_slope  / sqrt(subjects_per_condition) ...
    );

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    for EMPHASIS = 0:1
        samples = subjects(subjects(:, 3) == EMPHASIS, 4);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples));
        M = M * RT_slope + RT_intercept;
        SEM = SEM * RT_slope;
        Ms(EMPHASIS + 1) = M;
        SEMs(EMPHASIS + 1) = SEM;
        fprintf('emphasis = %d, %.3f +- %.3f\n', EMPHASIS, M, SEM);
    end
    
    subplot(3, 2, 3);
    barweb([1190.11 1291.26], [0 0]/sqrt(subjects_per_condition), 1, {}, ...
        '', 'PM Condition', 'Ongoing RT (msec)');
    h = legend({'Low', 'High'});
    set(h, 'FontSize', 10);
    ylim([1000 1400]);

    subplot(3, 2, 4);
    barweb(Ms, SEMs, 1, {}, ...
        '', 'PM Condition');
    ylim([1000 1400]);
end


% ----------------- OG RT: PM vs. No PM



[p, table] = anovan(OG_RTs, {subjects(:, 1)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: PM vs. No PM ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 131.66\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

fprintf('  ====>  %.4f +- %.4f vs. %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 1) == 0, 4)) * RT_slope  / sqrt(subjects_per_condition), ...
    mean(subjects(subjects(:, 1) == 1, 4)) * RT_slope + RT_intercept, ...
    std(subjects(subjects(:, 1) == 1, 4)) * RT_slope  / sqrt(subjects_per_condition) ...
    );

for OG_ONLY = 0:1
    samples = subjects(subjects(:, 1) == OG_ONLY, 4);
    M = mean(samples);
    SD = std(samples);
    SEM = SD / sqrt(length(samples));
    M = M * RT_slope + RT_intercept;
    SEM = SEM * RT_slope;
    Ms(EMPHASIS + 1) = M;
    SEMs(EMPHASIS + 1) = SEM;
    fprintf('og only = %d, %.3f +- %.3f\n', OG_ONLY, M, SEM);
end



% ----------------- OG RT: 2x2x2 ANOVA


[p, table] = anovan(OG_RTs, {subjects(:, 1) subjects(:, 2) subjects(:, 3)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: 2x2x2 ANOVA ------\n');
table(1:8,6)
p
fprintf('===> shit...theyre all significant.. .fuck fuck fuck\n');

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
    titles = {'', ''};
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
                M_ogonly = stats(stats(:,1) == 1 & stats(:, 2) == FOCAL & stats(:, 3) == EMPHASIS, 4);
                M_ogonly = M_ogonly * RT_slope + RT_intercept; 
                fprintf('focal = %d, emphasis = %d, %.3f +- %.3f, cost = %.3f\n', FOCAL, EMPHASIS, M, SEM, M - M_ogonly);
            end
        end

        barweb(Ms, SEMs, 1, {'Focal', 'Nonfocal'}, ...
            titles{s_id}, 'PM Condition');
        if s_id == 1
            h = legend({'Low', 'High'});
            set(h, 'FontSize', 10);
            ylabel('Ongoing RT (msec)');
        end
        ylim([1000 1700]);
    end
end


fprintf('\n\n');
fprintf('OG accuracy in ALL = %.4f +- %.4f\n', ...
    mean(subjects(:, 5)), ...
    std(subjects(:, 5)) / sqrt(subjects_per_condition));

fprintf('\n\n');
fprintf('OG RT cost in Focal, low emphasis = %.4f\n', mean(RT_costs(subjects(:, 2) == 1 & subjects(:, 3) == 0 & subjects(:, 1) == 0)) * RT_slope);
fprintf('OG RT cost in Focal, high emphasis = %.4f\n', mean(RT_costs(subjects(:, 2) == 1 & subjects(:, 3) == 1 & subjects(:, 1) == 0)) * RT_slope);
fprintf('OG RT cost in Nonfocal, low emphasis = %.4f\n', mean(RT_costs(subjects(:, 2) == 0 & subjects(:, 3) == 0 & subjects(:, 1) == 0)) * RT_slope);
fprintf('OG RT cost in Nonfocal, high emphasis = %.4f\n', mean(RT_costs(subjects(:, 2) == 0 & subjects(:, 3) == 1 & subjects(:, 1) == 0)) * RT_slope);


fprintf('\n\n');
fprintf('PM hit rate in Focal = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in Nonocal = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in High Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 3) == 1, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 3) == 1, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in Low Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 3) == 0, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 3) == 0, 7)) / sqrt(subjects_per_condition));

fprintf('\n\n');
fprintf('PM hit rate in Focal, Low Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 0, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 0, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in Focal, High Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 1, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 1, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in Nonfocal, Low Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 0, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 0, 7)) / sqrt(subjects_per_condition));
fprintf('PM hit rate in Nonfocal, High Emph = %.4f +- %.4f\n', ...
    mean(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 1, 7)), ...
    std(subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 1, 7)) / sqrt(subjects_per_condition));


fprintf('\n\n');
fprintf('OG RT in ALL = %.4f +- %.4f\n', ...
    nanmean(subjects(:, 4)) * RT_slope + RT_intercept, ...
    nanstd(subjects(:, 4)) * RT_slope / sqrt(subjects_per_condition));


fprintf('\n\n');
fprintf('PM RT in ALL = %.4f +- %.4f\n', ...
    nanmean(subjects(:, 6)) * RT_slope + RT_intercept, ...
    nanstd(subjects(:, 6)) * RT_slope / sqrt(subjects_per_condition));

% ------------- PM RT vs. OG RT
PM_RT = subjects(~isnan(subjects(:, 6)), 6);
OG_RT = subjects(:, 4);
%[p, table] = anova1([PM_RT OG_RT], {'PM RT', 'OG RT'}, 'off');
[p, table] = anovan([PM_RT; OG_RT], {[ones(length(PM_RT), 1); zeros(length(OG_RT), 1)]}, 'model','full', 'display', 'off');
fprintf('\n\n----- PM RT vs. OG OT ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


fprintf('\n\n');
fprintf('PM RT in Focal = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in Nonocal = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in High Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 3) == 1, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 3) == 1, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in Low Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 3) == 0, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 3) == 0, 6)) / sqrt(subjects_per_condition));

fprintf('\n\n');
fprintf('PM RT in Focal, Low Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 0, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 0, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in Focal, High Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 1, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 1 & subjects(:, 3) == 1, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in Nonfocal, Low Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 0, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 0, 6)) / sqrt(subjects_per_condition));
fprintf('PM RT in Nonfocal, High Emph = %.4f +- %.4f\n', ...
    nanmean(RT_intercept + RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 1, 6)), ...
    nanstd(RT_slope * subjects(subjects(:, 1) == 0 & subjects(:, 2) == 0 & subjects(:, 3) == 1, 6)) / sqrt(subjects_per_condition));


% ------------- PM RT in focal vs. nonfocal
PM_RT = subjects(~isnan(subjects(:, 6)), 6);
is_foc = subjects(~isnan(subjects(:, 6)), 2);
[p, table] = anovan(PM_RT, {is_foc}, 'model','full', 'display', 'off');
fprintf('\n\n----- PM RT in focal vs. nonfocal. ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

% ------------- PM RT in high vs. low emphasis
PM_RT = subjects(~isnan(subjects(:, 6)), 6);
is_high_emph = subjects(~isnan(subjects(:, 6)), 3);
[p, table] = anovan(PM_RT, {is_high_emph}, 'model','full', 'display', 'off');
fprintf('\n\n----- PM RT in high vs. low emphasis. ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));


% ------------- PM RT -- interaction
PM_RT = subjects(~isnan(subjects(:, 6)), 6);
is_high_emph = subjects(~isnan(subjects(:, 6)), 3);
is_foc = subjects(~isnan(subjects(:, 6)), 2);
[p, table] = anovan(PM_RT, {is_high_emph is_foc}, 'model','full', 'display', 'off');
fprintf('\n\n----- PM RT interaction ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{4,6}, p(3));


% ------------- Intention superiority effect
PM_miss_OG_RT = subjects(~isnan(subjects(:, 11)), 11);
OG_RT = subjects(:, 4);
[p, table] = anovan([PM_miss_OG_RT; OG_RT], {[ones(length(PM_miss_OG_RT), 1); zeros(length(OG_RT), 1)]}, 'model','full', 'display', 'off');
fprintf('\n\n----- Intention superiority effect ------\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));





% ---------------------------------------------------
% ----------------------------- TABLE 1 -------------
% ---------------------------------------------------





if DO_PLOT
    figure;

    subplot(3, 2, 1);
    title('Empirical Data');
    ylabel('OG RT (msec)');
    plot_all_conditions_exp1(empirical_stats(:, [1:3 4 5]), 1000, 1700, 1, 0, true);

    subplot(3, 2, 2);
    title('Simulation Data');
    %ylabel(OG_RT_label_cycles_to_msec);
    plot_all_conditions_exp1(simulation_stats(:, [1:3 4 5]), 1000, 1700, RT_slope, RT_intercept, false);

    subplot(3, 2, 3);
    ylabel('OG Accuracy (%)');
    plot_all_conditions_exp1(empirical_stats(:, [1:3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 4);
    plot_all_conditions_exp1(simulation_stats(:, [1:3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 5);
    ylabel('PM Hit Rate (%)');
    plot_all_conditions_exp1(empirical_stats(:, [1:3 10 11]), 40, 100, 1, 0, false);

    subplot(3, 2, 6);
    plot_all_conditions_exp1(simulation_stats(:, [1:3 10 11]), 40, 100, 1, 0, false);
end
