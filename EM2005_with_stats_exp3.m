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
 (see EM2005 exp 2)
%}

DO_PLOT = false;
subjects = data;


% -------------- define the empirical stats (Table 2 from E&M 2005)

%{
 stats variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS  (N/A in this case, just to make compatible with exp 1
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

subjects_per_condition = 32;

empirical_stats = [
    1 1 0, 4791, 618,        70, 11, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    low emph, 1 targets
    1 1 0, 4890, 508,        70, 11, NaN, NaN, NaN, NaN, 6;  % no-PM, focal,    low emph, 6 targets
    0 1 0, 4885, 591,        69, 10, NaN, NaN, 80,   28, 1;  % PM, focal,    low emph, 1 targets
    0 1 0, 5215, 422,        69, 10, NaN, NaN, 72,   25, 6;  % PM, focal,    low emph, 6 targets
];

% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);




% ------------- calculate simulation stats

simulation_stats = [];
EMPHASIS = 0;
FOCAL = 1;
EMPHASIS = 0;
% order here matters -- must be same as empirical_data above for line
% regression
for OG_ONLY = 1:-1:0
    for TARGETS = [1,6]
        stat = [OG_ONLY, FOCAL, EMPHASIS];
        for col = 4:7
            samples = subjects(subjects(:, 1) == OG_ONLY & subjects(:, 9) == TARGETS, col);
            samples = samples(~isnan(samples));
            M = mean(samples);
            SD = std(samples);
            % TODO B/C WE merge high and low emphasis
            %assert(length(samples) == subjects_per_condition);
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
    text(73, 5150, OG_RT_label_cycles_to_msec, 'fontsize', 12);
    lsline
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



PM_hit = subjects(:, 7);

% ----------------------- PM hit rate in focal vs. nonfocal ----

[p, table] = anovan(PM_hit, {subjects(:, 9)}, 'model','full', 'display', 'off');

fprintf('\n\n----- PM Performance: 1 target vs. 6 targets ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 1.38, p = 0.24\n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));

if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    targets_idx = 0;
    for TARGETS = [1,6]
        targets_idx = targets_idx + 1;
        samples = subjects(subjects(:, 1) == 0 & subjects(:, 9) == TARGETS, 7);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples) / 2); % TODO notice / 2 => b/c we technically are merging subjects in high/low emphasis conditions
        Ms(targets_idx) = M;
        SEMs(targets_idx) = SEM;
    end
    
    figure;
    
    subplot(3, 2, 1);
    barweb([80 72], [28 25]/sqrt(subjects_per_condition), 1, {}, ...
        'Empirical Data', 'PM Condition', 'PM Hit rate (%)');
    ylim([40 100]);

    subplot(3, 2, 2);
    barweb(Ms, SEMs, 1, {}, ...
        'Simulation Data', 'PM Condition');
    legend({'1 Target', '6 Targets'});
    ylim([40 100]);
end


% ---------------------------------------------------
% ----------------------------- OG PERFORMANCE ------
% ---------------------------------------------------


% ------------------ OG accuracy


OG_hit = subjects(:, 5);

[p, table] = anovan(OG_hit, {subjects(:, 1) subjects(:, 9)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG accuracy: 2x2 ANOVA ------\n');
table(1:4,6)
p
fprintf('===> shit, some of these are significant... oh well....\n');


if DO_PLOT
    Ms = zeros(1, 2);
    SEMs = zeros(1, 2);
    targets_idx = 0;
    for TARGETS = [1,6]
        targets_idx = targets_idx + 1;
        samples = subjects(subjects(:, 1) == 0 & subjects(:, 9) == TARGETS, 5);
        M = mean(samples);
        SD = std(samples);
        SEM = SD / sqrt(length(samples) / 2); % TODO notice / 2 => b/c we technically are merging subjects in high/low emphasis conditions
        Ms(targets_idx) = M;
        SEMs(targets_idx) = SEM;
    end
    
    subplot(3, 2, 3);
    barweb([69 70], [10 11]/sqrt(subjects_per_condition), 1, {}, ...
        '', 'PM Condition', 'OG Accuracy (%)');
    ylim([40 100]);
    
    subplot(3, 2, 4);
    barweb(Ms, SEMs, 1, {}, ...
        '', 'PM Condition');
    legend({'1 Target', '6 Targets'});
    ylim([40 100]);
end

% ----------------- OG RTs: 2x2 ANOVA


OG_RTs = subjects(:, 4);

[p, table] = anovan(OG_RTs, {subjects(:, 1), subjects(:, 9)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: 2x2 ANOVA ------\n');
table(1:4,6)
p
fprintf('===> all are significant.. wtf does that mean \n');


% ----------------- OG RT: PM vs. No PM


[p, table] = anovan(OG_RTs, {subjects(:, 1)}, 'model','full', 'display', 'off');

fprintf('\n\n----- OG RTs: PM vs. No PM ------\n');
fprintf('\n  Empirical Data -------\n');
fprintf('                 F = 27.36, \n');
fprintf('\n  Simulation Data -------\n');
fprintf('                 F = %.4f, p = %f\n', table{2,6}, p(1));



% ----------------- OG RT: 1 target vs. 6 targes, PM vs. No PM


if DO_PLOT
    titles = {'', ''};
    sources = {empirical_stats, simulation_stats};
    for s_id = 1:2
        subplot(3, 2, s_id + 4);

        stats = sources{s_id};
        Ms = zeros(2);
        SEMs = zeros(2);
        targets_idx = 0;
        for TARGETS = [1,6]
            targets_idx = targets_idx + 1;
            for OG_ONLY = 1:-1:0
                M = stats(stats(:,1) == OG_ONLY & ...
                    stats(:, 12) == TARGETS, 4);
                SEM = stats(stats(:,1) == OG_ONLY & ...
                    stats(:, 12) == TARGETS, 5);
                if s_id == 2
                    M = M * RT_slope + RT_intercept;
                    SEM = SEM * RT_slope;
                end
                Ms(2 - OG_ONLY, targets_idx) = M;
                SEMs(2 - OG_ONLY, targets_idx) = SEM;
                
                if s_id == 2
                    fprintf('targets = %d, og only = %d: $%.2f \\pm %.2f$\n', TARGETS, OG_ONLY,  M, SEM);
                end
            end
        end

        barweb(Ms, SEMs, 1, {'1 Target', '6 Targets'}, ...
            titles{s_id}, 'PM Condition');
        if s_id == 1
            ylabel('Ongoing RT (msec)');
        else
            legend({'No-PM', 'PM'});
        end
        ylim([4500 5500]);
    end
end

OG_6_target_slowing = simulation_stats(simulation_stats(:,1) == 0 & simulation_stats(:, 12) == 6, 4) ...
    -  simulation_stats(simulation_stats(:,1) == 1 & simulation_stats(:, 12) == 6, 4);
OG_6_target_slowing = OG_6_target_slowing * RT_slope;
fprintf('\n\n !!!! OG RT slowing in 6-target condition (must be ~300 msec) = \n   %.3f \n\n', OG_6_target_slowing);




if DO_PLOT
    figure;

    subplot(3, 2, 1);
    title('Empirical Data');
    ylabel('OG RT (msec)');
    plot_all_conditions_exp3(empirical_stats(:, [1 12 3 4 5]), 4500, 5500, 1, 0, false);

    subplot(3, 2, 2);
    title('Simulation Data');
    %ylabel(OG_RT_label_cycles_to_msec);
    plot_all_conditions_exp3(simulation_stats(:, [1 12 3 4 5]), 4500, 5500, RT_slope, RT_intercept, true);

    subplot(3, 2, 3);
    ylabel('OG Accuracy (%)');
    plot_all_conditions_exp3(empirical_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 4);
    plot_all_conditions_exp3(simulation_stats(:, [1 12 3 6 7]), 40, 100, 1, 0, false);

    subplot(3, 2, 5);
    ylabel('PM Hit Rate (%)');
    plot_all_conditions_exp3(empirical_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);

    subplot(3, 2, 6);
    plot_all_conditions_exp3(simulation_stats(:, [1 12 3 10 11]), 40, 100, 1, 0, false);
end
