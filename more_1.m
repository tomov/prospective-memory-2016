
%{
 subjects[subj_id, :] ==
 samples variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
 statistics order = 
    4 - OG_RT, 
    5 - OG_Hit, 
    6 - PM_RT, 
    7 - PM_Hit,
    8 - PM_miss_OG_hit
 (see EM2005)
%}

%EM2005_with_stats_exp1

A = [
    1159.21,  3.89;
    1223.7,  3.8;
    1420.5, 3.0;
    1597.7, 5.91
    ];

Z = [
    834.8750, 39.2831;
    905.5000, 45.2172;
    1268.4054, 55.9196;
    1321.7500, 62.0524;
    ];

figure;
stats = simulation_stats;
Ms = zeros(4,2);
SEMs = zeros(4,2);
idx = 0;
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        idx = idx + 1;
        Ms(idx, 1) = A(idx, 1);
        SEMs(idx, 1) = A(idx, 2);
        Ms(idx, 2) = Z(idx, 1);
        SEMs(idx, 2) = Z(idx, 2);
        fprintf('focal = %d, emphasis = %d, $%.1f \\pm %.1f$ ms    $%.1f \\pm %.1f$ ms\n', FOCAL, EMPHASIS, A(idx, 1), A(idx, 2), Z(idx, 1), Z(idx, 2));
    end
end

barweb(Ms, SEMs, 1, {'Foc,Lo', 'Foc,Hi', 'NFoc,Lo', 'NFoc,Hi'}, ...
    'Intention Superiority Effect', 'PM Condition');
    h = legend({'OG', 'PM miss'});
    set(h, 'FontSize', 15);
    h = ylabel('RTs (msec)');
    set(h, 'FontSize', 15);

ylim([700 1700]);








for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        ids = subjects(:, 2) == FOCAL & subjects(:, 3) == EMPHASIS;
        
        PM_miss_RT = subjects(ids, 8);
        PM_miss_RT = PM_miss_RT(~isnan(PM_miss_RT));

        OG_RT = subjects(ids, 4);
        OG_RT = OG_RT(~isnan(OG_RT));

        all = [PM_miss_RT; OG_RT];
        groups = [repmat({'PM miss'}, length(PM_miss_RT), 1); repmat({'OG'}, length(OG_RT), 1)];
        
        M_og = mean(OG_RT) * RT_slope + RT_intercept;
        SE_og = std(OG_RT) * RT_slope / sqrt(length(OG_RT));
        
        M_pm = mean(PM_miss_RT) * RT_slope + RT_intercept;
        SE_pm = std(PM_miss_RT) * RT_slope / sqrt(length(PM_miss_RT));

        [p, table] = anova1(all, groups, 'off')

        %[p, table] = anovan(samples(:,8), {samples(:, 1)}, 'model','full', 'display', 'off');
        fprintf('\n\n----- ISI: %s, %s ------\n', focal_titles{FOCAL+1}, emphasis_titles{EMPHASIS+1});
        fprintf('\n  Simulation Data -------\n');
        fprintf('                 F = %.4f, p = %f\n', table{2,5}, p(1));
        fprintf('                mean OG RT = $%.1f \\pm %.1f$ ms, mean PM miss RT = $%.1f \\pm %.1f$ ms\n', ...
            M_og, SE_og, M_pm, SE_pm);
    end
end



PM_miss_RT = subjects(:, 8);
PM_miss_RT = PM_miss_RT(~isnan(PM_miss_RT));

mean(PM_miss_RT)
std(PM_miss_RT) / sqrt(length(PM_RT_high))

mean(PM_miss_RT) * RT_slope + RT_intercept
std(PM_miss_RT) * RT_slope / sqrt(length(PM_RT_high))



OG_RT = subjects(:, 4);
OG_RT = OG_RT(~isnan(OG_RT));

all = [PM_miss_RT; OG_RT];
groups = [repmat({'PM miss'}, length(PM_miss_RT), 1); repmat({'OG'}, length(OG_RT), 1)];

[p, table] = anova1(all, groups, 'off')



%{

PM_RT = subjects(:, 6);
PM_RT_high = PM_RT(subjects(:, 3) == 1);
PM_RT_low = PM_RT(subjects(:, 3) == 0);
PM_RT_high = PM_RT_high(~isnan(PM_RT_high));
PM_RT_low = PM_RT_low(~isnan(PM_RT_low));
all = [PM_RT_high; PM_RT_low];
groups = [repmat({'PM high'}, length(PM_RT_high), 1); repmat({'PM low'}, length(PM_RT_low), 1)];

mean(PM_RT_high) * RT_slope + RT_intercept
std(PM_RT_high) * RT_slope / sqrt(length(PM_RT_high))

mean(PM_RT_low) * RT_slope + RT_intercept
std(PM_RT_low) * RT_slope / sqrt(length(PM_RT_low))


[p, table] = anova1(all, groups, 'off')

%}


%{


PM_RT = subjects(:, 6);
PM_RT_focal = PM_RT(subjects(:, 2) == 1);
PM_RT_nonfocal = PM_RT(subjects(:, 2) == 0);
PM_RT_focal = PM_RT_focal(~isnan(PM_RT_focal));
PM_RT_nonfocal = PM_RT_nonfocal(~isnan(PM_RT_nonfocal));
all = [PM_RT_focal; PM_RT_nonfocal];
groups = [repmat({'PM focal'}, length(PM_RT_focal), 1); repmat({'PM nonfocal'}, length(PM_RT_nonfocal), 1)];

mean(PM_RT_focal) * RT_slope + RT_intercept
std(PM_RT_focal) * RT_slope / sqrt(length(PM_RT_focal))

mean(PM_RT_nonfocal) * RT_slope + RT_intercept
std(PM_RT_nonfocal) * RT_slope / sqrt(length(PM_RT_nonfocal))


[p, table] = anova1(all, groups, 'off')


%}



%{

% --------------------------------------------------

PM_hit_RT = subjects(:, 6);
PM_hit_RT = PM_hit_RT(~isnan(PM_hit_RT));

OG_RT = subjects(:, 4);
OG_RT = OG_RT(~isnan(OG_RT));

all = [PM_hit_RT; OG_RT];
groups = [repmat({'PM'}, length(PM_hit_RT), 1); repmat({'OG'}, length(OG_RT), 1)];

[p, table] = anova1(all, groups, 'off')
%}

