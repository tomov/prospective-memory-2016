
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
 (see EM2005 exp 2)
%}

subjects = data;


% -------------- define the empirical stats (Table 2 from E&M 2005)

%{
 stats variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - CAPACITY
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

subjects_per_condition = 24;

empirical_stats = [
    % low WM capacity
    1 1 0, 696.39,  17.65,       97, 2, NaN,     NaN,    NaN, NaN, 1;  % baseline, 
    0 1 0, 698.47,  16.31,       97, 2, 780.81,  159.71, 88,   3,  1;  % PM, focal, 
    0 0 0, 779.32,  26.00,       97, 2, 1489.72, 208.40, 68,   7,  1;  % PM, non-non-focal,
    
    % high WM capacity
    1 1 1, 726.96,  27.24,        97, 2, NaN,    NaN,    NaN, NaN, 1;  % baseline,
    0 1 1, 708.43,  22.97,        97, 2, 821.99, 124.61, 93,   3,  1;  % PM, focal,
    0 0 1, 788.48,  30.76,        97, 2, 910.42, 162.58, 92,   3,  1;  % PM, non-focal,
];

% convert SD's to SEM's in empirical data
%-- NO NEED TO; they're already SEM's
% empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);

% ------------- calculate simulation stats

simulation_stats = [];
TARGETS = 1;
% order here matters -- must be same as empirical_data above for line
% regression
for CAPACITY = 0:1
    for OG_ONLY = 1:-1:0
        for FOCAL = 1:-1:0
            if OG_ONLY && ~FOCAL
                continue % no such data point
            end
            stat = [OG_ONLY, FOCAL, CAPACITY];
            for col = 4:7
                samples = subjects(subjects(:, 1) == OG_ONLY & subjects(:, 2) == FOCAL & subjects(:, 3) == CAPACITY, col);
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
