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
    10 - first PM RT
    11 - PM_miss_OG_RT
    12 - third task ex-target RT's
    13 - third task ex-target hits
    14 - third task ex-nontarget RT's
    15 - third task ex-nontarget hits
    16 - third task ex-target wrong answer PM hits
 (see EM2005)
%}

DO_PLOT = false;
subjects = data;


% -------------- define the empirical stats (Table 4 from E&M 2005)

%{
 stats variable order = 
    1 - INTER_TASK,        !!! different from usual
    2 - FOCAL, 
    3 - EMPHASIS
    12 - IS_TARGET_IN_INTER_TASK
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
empirical_subjects_per_condition = 72;
simulation_subjects_per_condition = 372;

% most table 4 from E&M
%
empirical_stats = [
    0 1 0,   NaN,   NaN, NaN, NaN, NaN, NaN,  76,  26, NaN;  % not much info on the OG + PM task (word rating in E&M)
    1 1 0,  576.50,  87.43, 96, 2, NaN, NaN, NaN, NaN, 0;  % inter task, non-targets (i.e. "previously presented" items in E&M)
    1 1 0,  631.13, 133.19, 96, 2, NaN, NaN, NaN, NaN, 1;  % inter task (lexical decision task in E&M), target items
];

% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(empirical_subjects_per_condition);


% ------------- calculate simulation stats (Table 1 from E&M 2005)

simulation_stats = [];

FOCAL = 1;
EMPHASIS = 0;

% First put the OG + PM stats
%
INTER_TASK = 0;
IS_TARGET_IN_INTER_TASK = NaN;
stat = [INTER_TASK, FOCAL, EMPHASIS];
for col = 4:7
    samples = subjects(:, col);
    samples = samples(~isnan(samples));

    M = mean(samples);
    SD = std(samples);
    stat = [stat, M, SD];
end
stat = [stat, IS_TARGET_IN_INTER_TASK];
simulation_stats = [simulation_stats; stat];


% Then put the Inter task stats
%
which_cols = [
    14 15 6 7; % nontarget items in Inter task
    12 13 6 7; % target items in Inter task
];
INTER_TASK = 1;
for IS_TARGET_IN_INTER_TASK = 0:1
    cols = which_cols(IS_TARGET_IN_INTER_TASK + 1, :);
    stat = [INTER_TASK, FOCAL, EMPHASIS];
    for col = cols
        samples = subjects(:, col);
        samples = samples(~isnan(samples));
        M = mean(samples);
        SD = std(samples);
        %assert(length(samples) == subjects_per_condition);
        stat = [stat, M, SD];
    end
    stat = [stat, IS_TARGET_IN_INTER_TASK];
    simulation_stats = [simulation_stats; stat];
end


% convert SD's to SEM's in simulation data
simulation_stats(:, SD_cols) = simulation_stats(:, SD_cols) / sqrt(simulation_subjects_per_condition);




% -------------- run linear regression to find slope and intercept for RT's

% TODO this is useless -- there are only 2 data points...

empirical_RTs = empirical_stats([2 3], 4);
simulation_cycles = simulation_stats([2 3], 4);

p = polyfit(simulation_cycles, empirical_RTs, 1);
RT_slope = p(1);
RT_intercept = p(2);
yfit = polyval(p, simulation_cycles);

yresid = empirical_RTs - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(empirical_RTs)-1) * var(empirical_RTs);
rsq = 1 - SSresid/SStotal;

OG_RT_label_cycles_to_msec = sprintf('OG RT (msec) = cycles * %.1f + %.1f', RT_slope, RT_intercept)


% do some analysis


Ms = [];
SEMs = [];


% from experiment 1
RT_slope = 8.9868;
RT_intercept = 360;

GroupA = subjects(:, 12); 

m = mean(GroupA) * RT_slope + RT_intercept
s = std(GroupA) * RT_slope / sqrt(length(GroupA))

Ms = [Ms; m];
SEMs = [SEMs; s];

GroupB = subjects(:, 14);

m = mean(GroupB) * RT_slope + RT_intercept
s = std(GroupB) * RT_slope / sqrt(length(GroupB))

Ms = [Ms; m];
SEMs = [SEMs; s];

all = [GroupA; GroupB];
groups = [repmat({'1'}, length(GroupA), 1); repmat({'4'}, length(GroupB), 1)];

[p, table] = anova1(all, groups, 'off')

%{
figure;
barweb(Ms, SEMs, 1, {}, ...
    'Simulation Data', 'Third Task Trial Type');
h = legend({'Target', 'Nontarget'});
set(h, 'FontSize', 15);
ylabel('Third Task RT (ms)');
ylim([1100, 1200]);
 %}       


