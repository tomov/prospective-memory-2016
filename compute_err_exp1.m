function error = compute_err_exp1(data, extra) % takes in the output of EM2005
% Computes the error for a particular simulation of EM2005

% get stats from subjects and analyze so you can fit with fitparam()
%{
 data[subj_id, :] ==
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

%%% CONFIG %%%
M_cols = [4,6,8,10]; % means. we use those to calculate the error
SD_cols = [5,7,9,11]; % SDs. we convert those to SEM by dividing by subjects_per_condition
RT_mean_cols = [4,8]; % RTs. we regress those against model cycles means to scale
RT_SD_cols = [5,9]; % RT SDs. we regress those against model cycles SDs to scale

use_M_cols = M_cols; % use all means (OG RT, OG accuracy, PM RT, PM hit rates) (notice we ignore the ones we don't know i.e. the NaNs)
use_SD_cols = [5 9 11]; % use SD for OG RT, PM RT and PM hit rates (we don't know the OG accuracy SD's)
use_cols    = [use_M_cols use_SD_cols]; % which columns to use in the error computation
err_weights = [2 2 2 2 1 1 1];  % means weight more than SD's when fitting (keep in mind we're scaling them to %'s so they're all on the same scale by default)
assert(length(use_cols) == length(err_weights));

subjects_per_condition = 24;

% comes from E&M 2005 paper
empirical_stats = [
   %       OG RT's          OG hit PM RT     PM hit
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

% resize weights of errors to # of conditions
err_scalers = repmat(err_weights, size(empirical_stats, 1), 1);
diff_err_scalers = [3; 3];

%
% ------------- calculate simulation stats (Table 1 from E&M 2005)
%

simulation_stats = [];
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        for OG_ONLY = 1:-1:0            
            stat = [OG_ONLY, FOCAL, EMPHASIS];
            for col = 4:7
                samples = data(data(:, 1) == OG_ONLY & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS, col);
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


%
% -------------- run linear regression to find slope and intercept for RTs
%

empirical_RTs = empirical_stats(:, 4);
simulation_cycles = simulation_stats(:, 4);

p = polyfit(simulation_cycles, empirical_RTs, 1);
RT_slope = p(1);
RT_intercept = p(2);
yfit = polyval(p, simulation_cycles);

% replace model cycles with fit RT values
simulation_stats(:,RT_mean_cols) = polyval(p, simulation_stats(:,RT_mean_cols));
simulation_stats(:,RT_SD_cols) = simulation_stats(:,RT_SD_cols) * RT_slope;


% ----------------- Cohen's D -- not sure how to use it though...
%
%focalD = @(stats) CohensD(stats(2, 4), stats(2, 5) * sqrt(data_per_condition), data_per_condition, stats(4, 4), stats(4, 5) * sqrt(data_per_condition), data_per_condition);
%focalD(empirical_stats)
%focalD(simulation_stats)

% --------------- measure the effect significance -- fit that too
%
empirical_Fs = [];
simulation_Fs = [];

PM_hit = data(:, 7); % OG_ONLY (i.e. NaN's) are automatically ignored
OG_hit = data(:, 5);
OG_RTs = data(:, 4);
% ----------------------- PM hit rate in focal vs. nonfocal ----
PM_hit_focal = PM_hit(data(:, 2) == 1);
PM_hit_nonfocal = PM_hit(data(:, 2) == 0);
[p, table] = anova1([PM_hit_focal PM_hit_nonfocal], {'Focal', 'Nonfocal'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 20.03];
% --------------------- PM hit rate in high emphasis vs. low emphasis ----
PM_hit_low = PM_hit(data(:, 3) == 0);
PM_hit_high = PM_hit(data(:, 3) == 1);
[p, table] = anova1([PM_hit_high PM_hit_low], {'High', 'Low'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 10.41];

% -------------- PM hit rate in high emph vs. low emph. for different focalities 
[p, table] = anovan(PM_hit, {data(:, 2) data(:, 3)}, 'model','interaction', 'display', 'off');
simulation_Fs = [simulation_Fs table{2, 6} table{3, 6} table{4, 6}];
empirical_Fs = [empirical_Fs 20.03 10.41 5.73];
%%% ----------- PM hit rate in high vs. low emphasis for focal condition only
PM_hit_low_focal = PM_hit(data(:, 3) == 0 & data(:, 2) == 1);
PM_hit_high_focal = PM_hit(data(:, 3) == 1 & data(:, 2) == 1);
[p, table] = anova1([PM_hit_high_focal PM_hit_low_focal], {'High, Focal', 'Low, Focal'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 1]; % it's F < 1 in paper
%%% ----------- PM hit rate in high vs. low emphasis for nonfocal condition only
PM_hit_low_nonfocal = PM_hit(data(:, 3) == 0 & data(:, 2) == 0);
PM_hit_high_nonfocal = PM_hit(data(:, 3) == 1 & data(:, 2) == 0);
[p, table] = anova1([PM_hit_high_nonfocal PM_hit_low_nonfocal], {'High, nonfocal', 'Low, nonfocal'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 5.90];
% ------------------ OG accuracy
[p, table] = anovan(OG_hit, {data(:, 1) data(:, 2) data(:, 3)}, 'model','full', 'display', 'off');
simulation_Fs = [simulation_Fs cell2mat(table(2:8, 6))'];
empirical_Fs = [empirical_Fs ones(1, 7)]; % it's F < 1 in paper and F = 1.48 for the last one
% ----------------- OG RT: focal vs. nonfocal
[p, table] = anovan(OG_RTs, {data(:, 2)}, 'model','full', 'display', 'off');
simulation_Fs = [simulation_Fs table{2, 6}];
empirical_Fs = [empirical_Fs 22.87];
% ----------------- OG RT: high vs. low emphasis
[p, table] = anovan(OG_RTs, {data(:, 3)}, 'model','full', 'display', 'off');
simulation_Fs = [simulation_Fs table{2, 6}];
empirical_Fs = [empirical_Fs 6.47];
% ----------------- OG RT: PM vs. No PM
[p, table] = anovan(OG_RTs, {data(:, 1)}, 'model','full', 'display', 'off');
simulation_Fs = [simulation_Fs table{2, 6}];
empirical_Fs = [empirical_Fs 131.665];
% ----------------- OG RT: cost, focal vs nonfocal
OG_RT_cost_focal = OG_RTs(data(:, 1) == 0 & data(:, 2) == 1) - OG_RTs(data(:, 1) == 1 & data(:, 2) == 1);
OG_RT_cost_nonfocal = OG_RTs(data(:, 1) == 0 & data(:, 2) == 0) - OG_RTs(data(:, 1) == 1 & data(:, 2) == 0);
[p, table] = anova1([OG_RT_cost_focal OG_RT_cost_nonfocal], {'Focal cost', 'Nonfocal cost'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 59.01]; % TODO maybe also include the M_cost ?
% ----------------- OG RT: cost, high emphasis vs low emphasis
OG_RT_cost_high = OG_RTs(data(:, 1) == 0 & data(:, 3) == 1) - OG_RTs(data(:, 1) == 1 & data(:, 3) == 1);
OG_RT_cost_low = OG_RTs(data(:, 1) == 0 & data(:, 3) == 0) - OG_RTs(data(:, 1) == 1 & data(:, 3) == 0);
[p, table] = anova1([OG_RT_cost_high OG_RT_cost_low], {'High emphasis cost', 'Low emphasis cost'}, 'off');
simulation_Fs = [simulation_Fs table{2, 5}];
empirical_Fs = [empirical_Fs 5.37]; % TODO maybe also include the M_cost ?
% ----------------- OG RT: cost, interaction focality & emphasis (not
% really calculating cost, but since baseline is the same, it's equivalent)
tEmp_Fs = [
    61.52, 127.96; % nonfocal: low, high
    1.73, 6.15     % focal: low, high
    ];
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        samples = data(data(:, 2) == FOCAL & data(:, 3) == EMPHASIS, :);
        [p, table] = anovan(samples(:,4), {samples(:, 1)}, 'model','full', 'display', 'off');
        simulation_Fs = [simulation_Fs table{2, 6}];
        empirical_Fs = [empirical_Fs tEmp_Fs(FOCAL + 1, EMPHASIS + 1)];
    end
end
fprintf('Fs\n');
simulation_Fs
empirical_Fs
F_errors = abs(simulation_Fs - empirical_Fs) ./ empirical_Fs .* 100; % TODO maybe err_scalers here too?
F_errors(isnan(F_errors)) = 10000; % NaNs here are REALLY bad 
F_error = sum(F_errors.^2) / length(F_errors);
fprintf('F error =                            %.4f\n', F_error);

%
% ------------------ compute the deviations and the final error
%

% -------------- Ball park for RT variances too, though I don't know if we really want to fit them
% Momchil: this is wrong; some of the SD's are RT's, the others are hit
% rates => can't regress
% ...removed
% TODO could scale the RT SD's by RT_slope

fprintf('M cols\n');
simulation_stats(:, M_cols)
empirical_stats(:,M_cols)


fprintf('SD cols\n');
simulation_stats(:,SD_cols)
empirical_stats(:,SD_cols)

deviations = (abs(simulation_stats(:,use_M_cols) - empirical_stats(:,use_M_cols))) ./ empirical_stats(:,use_M_cols) .* 100;
legit_deviations_count = sum(sum(~isnan(deviations)));
deviations(isnan(deviations)) = 0; % ignore NaN's e.g. in PM hit rates for OG_ONLY condition TODO fix and penalize real NaN's
deviations_error = sum(sum((deviations.^2))) / legit_deviations_count; % the squared error (Ms and SDs) scaled on each dim by err_scalers, averaged over all data points
fprintf('M deviations error =                 %.4f\n', deviations_error);

deviations = (abs(simulation_stats(:,use_SD_cols) - empirical_stats(:,use_SD_cols))) ./ empirical_stats(:,use_SD_cols) .* 100;
legit_deviations_count = sum(sum(~isnan(deviations)));
deviations(isnan(deviations)) = 0; % ignore NaN's e.g. in PM hit rates for OG_ONLY condition TODO fix and penalize real NaN's
deviations_error = sum(sum((deviations.^2))) / legit_deviations_count; % the squared error (Ms and SDs) scaled on each dim by err_scalers, averaged over all data points
fprintf('SD deviations error =               %.4f\n', deviations_error);

% now we can compute errors. let's do percent deviation from humans,
% squared (TODO or not squared?)
%
deviations = (abs(simulation_stats(:,use_cols) - empirical_stats(:,use_cols))) ./ empirical_stats(:,use_cols) .* 100;
legit_deviations_count = sum(sum(~isnan(deviations)));
deviations(isnan(deviations)) = 0; % ignore NaN's e.g. in PM hit rates for OG_ONLY condition TODO fix and penalize real NaN's
deviations_error = sum(sum((deviations.^2) .* err_scalers)) / legit_deviations_count; % the squared error (Ms and SDs) scaled on each dim by err_scalers, averaged over all data points
fprintf('M + SD deviations error (scaled) =  %.4f\n', deviations_error);

% this is an extra error term for differences between RT's in conditions
% -- makes figure 9 look prettier; kind of cheating but oh well
%
%diff_deviations = (empirical_stats([4 8],4) - empirical_stats([2 6],4)) - (simulation_stats([4 8],4) - simulation_stats([2 6],4));
%diff_deviations = diff_deviations ./ empirical_stats([4 8], 4) .* 100;

error = deviations_error + F_error / 1000; % + sum((diff_deviations.^2) .* diff_err_scalers);
fprintf('total error =                      %.4f\n', error);
