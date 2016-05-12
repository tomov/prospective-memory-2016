function error = compute_err_exp3(data, extra) % takes in the output of EM2005
% Computes the error for a particular simulation of EM2005

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
%%% CONFIG %%%
SD_cols = [5,7,9,11]; % SDs. we convert those to SEM by dividing by subjects_per_condition
RT_mean_cols = [4,8]; % RTs. we regress those against model cycles means to scale
RT_SD_cols = [5,9]; % RT SDs. we regress those against model cycles SDs to scale

use_cols    = [4, 6, 10]; % which columns to use in the error calculation -- OG RT's, OG accuracy, and PM hit rates
% weights of the errors -- keep in mind we're scaling them to %'s so they're all on the same scale
%
err_weights = [1, 0.0, 0.0];  % the OG RT's are most important; also they're big => the % deviations are small compared to the hit rate % deviations

subjects_per_condition = 32;

empirical_stats = [
    1 1 0, 4791, 618,        70, 11, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    low emph, 1 targets
    1 1 0, 4890, 508,        70, 11, NaN, NaN, NaN, NaN, 6;  % no-PM, focal,    low emph, 6 targets
    0 1 0, 4885, 591,        69, 10, NaN, NaN, 80,   28, 1;  % PM, focal,    low emph, 1 targets
    0 1 0, 5215, 422,        69, 10, NaN, NaN, 72,   25, 6;  % PM, focal,    low emph, 6 targets
];


% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);

% resize weights of errors to # of conditions
err_scalers = repmat(err_weights, size(empirical_stats, 1), 1);
diff_err_scaler = 0.1; % b/c the % deviations for the RT's are very small (b/c the Y-intercept is big)


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

% replace model cycles with fit RT values
simulation_stats(:,RT_mean_cols) = polyval(p, simulation_stats(:,RT_mean_cols));

% now we can compute errors. let's do percent deviation from humans,
% squared (TODO or not squared?)
%
deviations = (abs(simulation_stats(:,use_cols) - empirical_stats(:,use_cols))) ./ empirical_stats(:,use_cols) .* 100;
deviations(isnan(deviations)) = 0; % ignore NaN's

error = sum(sum((deviations.^2) .* err_scalers));

% the OG RT slowing for 6 targets (PM vs. non-PM) should be ~300 msec, add that to the error too
% TODO necessary?  Figure 20 is prettier
%
OG_6_target_slowing = simulation_stats(simulation_stats(:,1) == 0 & simulation_stats(:, 12) == 6, 4) ...
    -  simulation_stats(simulation_stats(:,1) == 1 & simulation_stats(:, 12) == 6, 4);
error = error + (((OG_6_target_slowing - 300) / 300) * 100) .^ 2 * diff_err_scaler;