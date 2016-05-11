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
SD_cols = [5,7,9,11]; % SDs. we convert those to SEM by dividing by subjects_per_condition
RT_mean_cols = [4,8]; % RTs. we regress those against model cycles means to scale
RT_SD_cols = [5,9]; % RT SDs. we regress those against model cycles SDs to scale

use_cols    = [4  6, 10]; % which columns to use in the error calculation -- OG RT's, OG accuracy, and PM hit rates
err_weights = [1, 1, 1];  % weights of the errors -- keep in mind we're scaling them to %'s so they're all on the same scale

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


% ------------- calculate simulation stats (Table 1 from E&M 2005)

simulation_stats = [];
for FOCAL = 1:-1:0
    for EMPHASIS = 0:1
        for OG_ONLY = 1:-1:0            
            stat = [OG_ONLY, FOCAL, EMPHASIS];
            for col = 4:7
                samples = data(data(:, 1) == OG_ONLY & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS, col);
                samples = samples(~isnan(samples));
                M = mean(samples, 'omitnan');
                SD = std(samples, 'omitnan');
                %assert(length(samples) == subjects_per_condition);
                stat = [stat, M, SD];
            end
            simulation_stats = [simulation_stats; stat];
        end
    end
end

% convert SD's to SEM's in simulation data
simulation_stats(:, SD_cols) = simulation_stats(:, SD_cols) / sqrt(subjects_per_condition);



% -------------- run linear regression to find slope and intercept for RTs

empirical_RTs = empirical_stats(:, 4);
simulation_cycles = simulation_stats(:, 4);

p = polyfit(simulation_cycles, empirical_RTs, 1);
RT_slope = p(1);
RT_intercept = p(2);
yfit = polyval(p, simulation_cycles);

% replace model cycles with fit RT values
simulation_stats(:,RT_mean_cols) = polyval(p, simulation_stats(:,RT_mean_cols));

% -------------- Ball park for RT variances too, though I don't know if we really want to fit them
% Momchil: this is wrong; some of the SD's are RT's, the others are hit
% rates => can't regress
% ...removed
% TODO could scale the RT SD's by RT_slope

% now we can compute errors. let's do percent deviation from humans,
% squared (TODO or not squared?)
%
deviations = (abs(simulation_stats(:,use_cols) - empirical_stats(:,use_cols))) ./ empirical_stats(:,use_cols) .* 100;
deviations(isnan(deviations)) = 0; % ignore NaN's

% this is an extra error term for differences between RT's in conditions
% -- makes figure 9 look prettier; kind of cheating but oh well
%
%diff_deviations = (empirical_stats([4 8],4) - empirical_stats([2 6],4)) - (simulation_stats([4 8],4) - simulation_stats([2 6],4));
%diff_deviations = diff_deviations ./ empirical_stats([4 8], 4) .* 100;

error = sum(sum((deviations.^2) .* err_scalers)); % + sum((diff_deviations.^2) .* diff_err_scalers);
