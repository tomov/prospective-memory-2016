function error = compute_err_exp2(data, extra) % takes in the output of EM2005
% Computes the error for a particular simulation of EM2005

%{
 subjects[subj_id, :] ==
 samples variable order = 
    1 - OG_ONLY,
    2 - FOCAL, 
    3 - EMPHASIS
    10 - block id
 statistics order = 
    4 - OG_RT, 
    5 - OG_Hit, 
    6 - PM_RT, 
    7 - PM_Hit,
    8 - PM_miss_OG_hit
    9 - subject id
 (see EM2005 exp 2)
%}

DO_PLOT = false;
blocks = data;

% -------------- define the empirical stats (Table 1 from E&M 2005)

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

%%% CONFIG %%%
SD_cols = [5,7,9,11]; % SDs. we convert those to SEM by dividing by subjects_per_condition
RT_mean_cols = [4,8]; % RTs. we regress those against model cycles means to scale
RT_SD_cols = [5,9]; % RT SDs. we regress those against model cycles SDs to scale

subjects_per_condition = 24;

% comes from E&M 2005 paper
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
%
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);

% which columns to use in the error calculation
%
use_cols = [4  6, 10]; % OG RT's, OG accuracy, and PM hit rates

% resize weights of errors to # of conditions
%
err_weights = [1, 1.5, 1.5];  % weights of the errors -- keep in mind we're scaling the deviations to %'s so they're comparable. Give more weights to the PM hit rates
err_scalers = repmat(err_weights, size(empirical_stats, 1), 1);
err_scalers(empirical_stats(:, 1) == 0 & empirical_stats(:, 2) == 0, 1) = 2; % enhance the OG RT's in the PM, nonfocal conditions for prettyness (to balance out the other conditions which are effectively the same thing but multiple times)
err_scalers(empirical_stats(:, 2) == 1, 1) = 0.1; % the OG RT's focal conditions b/c they are bogus (should be same as nonfocal, no-PM condition)


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
                M = mean(samples);
                SD = std(samples);
                %assert(length(samples) == subjects_per_condition);
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

%TODO FIXME use slope and intercept from experiment 1 ?
%RT_slope = 10;
%RT_intercept = 100;

% replace model cycles with fit RT values
simulation_stats(:,RT_mean_cols) = polyval(p, simulation_stats(:,RT_mean_cols));

% now we can compute errors. let's do percent deviation from humans,
% squared (TODO or not squared?)
%
deviations = (abs(simulation_stats(:,use_cols) - empirical_stats(:,use_cols))) ./ empirical_stats(:,use_cols) .* 100;
deviations(isnan(deviations)) = 0; % ignore NaN's

error = sum(sum((deviations.^2) .* err_scalers));
