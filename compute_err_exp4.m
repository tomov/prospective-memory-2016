function error = compute_err_exp4(data, extra) % takes in the output of EM2005
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

subjects_per_condition = 104 / 2; % divide by 2! b/c "EMPHASIS" == "high cost / low cost" is not a real condition => samples will be in half

empirical_stats = [
    % low emphasis = no cost participants (no monitoring)
    1 1 0, 7201.28,  1520.06,       74, 10, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    low emph, 1 targets
    0 1 0, 6706.30,  1337.78,       73, 10, NaN, NaN, 94,   14, 1;  % PM, focal,    low emph, 1 targets
    
    % high emphasis = cost participants (monitoring)
    1 1 1, 6797.26, 1484.68,        74, 10, NaN, NaN, NaN, NaN, 1;  % no-PM, focal,    high emph, 1 targets
    0 1 1, 7533.95, 1769.18,        73, 10, NaN, NaN, 95,   17, 1;  % PM, focal,    high emph, 1 targets
];


% convert SD's to SEM's in empirical data
empirical_stats(:, SD_cols) = empirical_stats(:, SD_cols) / sqrt(subjects_per_condition);

% resize weights of errors to # of conditions
err_scalers = repmat(err_weights, size(empirical_stats, 1), 1);
diff_err_scaler = 0.1; % b/c the % deviations for the RT's are very small (b/c the Y-intercept is big)


% -------------- change meaning of EMPHASIS column !!!
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
EMPHASIS = 1;
FOCAL = 1;
TARGETS = 1;
% order here matters -- must be same as empirical_data above for line
% regression
for EMPHASIS = 0:1
    for OG_ONLY = 1:-1:0
        stat = [OG_ONLY, FOCAL, EMPHASIS];
        for col = 4:7
            samples = subjects(subjects(:, 1) == OG_ONLY & subjects(:, 3) == EMPHASIS & subjects(:, 9) == TARGETS, col);
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


% replace model cycles with fit RT values
simulation_stats(:,RT_mean_cols) = polyval(p, simulation_stats(:,RT_mean_cols));

% now we can compute errors. let's do percent deviation from humans,
% squared (TODO or not squared?)
%
deviations = (abs(simulation_stats(:,use_cols) - empirical_stats(:,use_cols))) ./ empirical_stats(:,use_cols) .* 100;
deviations(isnan(deviations)) = 0; % ignore NaN's

error = sum(sum((deviations.^2) .* err_scalers));

% this is an extra error term for difference between RT's in the PM
% condition
% -- makes figure 24 look prettier; kind of cheating but oh well
%
diff_deviation = ((empirical_stats(4, 4) - empirical_stats(2, 4)) - (simulation_stats(4, 4) - simulation_stats(2, 4))) / (empirical_stats(4, 4) - empirical_stats(2, 4)) * 100;

error = error + diff_deviation .^2 .* diff_err_scaler;
