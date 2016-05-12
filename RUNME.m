warning('off', 'MATLAB:ClassInstanceExists');
clear classes % ! super important ! if you don't do this, MATLAB won't reload your classes

% best parameters so far...
% KEEP I_WM equal for both task and feature units

%{
params = 
  [WM units focal, low emph ...
   WM units focal, high emph ...
   WM units nonfocal, low emph ...
   WM units nonfocal, high emph ...
   WM bias,
   WM bias];
where
  WM units = [
    OG Task     PM Task     OG features     Monitor tortoise    Monitor tor
  ];
%}

debug_mode = false; % only run 1 subject per condition and show progress ; must change parfor to for in EM2005
fitting_mode = false; % used when fitting the parameters ; uses a more efficient setup that produces similar results
experiment = 3;

if experiment == 1
    % all 4 conditions; RT's and hit rates.
    % free params = startpar
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ... % biases
		        0 0, ...   % no noise...
                0.0004];   % gamma
            
            
    startpar([2 4 6 8 10 12 14 16 22])  =  [ 0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180 / 1000];
         % [0.2453    0.3713    0.5446    0.4261    0.8402    0.7132
         % 0.7766    0.8729]; % best par from solve.m for experiment 1
elseif experiment == 2
    % decay of monitoring => PM hit rate lowers over time in nonfocal
    % free param -- gamma = 0.0004
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.78, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ... % biases
		        0 0, ...   % no noise...
                0.0004];   % gamma
    startpar([2 4 6 8 10 12 14 16 22])  =  [ 0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180 / 1000];
            
elseif experiment == 3 
    % 6 targets => slower OG in 6 vs. 1 target
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.0    1    0.6, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.0    1    0.7, ...     % focal, high emph      % exp1_v16
                1  0.3    1    0.5, ...    % INVALID? nonfocal, low emph   % exp2_v11
                1  0.6    1    0.5, ...    % INVALID? nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...   % biases
         		0 0, ...     % no noise
                0.0004];
elseif experiment == 4  
    % cross-subject monitoring strategies => half of subjects slower OG than
    % other half
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.1    1    0.1, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % INVALID? focal, high emph      % exp1_v16
                1  0.8    1    0.75, ...    % INVALID? nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % INVALID? nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...   % biases
         		1 0.8, ...   % uniform noise
                0.0004];     % gamma
end

tic
[data, extra] = EM2005(startpar, experiment, fitting_mode, debug_mode, true);
toc

if debug_mode
	m = Model(true, startpar([1 2 3 4 17 18 19 20 21]), false);
    wm_ids = m.wm_ids;
    context_ids = m.context_ids;
    act = extra{1, 8};
    nets = extra{1, 12};
   % figure;
   % plot([act(1:100, context_ids), nets(1:100, context_ids)]);
else
    save('exp4-data-newww.mat');
end


        
%{
[data, ~] = EM2005(startpar, 1);
data
save('rondo-run-data-exp-1.mat');
EM2005_with_stats_exp1
save('rondo-run-data-exp-1-with-stats.mat');


data_exp1 = data;


[data, ~] = EM2005(startpar, 2);
data
save('rondo-run-data-exp-12.mat');
EM2005_with_stats_exp2
save('rondo-run-data-exp-12-with-stats.mat');


data_exp2 = data;


[data, ~] = EM2005(startpar, 3);
data
save('rondo-run-data-exp-123.mat');
%EM2005_with_stats_exp1
%save('rondo-run-data-exp-123-with-stats.mat');


data_exp3 = data;

save('goodmorning.mat');
            
            
%}
            

%% ----- BIG TODOs ----

%{
make experiment 5 gather relevant data



make number of subjects normal

go and make sure there are no hacks left (FIXME) and no hardcoded numbers
  
send them all as jobs 
%}

