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


    % experiment 1
    startpar = [1  0.35   1    0.3, ...      % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...    % focal, high emph      % exp1_v16
                1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ... % biases
		0 0];      % noise -- no noise...
 
debug_mode = true;
backprop = false;
experiment = 1;

tic
[data, extra] = EM2005(startpar, experiment, debug_mode, backprop);
toc

%save('exp1-data-newww.mat');


%{
data

if debug_mode
	m = Model(startpar, false);
    wm_ids = m.wm_ids;
    context_ids = m.context_ids;
    act = extra{1, 8};
    nets = extra{1, 12};
   % figure;
   % plot([act(1:100, context_ids), nets(1:100, context_ids)]);
else
    filename = sprintf('/mnt/cd/people/mtomov/data/exp1-data.mat');
    save('exp1-data-new.mat');
    if experiment == 1
        EM2005_with_stats_exp1
    elseif experiment == 2
        EM2005_with_stats_exp2
    elseif experiment == 3
        EM2005_with_stats_exp3
    end
end
%}

        
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

