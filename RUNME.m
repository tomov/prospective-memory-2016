warning('off', 'MATLAB:ClassInstanceExists');
clear classes % ! super important ! if you don't do this, MATLAB won't reload your classes

rng('shuffle');

runhash = randstr(10); % so we can tag the output files
fprintf('runhash = %s\n', runhash);

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

debug_mode = false; % only run 1 subject per condition and show progress ; !!!!!IMPORTANT!!!!!!!!! must change parfor to for in EM2005



experiment = 6; % <------------------------------- HERE --------------------------------------



if experiment == 1
    % all 4 conditions; RT's and hit rates.
    % free params = startpar
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph
                1  0.6    1    0.6, ...     % focal, high emph
                1  0.8    1    0.7, ...    % nonfocal, low emph
                1  0.9    1    0.9, ...    % nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
		        0.2 0.2, ... % cross-subject init wm noise sigma -- PM task, target
                0.0004,  ... % gamma
                0.1 0.01, ...% ffwd noise, wm noise sigma
                0.1, ...     % wm bias noise sigma
                0.0, ...     % OG weights noise
                NaN NaN NaN, ... % low WM bias => not applicable here
                NaN, ... % PM Context during the third task => not applicable here
                0 0]; % extra weight to add to EM connections => not applicable here (leave as zero)
            
elseif experiment == 2
    % decay of monitoring => PM hit rate lowers over time in nonfocal
    % free param -- gamma = 0.0004
    % same design as experiment 1
    %
    % OG task, PM task, OG features, target(s)    
    startpar = [1  0.35   1    0.3, ...     % focal, low emph
                1  0.6    1    0.6, ...     % focal, high emph
                1  0.8    1    0.7, ...    % nonfocal, low emph
                1  0.9    1    0.9, ...    % nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
		        0.2 0.2, ... % cross-subject init wm noise sigma -- PM task, target
                0.0004,  ... % gamma
                0.1 0.01, ...% ffwd noise, wm noise sigma
                0.1, ...     % wm bias noise sigma
                0.0, ...      % OG weights noise
                NaN NaN NaN, ... % low WM bias => not applicable here            
                NaN, ... % PM Context during the third task => not applicable here
                0 0]; % extra weight to add to EM connections => not applicable here (leave as zero)
            
elseif experiment == 3 
    % 6 targets => slower OG in 6 vs. 1 target
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.0    1    0.6, ...     %  0.0  0.6  focal, low emph
                NaN NaN   NaN  NaN, ...    % INVALID focal, high emph
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, low emph
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
         		0.2 0.2, ...     % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...  %xssssssssssssssssaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa  gamma
                0.1 0.01, ...   % ffwd noise, wm noise sigma
                0.1, ...    % wm bias noise sigma
                0.4, ...    % OG weights noise (e.g. 0.4)
                NaN NaN NaN, ... % low WM bias => not applicable here
                NaN, ... % PM Context during the third task => not applicable here
                0 0]; % extra weight to add to EM connections => not applicable here (leave as zero)
                        
elseif experiment == 4  
    % cross-subject monitoring strategies => half of subjects slower OG than
    % other half
    % same design as experiment 3
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.0    1    0.6, ...     % focal, low emph
                NaN NaN   NaN  NaN, ...    % INVALID focal, high emph
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, low emph
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
         		0.4 0.4, ... %1 0.8, ...   % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...  % gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.1, ...     % wm bias noise sigma
                0.4, ...      % OG weights noise
                NaN NaN NaN, ... % low WM bias => not applicable here
                NaN, ... % PM Context during the third task => not applicable here
                0 0]; % extra weight to add to EM connections => not applicable here (leave as zero)

elseif experiment == 5
    % interleave third task, see aftereffects of intention
    % free params = pm context turned off ?
    % same as experiment 1, low emphasis, focal
    %
    % OG task, PM task, OG features, target(s)
    
    startpar = [1  0.35   1    0.3, ...     % focal, low emph
                1  0.6    1    0.6, ...     % focal, high emph
                1  0.8    1    0.7, ...    % nonfocal, low emph
                1  0.9    1    0.9, ...    % nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
		        0.2 0.2, ... % cross-subject init wm noise sigma -- PM task, target
                0.0004,  ... % gamma
                0.1 0.01, ...% ffwd noise, wm noise sigma
                0.1, ...     % wm bias noise sigma
                0.0, ...     % OG weights noise
                NaN NaN NaN, ... % low WM bias => not applicable here
                0.5, ... % PM Context during the third task
                0 0]; % extra weight to add to EM connections => not applicable here (leave as zero)
            
elseif experiment == 6
    % Brewer et al. 2010 -- low vs. high wm capacity
    % free params = biases
    % same design as experiment 1
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph
                NaN  NaN  NaN  NaN, ...     % focal, high emph
                1  0.8    1    0.7, ...    % nonfocal, low emph
                NaN  NaN  NaN  NaN, ...    % nonfocal, high emph
                4.2  4.2  4.2, ...  % biases, high wm capacity -- tasks, attention, context
		        0.1 0.1, ...        % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...     % gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.1, ...    % wm bias noise sigma
                0.0, ...    % OG weights noise
                3.8 3.8 3.8, ...       % biases, low wm capacity
                NaN, ... % PM Context during the third task => not applicable here
                0 0]; % extra weight to add to EM connections => set to 50, 10 for EM strength effect
            
     % !!!!!IMPORTANT!!!!!!
     startpar([5 6 7 8]) = startpar([1 2 3 4]); % high emphasis = low emphasis (b/c we use it to mean "wm capacity" #hacksauce)
     startpar([13 14 15 16]) = startpar([9 10 11 12]); % high emphasis = low emphasis (b/c we use it to mean "wm capacity" #hacksauce)
elseif experiment == 7
    % explore aftereffects of intention with third task only (not interleaved)
    % very similar to experiment 5
    % free params = pm context turned off ?
    % same params as experiment 1, low emphasis, focal
    %
    % OG task, PM task, OG features, target(s)
    
    startpar = [1  0.35   1    0.3, ...     % focal, low emph
                1  0.6    1    0.6, ...     % focal, high emph
                1  0.8    1    0.7, ...    % nonfocal, low emph
                1  0.9    1    0.9, ...    % nonfocal, high emph
                4 4 4, ...   % biases -- tasks, attention, context
		        0.2 0.2, ... % cross-subject init wm noise sigma -- PM task, target
                0.0004,  ... % gamma
                0.1 0.01, ...% ffwd noise, wm noise sigma
                0.1, ...     % wm bias noise sigma
                0.0, ...     % OG weights noise
                NaN NaN NaN, ... % low WM bias => not applicable here
                0.5, ... % PM Context during the third task
                0 0]; % extra weight to add to EM connections => for strong EM effect, set to 15 0 with PM Context = 0.5
            
end

tic
[data, extra] = EM2005(startpar, experiment, debug_mode, true, 1, 150);
data = data{1};
extra = extra{1};
toc

if ~debug_mode
    filename = sprintf('exp%d-runhash-%s.mat', experiment, runhash);
    fprintf('\n\n\nOUTPUT SAVED TO FILE %s\n', filename);
    save(filename);
end
