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
fitting_mode = false; % used when fitting the parameters ; uses a more efficient setup that produces similar results




experiment = 1; % <------------------------------- HERE --------------------------------------





if experiment == 1
    % all 4 conditions; RT's and hit rates.
    % free params = startpar
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...   % biases -- tasks, attention, context
		        0.3 0.3, ... % cross-subject init wm noise sigma -- PM task, target
                0.0004,  ... % gamma
                0.1 0.01, ...% ffwd noise, wm noise sigma
                0.0];        % wm bias noise sigma
            
  %   startpar([2 4 6 8 10 12 14 16 22])  =  [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131 / 1000];
  %   startpar([2 4 6 8 10 12 14 16 20 21 22])  =  [0.3027    0.2594    0.5188    0.3459    0.6918    0.6486    0.7783 0.7177    0.2121    0.1297    0.3459 / 1000] % M * 2 + SD + F / 1000
  %   startpar([2 4 6 8 10 12 14 16 22]) = [0.3744    0.2865    0.5731 ...
  %   0.4222    0.7641    0.7163    0.8997 0.7927    0.4222 / 1000];  % fit
  %   Ms <--- BAD
  %   startpar([2 4 6 8 10 12 14 16 22]) = [0.3428    0.2938    0.5884 ...
  %   0.3918    0.7843    0.7346    0.8815 0.8129    0.3926 / 1000]; % fit
  %   Ms  <-- BAD
  %   startpar([2 4 6 8 10 12 14 16 22]) = [0.3678    0.3193    0.5822 ...
  %   0.3881    0.8047    0.7285    0.8732 0.8053    0.3881 / 1000]; % fit
  %   Ms  <-- BAD
     
   % startpar([2 4 6 8 10 12 14 16 22])  =  [ 0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180 / 1000];
         % [0.2453    0.3713    0.5446    0.4261    0.8402    0.7132
         % 0.7766    0.8729]; % best par from solve.m for experiment 1
        % from git hash a460e3a5c78a0492811b42a831f37c8b7023e364  ---> [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131]
elseif experiment == 2
    % decay of monitoring => PM hit rate lowers over time in nonfocal
    % free param -- gamma = 0.0004
    % same design as experiment 1
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.78, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...    % biases -- tasks, attention, context
		        0.2 0.2, ...  % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...   % gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.0];         % wm bias noise sigma

     startpar([2 4 6 8 10 12 14 16 22])  =  [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131 / 1000];
%    startpar([2 4 6 8 10 12 14 16 22])  =  [ 0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180 / 1000];
 % from git hash a460e3a5c78a0492811b42a831f37c8b7023e364  ---> [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131]
elseif experiment == 6
    % Brewer et al. 2010 -- low vs. high wm capacity
    % free params = biases
    % same design as experiment 1
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.70, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4.5   4.5   4.5, ...  % biases, high wm capacity -- tasks, attention, context
		        0 0, ...        % cross-subject init wm noise sigma -- PM task, target
                0.0002, ...     % gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.0, ...    % wm bias noise sigma
                4 4 4];       % biases, low wm capacity
            
     % biases = 3.5 3.3 4 not bad result ; nvm bad result -> features must
     % be high.... plots 1 and 2 ok
     % biases = 3.5 4 4 => plots 2 and 4 ok
     % 
     startpar([2 4 6 8 10 12 14 16 22])  =  [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131 / 1000]; % from experiment 1
     
     startpar([5 6 7 8]) = startpar([1 2 3 4]); % high emphasis = low emphasis (b/c we use it to mean "wm capacity" #hacksauce)
     startpar([13 14 15 16]) = startpar([9 10 11 12]); % high emphasis = low emphasis (b/c we use it to mean "wm capacity" #hacksauce)
 %   startpar([2 4 6 8 10 12 14 16 22])  =  [ 0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180 / 1000];
         % [0.2453    0.3713    0.5446    0.4261    0.8402    0.7132
         % 0.7766    0.8729]; % best par from solve.m for experiment 1
    
elseif experiment == 3 
    % 6 targets => slower OG in 6 vs. 1 target
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.0    1    0.6, ...     %  0.0  0.6  focal, low emph     % exp1_v16, exp2_v19
                NaN NaN   NaN  NaN, ...    % INVALID focal, high emph      % exp1_v16
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, low emph   % exp2_v11
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...   % biases -- tasks, attention, context
         		0 0, ...     % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...  % gamma
                0.1 0.01, ...   % ffwd noise, wm noise sigma
                0.0];     % wm bias noise sigma
            
    startpar([2 4 20 21]) = [0    0.1368    0.7049    0.5319];
     % startpar([2 4 20 21 22]) = [ 0.0155         0    1.0000         0         0 / 1000 ]; % probabilistic -- sometimes work sometimes not so well (b/c of the randomness in monitoring)
% from git hash a460e3a5c78a0492811b42a831f37c8b7023e364  --->     [0.2428    1.0000    1.0000    1.0000]
elseif experiment == 4  
    % cross-subject monitoring strategies => half of subjects slower OG than
    % other half
    % same design as experiment 3
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.1    1    0.1, ...     % focal, low emph     % exp1_v16, exp2_v19
                NaN NaN   NaN  NaN, ...    % INVALID focal, high emph      % exp1_v16
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, low emph   % exp2_v11
                NaN NaN   NaN  NaN, ...    % INVALID nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ...   % biases -- tasks, attention, context
         		1 0.8, ...   % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...  % gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.0];     % wm bias noise sigma

    startpar([2 4 20 21]) = [0    0.1368    0.7049    0.5319];            
%    startpar([2 4 20 21 22]) = [ 0.0155         0    1.0000         0         0 / 1000 ];
% from git hash a460e3a5c78a0492811b42a831f37c8b7023e364  --->     [0.2428
% 1.0000    1.0000    1.0000] ..... also wtf [0    0.0818         0
% 0.0053]
elseif experiment == 5
    % interleave third task, see aftereffects of intention
    % free params = pm context turned off ?
    % same as experiment 1, low emphasis, focal
    %
    % OG task, PM task, OG features, target(s)
    startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
                1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
                1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
                1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
                4 4 4, ... % biases -- tasks, attention, context
		        0 0, ...   % cross-subject init wm noise sigma -- PM task, target
                0.0004, ...% gamma
                0.1 0.01, ... % ffwd noise, wm noise sigma
                0.0];     % wm bias noise sigma
            
            
     startpar([2 4 6 8 10 12 14 16 22])  =  [0.3425    0.2937    0.5865    0.4131    0.7830  0.7332    0.9019 0.8114    0.4131 / 1000]; % from experiment 1

end

tic
[data, extra] = EM2005(startpar, experiment, fitting_mode, debug_mode, true, 1);
data = data{1};
toc

if ~debug_mode
    filename = sprintf('exp%d-runhash-%s.mat', experiment, runhash);
    fprintf('\n\n\nOUTPUT SAVED TO FILE %s\n', filename);
    save(filename);
end
