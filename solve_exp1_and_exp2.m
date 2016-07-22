% Find the best set of free parameters
%

         % PM Task  PM target(s) initial WM activations
init_par = [0.35   0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
            0.6    0.4, ...     % focal, high emph      % exp1_v16
            0.8    0.75, ...    % nonfocal, low emph   % exp2_v11
            0.9    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
            0.4];   % gamma * 10^-3
        
        
min_par =  [0 0 0 0 0 0 0 0 0];
max_par =  [1 1 1 1 1 1 1 1 1];
assert(length(min_par) == length(init_par));

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter', 1000, 'DiffMinChange', 0.001);
best_par = fmincon(@fit_exp1_and_exp2, init_par, [], [], [], [], min_par, max_par, [], options);

best_par

% Run experiment 1 and save the actual data
%
experiment = 1;

rng('shuffle');
runhash = randstr(10); % so we can tag the output files

[data, extra] = EM2005(best_par, experiment, false, false, true, 1); 

filename = sprintf('exp%d-runhash-%s.mat', experiment, runhash);
fprintf('\n\n\nOUTPUT SAVED TO FILE %s\n', filename);
save(filename);



% SOLUTIONS FOR EXPERIMENT 1 ---
%

%  fitting [ 0.27808 0.42082 0.48387 0.48294 0.84110 0.72067 0.88019
%  0.86090] --> 740
%  0.28977 0.41855 0.50422 0.50325 0.83942 0.70892 0.91721 0.85925]. -->
%  730

%   [ 0.24536 0.37131 0.54459 0.42612 0.84025 0.71321 0.77664 0.87286]. -->
%   586

% itting [ 0.24536 0.37131 0.54459 0.42612 0.84025 0.71821 0.77664
% 0.87286]. --> 673

% [ 0.23122 0.36407 0.57085 0.40156 0.80379 0.71066 0.73187 0.88019]  -->
% 660

%  [ 0.24420 0.37071 0.54675 0.42410 0.83725 0.71300 0.77295 0.87346] -->
%  628

% best par = 0.2453    0.3713    0.5446    0.4261    0.8402    0.7132    0.7766    0.8729


% EXPERIMENT 2

%    ======================>>>>> fitting [  0.36285  0.31384  0.60791  0.39209  0.80395  0.73517  0.88220  0.83836  0.41186]......
% , error = 767.963778 * 3 + 3652.177395 = 5956.068727
%
% best par =  0.3728    0.2708    0.5854    0.4157    0.7403    0.7339    0.8124    0.8061    0.4180
