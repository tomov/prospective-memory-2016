% Find the best set of free parameters
% TODO only works for experiment 1
%


init_par = [0.35    0.3, ...     % focal, low emph 
            0.6     0.4, ...     % focal, high emph 
            0.8     0.75, ...    % nonfocal, low emph 
            0.9     0.83];
min_par =  [0 0 0 0 0 0 0 0];
max_par =  [1 1 1 1 1 1 1 1];

best_par = fmincon(@fitme, init_par, [], [], [], [], min_par, max_par);

best_par