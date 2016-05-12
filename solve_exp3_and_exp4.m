% Find the best set of free parameters
% TODO only works for experiment 1
%

         % PM Task  PM target(s) initial WM activations
init_par = [0       0.6, ...     % focal, low emph 
            0       0, ...       % noise
            0.4];         % gamma * 10^-3
min_par =  [0 0 0 0 0];
max_par =  [1 1 1 1 1];

options = optimoptions(@fmincon,'Algorithm','sqp','MaxIter', 1000, 'DiffMinChange', 0.001);
best_par = fmincon(@fit_exp3_and_exp4, init_par, [], [], [], [], min_par, max_par, [], options);

best_par
