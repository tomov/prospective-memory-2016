function error = fit_exp1_and_exp2( free_params )
% the error function to fit with fmincon or whatever which takes only a set of "free" parameters
% and returns the computed error 

fprintf('\n\n\n     ======================>>>>> fitting [');
fprintf('%9.5f', free_params);
fprintf(']......\n');

debug_mode = false;
fitting_mode = false;
runs = 1; % how many times to run the experiment for each set of parameters; cost f'n is averaged

startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
            1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
            1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
            1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
            4 4 4, ... % biases -- tasks, attention, context
            0.3 0.3, ... % cross-subject init wm noise sigma -- PM task, target
            0.0004,  ... % gamma
            0.1 0.01, ...   % ffwd noise, wm noise sigma
            0.3];       % wm bias noise sigma

startpar([2 4 6 8 10 12 14 16 22]) = free_params; % set the params we're fitting
startpar(22) = free_params(9) * 10^(-3);

% experiment 1
[data, extra] = EM2005(startpar, 1, fitting_mode, debug_mode, false, runs);
errors_exp1 = zeros(runs, 1);
for run = 1:runs
    errors_exp1(run) = compute_err_exp1(data{run}, extra);
end
error_exp1 = mean(errors_exp1);
error_exp1_std = std(errors_exp1);
fprintf('\n                   .............. AVERAGE ERROR = %.4f (std = %.4f)\n', error_exp1, error_exp1_std);

% experiment 2
% TODO UNDO
%[data, extra] = EM2005(startpar, 2, fitting_mode, debug_mode, false);
error_exp2 = 0; %error_exp2 = compute_err_exp2(data, extra);

error = error_exp1 + error_exp2;
%fprintf('\n\n\n     ======================>>>>> fitting [');
%fprintf('%8.2f', free_params);
%fprintf(', error = %f * 3 + %f = %f\n\n\n', error_exp1, error_exp2, error);

end

