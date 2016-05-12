function error = fit_exp3_and_exp4( free_params )
% the error function to fit with fmincon or whatever which takes only a set of "free" parameters
% and returns the computed error 

fprintf('\n\n\n     ======================>>>>> fitting [');
fprintf('%5.5f', free_params);
fprintf(']......\n');

debug_mode = false;
fitting_mode = true;

startpar = [1  0.0    1    0.6, ...     %  0.0  0.6  focal, low emph     % exp1_v16, exp2_v19
            NaN NaN   NaN  NaN, ...    % INVALID focal, high emph      % exp1_v16
            NaN NaN   NaN  NaN, ...    % INVALID nonfocal, low emph   % exp2_v11
            NaN NaN   NaN  NaN, ...    % INVALID nonfocal, high emph  % exp1_v16 -- sorta
            4 4 4, ...   % biases
            0 0, ...     % no noise
            0.0004];

startpar([2 4 20 21 22]) = free_params; % set the params we're fitting
startpar(22) = free_params(5) * 10^(-3);

% experiment 3
[data, extra] = EM2005(startpar, 3, fitting_mode, debug_mode, false);
error_exp3 = compute_err_exp3(data, extra);

% experiment 4
[data, extra] = EM2005(startpar, 4, fitting_mode, debug_mode, false);
error_exp4 = compute_err_exp4(data, extra);

error = error_exp3 * 2 + error_exp4;
%fprintf('\n\n\n     ======================>>>>> fitting [');
%fprintf('%5.2f', free_params);
fprintf(', error = %f * 2 + %f = %f\n\n\n', error_exp3, error_exp4, error);

end

