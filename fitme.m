function error = fitme( free_params )
% the error function to fit with fmincon or whatever which takes only a set of "free" parameters
% and returns the computed error 

fprintf('\n\n\n     ======================>>>>> fitting [');
fprintf('%9.5f', free_params);
fprintf(']......\n');

debug_mode = false;
fitting_mode = true;

startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
            1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
            1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
            1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
            4 4 4, ... % biases
            0 0, ...
            0.0004];      % no noise...

startpar([2 4 6 8 10 12 14 16 22]) = free_params; % set the params we're fitting
startpar(22) = free_params(9) * 10^(-3);

% experiment 1
[data, extra] = EM2005(startpar, 1, fitting_mode, debug_mode, false);
error_exp1 = compute_err_exp1(data, extra);

% experiment 2
[data, extra] = EM2005(startpar, 2, fitting_mode, debug_mode, false);
error_exp2 = compute_err_exp2(data, extra);

error = error_exp1 * 3 + error_exp2;
%fprintf('\n\n\n     ======================>>>>> fitting [');
%fprintf('%8.2f', free_params);
fprintf(', error = %f * 3 + %f = %f\n\n\n', error_exp1, error_exp2, error);

end

