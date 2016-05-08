function error = fitme( free_params )
% the error function to fit with fmincon or whatever which takes only a set of "free" parameters
% and returns the computed error 

fprintf('\n\n\n     ======================>>>>> fitting [');
fprintf('%8.5f', free_params);
fprintf(']......\n');

debug_mode = false;
experiment = 1;

startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
            1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
            1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
            1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
            4 4 4, ... % biases
            0 0];      % no noise...

startpar([2 4 6 8 10 12 14 16]) = free_params; % set the params we're fitting

[data, extra] = EM2005(startpar, experiment, debug_mode, false);

error = compute_err(data, extra);

%fprintf('\n\n\n     ======================>>>>> fitting [');
%fprintf('%8.2f', free_params);
fprintf(', error = %f\n\n\n', error);

end

