function error = fit_exp1_and_exp2( free_params, runs )
% the error function to fit with fmincon or whatever which takes only a set of "free" parameters
% and returns the computed error 

fprintf('\n\n\n     ======================>>>>> fitting [');
fprintf('%9.5f', free_params);
fprintf(']......\n');

debug_mode = false;
parallel_runs = false;

startpar = free_params_to_params(free_params);

% experiment 1
if parallel_runs
    [data, extra] = EM2005(startpar, 1, debug_mode, false, runs);
end
errors_exp1 = zeros(runs, 1);
for run = 1:runs
    if parallel_runs
        errors_exp1(run) = compute_err_exp1(data{run}, extra);
    else
        % serial runs
        fprintf('                            ... run = %d\n', run);
        [data, extra] = EM2005(startpar, 1, debug_mode, false, 1);
        errors_exp1(run) = compute_err_exp1(data{1}, extra);
    end
end
error_exp1 = mean(errors_exp1);
error_exp1_std = std(errors_exp1);
fprintf('\n                   .............. AVERAGE ERROR = %.4f (std = %.4f)\n', error_exp1, error_exp1_std);

% experiment 2
% TODO UNDO
%[data, extra] = EM2005(startpar, 2, debug_mode, false);
error_exp2 = 0; %error_exp2 = compute_err_exp2(data, extra);

error = error_exp1 + error_exp2;
%fprintf('\n\n\n     ======================>>>>> fitting [');
%fprintf('%8.2f', free_params);
%fprintf(', error = %f * 3 + %f = %f\n\n\n', error_exp1, error_exp2, error);

end

