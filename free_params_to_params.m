function params = free_params_to_params( free_params )
% takes only the free params that we fit with FMINCON
% returns the full set of params you can pass to EM2005

startpar = [1  0.35   1    0.3, ...     % focal, low emph     % exp1_v16, exp2_v19
            1  0.6    1    0.4, ...     % focal, high emph      % exp1_v16
            1  0.8    1    0.75, ...    % nonfocal, low emph   % exp2_v11
            1  0.9    1    0.83, ...    % nonfocal, high emph  % exp1_v16 -- sorta
            4 4 4, ... % biases -- tasks, attention, context
            0.3 0.3, ... % cross-subject init wm noise sigma -- PM task, target
            0.0004,  ... % gamma
            0.1 0.01, ...   % ffwd noise, wm noise sigma
            0.0];       % wm bias noise sigma

% COPY-PASTE in solve_exp1_and_exp2.m, at the very end
startpar([2 4 6 8 10 12 14 16 22]) = free_params; % set the params we're fitting
startpar(22) = free_params(9) * 10^(-3);

end

