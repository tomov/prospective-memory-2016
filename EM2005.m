function [data, extra] = EM2005( params, exp_id, fitting_mode, debug_mode, do_print, experiment_runs)
% run a simulation of the E&M with certain parameters and spit out the data
% for all subjects


% create parallel pool
%
if ~debug_mode
    if (strfind(version('-date'), '2013')) % rondo lives in 2013
        if matlabpool('size') == 0
            matlabpool;
        end
        fprintf('num of 2013 parallel workers = %d\n', matlabpool('size'));
    else
        poolobj = gcp;
        fprintf('num of parallel workers = %d\n', poolobj.NumWorkers);
    end
end

% parse parameters

if do_print, params
end

assert(exp_id == 1 || exp_id == 2 || exp_id == 3 || exp_id == 4 || exp_id == 5 || exp_id == 6);

if do_print, fprintf('\n\n--------========= RUNNING E&M EXPERIMENT %d ======-------\n\n', exp_id); end

% from E&M Experiment 1 & 2 methods
subjects_per_condition = [24 24 32 104 1000 30]; % experiment 5 is 72 subjects but that's not significant...
blocks_per_condition = [8 4 1 1 10 1];
trials_per_block = [24 40 110 110 18 110];
pm_blocks_exp1 = [1 3 6 7];
pm_trials_exp2 = [40 80 120 160];
pm_trials_exp3 = [26 52 78 104];
pm_trials_exp6 = [25 50 75 100];

% since we're doing only 1 experiment at a time
blocks_per_condition = blocks_per_condition(exp_id);
trials_per_block = trials_per_block(exp_id);
subjects_per_condition = subjects_per_condition(exp_id);

og_range = 0:1;
focal_range = 1:-1:0;
emphasis_range = 0:1;
target_range = [1, 6];

if exp_id == 1
    target_range = 1;
elseif exp_id == 2
    emphasis_range = 0;
    target_range = 1;
elseif exp_id == 3
    emphasis_range = 0;
    focal_range = 1;
elseif exp_id == 4
    focal_range = 1;
    target_range = 1;
    emphasis_range = 0;
elseif exp_id == 5
    og_range = 0; % TODO it's hardcoded -- there must be a PM task
    focal_range = 1; % TODO hardcoded too
    emphasis_range = 0;
    target_range = 1;
elseif exp_id == 6
    target_range = 1;
    % here emphasis means low/high wm capacity
end


if debug_mode
    % we show ~16 figures per subject. You don't want more than one subject
    %
    subjects_per_condition = 1;
    og_range = 1;
    focal_range = 1;
    emphasis_range = 1;
    target_range = [1];
    trials_per_block = 20;
    blocks_per_condition = 10;
elseif fitting_mode
    % when fitting, use less subjects for speed
    %
    if exp_id == 1 || exp_id == 2 || exp_id == 6
        subjects_per_condition = 4;
    else
        assert(exp_id == 3 || exp_id == 4);
        subjects_per_condition = 16;
    end
end

% List out all PM conditions
%
conditions = [];
og_only_condition = [];
for run = 1:experiment_runs
    for OG_ONLY = og_range
        if OG_ONLY
            % optimization -- the OG_ONLY case is essentially the same
            % regardless of the other conditions; so run only 1 simulation for
            % OG_ONLY and then resuse it for all other OG_ONLY conditions
            %
            og_only_condition = [OG_ONLY focal_range(1) emphasis_range(1) target_range(1)];
            conditions = [conditions; run og_only_condition];
        else
            for FOCAL = focal_range
                for EMPHASIS = emphasis_range
                    for TARGETS = target_range
                        conditions = [conditions; run OG_ONLY FOCAL EMPHASIS TARGETS];
                    end
                end
            end
        end
    end
end

data = []; 
extra = [];
run_ids = [];

% ----------> THIS IS WHERE IT'S AT <---------
% this is where we run the actual simulations...
% For each condition (in parallel), run all subjects through the
% experiment and append data
%
if ~debug_mode
    % parallelize like a boss
    %
    parfor cond_id = 1:size(conditions, 1)
        condition = conditions(cond_id, :);

        [data_for_cond, extra_for_cond, run_ids_for_cond] = ...
            EM2005_condition(params, exp_id, condition, subjects_per_condition, blocks_per_condition, trials_per_block, debug_mode, fitting_mode, do_print, pm_blocks_exp1, pm_trials_exp2, pm_trials_exp3, pm_trials_exp6);

        data = [data; data_for_cond];
        extra = [extra; extra_for_cond];
        run_ids = [run_ids; run_ids_for_cond];
    end
else
    % ...except in debug_mode (otherwise can't show figures)
    %
    for cond_id = 1:size(conditions, 1)
        condition = conditions(cond_id, :);

        [data_for_cond, extra_for_cond, run_ids_for_cond] = ...
            EM2005_condition(params, exp_id, condition, subjects_per_condition, blocks_per_condition, trials_per_block, debug_mode, fitting_mode, do_print, pm_blocks_exp1, pm_trials_exp2, pm_trials_exp3, pm_trials_exp6);

        data = [data; data_for_cond];
        extra = [extra; extra_for_cond];
        run_ids = [run_ids; run_ids_for_cond];
    end
end

% POST-SIMULATION stuffs
% resuse the OG_ONLY simulation for all conditions where OG_ONLY is true
%
for run = 1:experiment_runs
    for OG_ONLY = og_range
        if OG_ONLY
            og_only_rows = data(:, 1) == OG_ONLY & run_ids(:) == run;
            data_og_only = data(og_only_rows, :);
            extra_og_only = extra(og_only_rows, :);
            run_ids_og_only = run_ids(og_only_rows, :);
            for FOCAL = focal_range
                for EMPHASIS = emphasis_range
                    for TARGETS = target_range
                        if sum([OG_ONLY FOCAL EMPHASIS TARGETS] == og_only_condition) == length(og_only_condition)
                            continue % skip the one OG_ONLY condition we actually simulated
                        end
                        % set the other condition variables
                        data_og_only(:, 1) = OG_ONLY;
                        data_og_only(:, 2) = FOCAL;
                        data_og_only(:, 3) = EMPHASIS;
                        data_og_only(:, 9) = TARGETS;
                        % append the OG_ONLY responses as a new condition
                        data = [data; data_og_only];
                        extra = [extra; extra_og_only];
                        run_ids = [run_ids; run_ids_og_only];
                    end
                end
            end
        end
    end
end

data_reshaped = cell(experiment_runs, 1);
extra_reshaped = cell(experiment_runs, 1);
for run = 1:experiment_runs
    data_reshaped{run} = data(run_ids(:) == run, :);
    extra_reshaped{run} = extra(run_ids(:) == run, :);
end
data = data_reshaped;
extra = extra_reshaped;