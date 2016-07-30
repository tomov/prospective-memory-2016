function [data, extra] = EM2005( params, exp_id, debug_mode, do_print, experiment_runs, max_subjects_per_pool)
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
%

if do_print, params
end

focal_low_init_wm = params(1:4);
focal_high_init_wm = params(5:8);
nonfocal_low_init_wm = params(9:12);
nonfocal_high_init_wm = params(13:16);
bias_for_task = params(17);
bias_for_attention = params(18);
bias_for_context = params(19);
init_pm_task_noise_sigma = params(20);
init_pm_target_noise_sigma = params(21);
gamma = params(22);
noise_sigma_ffwd = params(23);
noise_sigma_wm = params(24);
wm_bias_noise_sigma = params(25);
og_weights_noise_factor = params(26);
bias_for_task_low_wm = params(27);
bias_for_attention_low_wm = params(28);
bias_for_context_low_wm = params(29);

assert(exp_id == 1 || exp_id == 2 || exp_id == 3 || exp_id == 4 || exp_id == 5 || exp_id == 6);

if do_print, fprintf('\n\n--------========= RUNNING E&M EXPERIMENT %d ======-------\n\n', exp_id); end

% from E&M Experiment 1 & 2 methods
subjects_per_condition = [24 24 32 104 72 30]; % experiment 5 is 72 subjects but that's not significant...
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

data = []; 
extra = [];
run_ids = [];
subject_pool_ids =[];

% Set up different experimental conditions
%

focal_range = 1:-1:0;
emphasis_range = 0:1;
target_range = [1, 6];
og_range = 0:1;

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
    % we show 1 figure per subject per condition. Don't go overboard
    %
    subjects_per_condition = 1;
    og_range = 0;
    focal_range = 0;
    emphasis_range = 0;
    target_range = [1];
    trials_per_block = 19;
    blocks_per_condition = 1;
end

% For each condition, split subjects into separate subject pools 
% which we can parallelize for faster computation. Special care must be taken
% in the end when we join the pools together for each condition.
%
n_subject_pools = double(idivide(int32(subjects_per_condition), max_subjects_per_pool, 'floor'));
subject_pool_sizes = ones(n_subject_pools, 1) * max_subjects_per_pool;
remainder = mod(subjects_per_condition, max_subjects_per_pool);
if remainder > 0
    n_subject_pools = n_subject_pools + 1;
    subject_pool_sizes = [subject_pool_sizes; remainder];
end
subject_pool_sizes


% List out all PM conditions. We parallelize them
%
conditions = [];
for run = 1:experiment_runs
    for FOCAL = focal_range
        for EMPHASIS = emphasis_range
            for TARGETS = target_range
                for subject_pool_id = 1:n_subject_pools
                    conditions = [conditions; run FOCAL EMPHASIS TARGETS subject_pool_id];
                end
            end
        end
    end
end

%
% Set up the sequence of stimuli and correct responses
% based on which experiment we're doing and other params
%

% init OG trial pool -- the 1 is the timeout in seconds
og_stimuli_pattern = [
    {'crocodile,an animal'}, 1;
    {'crocodile,a subject'}, 1;
    {'physics,an animal'}, 1;
    {'physics,a subject'}, 1;
    {'math,an animal'}, 1;
    {'math,a subject'}, 1;
];
og_correct_pattern = {'Yes'; 'No'; 'No'; 'Yes'; 'No'; 'Yes'}; % TODO rename 'correct' to 'expected' or something

% init PM trial pool
pm_targets_pattern = [
    {'tortoise,an animal'}, 1;
    {'tortoise,a subject'}, 1;
];
pm_og_correct_pattern = {'Yes'; 'No'};
pm_correct_pattern = {'PM'; 'PM'};

% generate OG block
og_block = repmat(og_stimuli_pattern, trials_per_block, 1);
og_block_correct = repmat(og_correct_pattern, trials_per_block, 1);
og_block = og_block(1:trials_per_block,:); % truncate
og_block_correct = og_block_correct(1:trials_per_block,:); % truncate

% generate trial sequence (all blocks concatenated) for both halves of the experiment
% cell index = OG_ONLY + 1
% note that initially, both halves look like the control half. We add the PM targets next
%
% PM half
stimuli = {repmat(og_block, blocks_per_condition, 1)};
correct = {repmat(og_block_correct, blocks_per_condition, 1)};
og_correct = correct;
is_target = {zeros(blocks_per_condition * trials_per_block, 1)};
is_inter_task = {[]};
% control (OG only) half
stimuli{2} = stimuli{1};
correct{2} = correct{1};
og_correct{2} = og_correct{1};
is_target{2} = is_target{1};
is_inter_task{2} = is_inter_task{1};

% insert the PM targets in the PM half of the experiment (note the OG half has 0 targets)
%
if debug_mode                        
    % in debug mode, PM trials are more frequent. This is only for
    % testing; not used in any of E&M's experiments
    %
    for i = 1:length(stimuli{1})
        if mod(i,5) == 0
            target_id = mod(i, size(pm_targets_pattern, 1)) + 1;
            middle = i;
            stimuli{1}(middle,:) = pm_targets_pattern(target_id, :);
            correct{1}(middle) = pm_correct_pattern(target_id);
            og_correct{1}(middle) = pm_og_correct_pattern(target_id);
            is_target{1}(middle) = 1;
        end
    end
else % these are based on E&M 2005
    if exp_id == 1
        % insert one PM target in each of the PM blocks
        % in experiment 1, there is a target in blocks 1, 3, 6, 7
        for i = 1:length(pm_blocks_exp1)
            b = pm_blocks_exp1(i);
            block_start = (b - 1) * trials_per_block + 1;
            block_end = b * trials_per_block;
            middle = int32((block_start + block_end) / 2);
            target_id = mod(i, size(pm_targets_pattern, 1)) + 1;

            stimuli{1}(middle,:) = pm_targets_pattern(target_id, :);
            correct{1}(middle) = pm_correct_pattern(target_id);
            og_correct{1}(middle) = pm_og_correct_pattern(target_id);
            is_target{1}(middle) = 1;
        end
    elseif exp_id == 2 || exp_id == 3 || exp_id == 4 || exp_id == 6
        % in experiment 2, trials 40, 80, 120, and 160 are
        % targets
        % experiment 3 also has 4 target trials
        if exp_id == 2
            pm_trials = pm_trials_exp2;
        elseif exp_id == 6
            pm_trials = pm_trials_exp6;
        else
            assert(exp_id == 3 || exp_id == 4)
            pm_trials = pm_trials_exp3;
        end
        for i = 1:length(pm_trials)
            target_id = mod(i, size(pm_targets_pattern, 1)) + 1;
            trial = pm_trials(i);
            stimuli{1}(trial,:) = pm_targets_pattern(target_id, :);
            correct{1}(trial) = pm_correct_pattern(target_id);
            og_correct{1}(trial) = pm_og_correct_pattern(target_id);
            is_target{1}(trial) = 1;                        
        end
    end
end % end inserting PM targets for exp. 1,2,3,4 and 6

% insert PM targets for exp. 5
%
if exp_id == 5
    % experiment 5 is special altogether
    % Note that the "inter task" is fake -- it's
    % actually just the OG task without the PM
    % instruction. This is fine b/c there's no priming
    % in our model.
    % TODO make num of targets same as in E&M (only 4 of the 10 blocks had targets)
    %
    stimuli_pattern = [
        {'switch back to OG and PM'}, 1; % task switch

        {'crocodile,an animal'}, 1;
        {'crocodile,a subject'}, 1;
        {'physics,an animal'}, 1;
        {'tortoise,an animal'}, 1; % PM target
        {'physics,a subject'}, 1;
        {'math,an animal'}, 1;
        {'math,a subject'}, 1;
        {'tortoise,a subject'}, 1; % PM target

        {'switch to Inter Task'}, 1; % task switch

        {'crocodile,an animal'}, 1;
        {'crocodile,a subject'}, 1;
        {'physics,an animal'}, 1;
        {'tortoise,an animal'}, 1; % formerly PM target
        {'physics,a subject'}, 1;
        {'math,an animal'}, 1;
        {'math,a subject'}, 1;
        {'tortoise,a subject'}, 1; % formerly PM target
    ];
    is_target_pattern = zeros(length(stimuli_pattern), 1);
    is_target_pattern([5 9]) = 1;
    is_or_was_target_pattern = zeros(length(stimuli_pattern), 1);
    is_or_was_target_pattern([5 9 14 18]) = 1;
    is_inter_task_pattern = [zeros(10, 1); ones(8, 1)]; % count switches as part of inter task
    og_correct_pattern = {'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; 'No'; 'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; 'No'};
    correct_pattern = og_correct_pattern;
    correct_pattern([5 9]) = {'PM'};

    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(correct_pattern));
    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(is_target_pattern));
    assert(length(stimuli_pattern) == length(is_or_was_target_pattern));
    assert(length(stimuli_pattern) == length(is_inter_task_pattern));

    % copy & trim 'em
    reps = blocks_per_condition * trials_per_block;

    stimuli{1} = repmat(stimuli_pattern, reps, 1);
    correct{1} = repmat(correct_pattern, reps, 1);
    og_correct{1} = repmat(og_correct_pattern, reps, 1);
    is_target{1} = repmat(is_target_pattern, reps, 1);
    is_or_was_target = repmat(is_or_was_target_pattern, reps, 1);
    is_inter_task{1} = repmat(is_inter_task_pattern, reps, 1);

    % truncate
    stimuli{1} = stimuli{1}(1:reps, :);
    correct{1} = correct{1}(1:reps, :);
    og_correct{1} = og_correct{1}(1:reps, :);
    is_target{1} = is_target{1}(1:reps, :);
    is_or_was_target = is_or_was_target(1:reps, :);
    is_inter_task{1} = is_inter_task{1}(1:reps, :);
else
    is_or_was_target = zeros(); % hacks to make parfor work
end % if exp_id == 5
assert(length(correct{1}) == length(stimuli{1}));
assert(length(is_target{1}) == length(stimuli{1}));


% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%

% ...THIS IS IT...
% For each condition (in parallel), run all subjects through the
% experiment
% TODO add more parallellism e.g. to run same simulation multiple times,
% just add the conditions several times
%
% PARFOR
parfor cond_id = 1:size(conditions, 1)
    condition = conditions(cond_id, :);
    run = condition(1);
    FOCAL = condition(2);
    EMPHASIS = condition(3);
    TARGETS = condition(4);
    subject_pool_id = condition(5);
    subject_pool_size = subject_pool_sizes(subject_pool_id);

    if do_print
        fprintf('--------------------------------- parfor condition %d %d %d %d %d (%d) ------------------------------\n', run, FOCAL, EMPHASIS, TARGETS, subject_pool_id, subject_pool_size);
    end

    % randomize order
    % ...not today tho
    %{
    idx = randperm(size(stimuli, 1))';
    stimuli = stimuli(idx, :);
    is_target = is_target(idx, :);
    correct = correct(idx, :);
    %}

    %
    % Get appropriate subject parameters depending on the condition
    %
    model_params = [];
    model_params(7) = 1; % init PM context
    model_params(8) = 0; % init Other context
    if exp_id == 6 && EMPHASIS == 0
        % Brewer et al. 2010, low WM capacity
        % note that we use EMPHASIS to denote high/low wm
        % capacity #hacksauce
        %
        model_params(5) = bias_for_task_low_wm;
        model_params(6) = bias_for_attention_low_wm;
        model_params(9) = bias_for_context_low_wm;
    else
        % Brewer et al. 2010, high WM capacity
        % ...or any other experiment
        %
        model_params(5) = bias_for_task;
        model_params(6) = bias_for_attention;
        model_params(9) = bias_for_context;
    end
    model_params(10) = gamma;
    model_params(11) = noise_sigma_ffwd;
    model_params(12) = noise_sigma_wm;
    model_params(13) = og_weights_noise_factor;

    % personalize perception & response monitoring for each subject, optionally
    % adding cross-subject variability
    %
    % Parameters that vary across subjects, in order:
    % 1. initial PM Task activation
    % 2. initial PM target activation
    % 3. WM tasks bias
    % 4. WM attention bias
    % 5. WM context bias
    which_params_we_vary_across_subjects = [2 4 5 6 9];
    subject_params_noise = [...
        normrnd(0, init_pm_task_noise_sigma, subject_pool_size, 1) ... % init PM task noise
        normrnd(0, init_pm_target_noise_sigma, subject_pool_size, 1) ... % init PM target noise
        repmat(normrnd(0, wm_bias_noise_sigma, subject_pool_size, 1), 1, 3) ... % WM bias noise IMPORTANT -- make sure noise term is the same for all 3 WM biases for a given subject
    ];

    % Run the PM half (OG_ONLY = 0) and control half (OG_ONLY = 1) of the experiment for each subject.
    % Note that, unlike E&M, we don't need to counterbalance the order here since
    % there are no priming effects and, more importantly, no PM targets in the control havles
    % => there can be no aftereffects of intention
    %
    for OG_ONLY = og_range

        % Readjust initial WM based on which half of the experiment we're in
        %
        if OG_ONLY
            model_params(1:4) = [1 0 1 0];
        else       
            if FOCAL
                if ~EMPHASIS
                    model_params(1:4) = focal_low_init_wm; % focal, low emphasis
                else
                    model_params(1:4) = focal_high_init_wm; % focal, high emphasis
                end
            else
                if ~EMPHASIS
                    model_params(1:4) = nonfocal_low_init_wm; % nonfocal, low emphasis
                else
                    model_params(1:4) = nonfocal_high_init_wm; % nonfocal, high emphasis
                end
            end
        end
       
        % Replicate parameters across the subjects and add cross-subject variability
        %
        subject_params = repmat(model_params(which_params_we_vary_across_subjects), subject_pool_size, 1);

        if ~debug_mode % no cross-subject noise in debug mode
            if ~OG_ONLY
                % add noise to initial WM activations across subjects
                subject_params(:, [1 2]) = subject_params(:, [1 2]) + subject_params_noise(:, [1 2]);
                % PM task cannot be > OG task (or too close to it)
                bad_ones = subject_params(:, 1) > model_params(1) - 0.07;
                subject_params(bad_ones, 1) = model_params(1) - 0.07;
                % PM target cannot be > OG features (or too close to it)
                bad_ones = subject_params(:, 2) > model_params(3) - 0.07;
                subject_params(bad_ones, 2) = model_params(3) - 0.07;
            end
            % add WM bias noise across subjects
            subject_params(:, [3 4 5]) = subject_params(:, [3 4 5]) + subject_params_noise(:, [3 4 5]);
        end

        % initialize simulator (for all subjects in the given condition)
        %
        sim = Simulator(FOCAL, model_params, subject_pool_size, subject_params);

        if do_print
            for s=1:subject_pool_size
                temp_params = model_params;
                temp_params(which_params_we_vary_across_subjects) = subject_params(s,:);
                fprintf('\nsubject %d: curpar(2,4) = %.2f %.2f curpar(5,6,9) = %.2f %.2f %.2f\n', s, temp_params(2), temp_params(4), temp_params(5), temp_params(6), temp_params(9));
            end
        end

        % PM instruction
        %
        if ~OG_ONLY
            if FOCAL
                if TARGETS == 6
                    % TODO FIXME DISCUSS -- we don't rely on holding all of
                    % these in WM (see Jon's comment in thesis). Instead, have
                    % 1 unit which means "Monitor for PM features" which
                    % coincide with the OG features -> they complement each
                    % other (just like the previous examples)
                    % STILL USEFUL e.g. FOR NONFOCAL -- there we would have ot
                    % keep them in WM
                    sim.instruction({'tortoise', 'dog', 'cat', 'kiwi', 'panda', 'monkey'}, true);
                    %sim.instruction({'tortoise'}, true);
                else
                    assert(TARGETS == 1);
                    %if exp_id == 5 TODO cleanup
                    %    sim.instruction({'tortoise'}, false);
                    %else
                        sim.instruction({'tortoise'}, true);
                    %end
                end
            else
                sim.instruction({'tor'}, true);
            end
        end

        % run the actual simulations for all subjects in the given condition
        %
        [responses, RTs, act, acc, onsets, offsets, nets] = sim.run(stimuli{OG_ONLY + 1}, false);
        
        % collect the relevant data
        %
        if exp_id == 1 || exp_id == 3 || exp_id == 4 || exp_id == 5 || exp_id == 6
            % for all experiments except for exp. 2, there are two samples per subject --
            % one for the control half (OG_ONLY = 1) and one for the PM half (OG_ONLY = 0) of the experiment.
            %
            [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_RT, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                responses, RTs, act, acc, onsets, offsets, nets, ...
                is_target{OG_ONLY + 1}, correct{OG_ONLY + 1}, og_correct{OG_ONLY + 1}, is_inter_task{OG_ONLY + 1}, ...
                false, do_print);

            
            % add each sample separately
            %
            for s = 1:subject_pool_size
                % IMPORTANT -- ordering here is critical. If you change stuff, you need to also change EM2005_with_stats*.m and B2010_with_stats.m
                %
                subject = [OG_ONLY, FOCAL, EMPHASIS, OG_RT(s,:)', OG_Hit(s,:)', PM_RT(s,:)', PM_Hit(s,:)', PM_miss_OG_hit(s,:)', TARGETS, first_PM_RT(s), PM_miss_OG_RT(s,:)'];

                if exp_id == 5
                    % extra analysis for experiment 5
                    %
                    assert(OG_ONLY == 0);
                    it_targets = logical(is_or_was_target) & logical(is_inter_task{1});
                    it_target_RTs = RTs(s, it_targets)'; % RT for (ex-)target items in the Inter task TODO only take the correct ones!
                    IT_TAR_RT = mean(it_target_RTs); 
                    IT_TAR_SEM = std(it_target_RTs) / sqrt(length(it_target_RTs));
                    
                    it_nontargets = logical(~is_or_was_target) & logical(is_inter_task{1});
                    it_nontarget_RTs = RTs(s, it_nontargets)';  % RT for (ex-)nontarget items in the Inter task TODO normalize! 'neural' items and 'targets' should be same # (for SD's)
                    IT_NONTAR_RT = mean(it_nontarget_RTs);
                    IT_NONTAR_SEM = std(it_nontarget_RTs) / sqrt(length(it_nontarget_RTs));
                    if do_print, fprintf(' bonus Exp 5: target RT = %.2f (%.2f), nontarget RT = %.2f (%.2f)\n', ...
                        IT_TAR_RT, IT_TAR_SEM, IT_NONTAR_RT, IT_NONTAR_SEM); end

                    it_tar_resp = responses(s, it_targets)';
                    it_tar_correct = correct{1}(it_targets); % TODO rename to it_tar_expected or something
                    num_correct_IT_responses_on_targets = sum(strcmp(it_tar_resp, it_tar_correct));
                    IT_TAR_HIT = num_correct_IT_responses_on_targets / length(it_tar_correct) * 100; % accuracy on (ex-)target items in the Inter task
                    if do_print, fprintf('            : accuracy on targets = %.2f\n', IT_TAR_HIT); end

                    IT_TAR_PM_HIT = sum(strcmp(it_tar_resp, 'PM')) / (length(it_tar_correct) - num_correct_IT_responses_on_targets) * 100; % how many of the wrong answers on the (ex-)targets were PM responses
                    if do_print, fprintf('            : PM hits on wrong-answer targets = %.2f\n', IT_TAR_PM_HIT); end

                    it_nontar_resp = responses(s, it_nontargets)';
                    it_nontar_correct = correct{1}(it_nontargets);
                    IT_NONTAR_HIT = sum(strcmp(it_nontar_resp, it_nontar_correct)) / length(it_nontar_correct) * 100; % accuracy on (ex-)nontarget items in the Inter task
                    if do_print, fprintf('            : accuracy on non-targets = %.2f\n', IT_NONTAR_HIT); end

                    subject = [subject, IT_TAR_RT, IT_TAR_HIT, IT_NONTAR_RT, IT_NONTAR_HIT, IT_TAR_PM_HIT];
                end

                data = [data; subject];
                run_ids = [run_ids; run];
                subject_pool_ids = [subject_pool_ids; subject_pool_id];

                temp_params = model_params;
                temp_params(which_params_we_vary_across_subjects) = subject_params(s,:);
                if debug_mode
                    subject_extra = {sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, responses(s,:)', RTs(s,:)', squeeze(act(s,:,:)), squeeze(acc(s,:,:)), onsets(s,:)', offsets(s,:)', squeeze(nets(s,:,:)), temp_params};
                    extra = [extra; subject_extra];
                else
                    extra = [extra; temp_params];
                end
            end
        elseif exp_id == 2
            % for experiment 2, there are 8 (= 2 x 4) samples per subject -- one sample for each of the four blocks, and
            % the two halves of the experiment -- the control half (OG_ONLY = 1) and PM half (OG_ONLY = 0)
            %
            for block_id = 1:blocks_per_condition
                block_start = (block_id - 1) * trials_per_block + 1;
                block_end = block_id * trials_per_block;                    
                [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_RT, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                    responses(:, block_start:block_end), RTs(:, block_start:block_end), [], [], [], [], [], ...
                    is_target{OG_ONLY + 1}(block_start:block_end), ...
                    correct{OG_ONLY + 1}(block_start:block_end), ...
                    og_correct{OG_ONLY + 1}(block_start:block_end), ...
                    [], ...
                    false, do_print);

                % put subject and block id's at the end to make it
                % compatible with the data from experiment 1
                %
                for s = 1:subject_pool_size
                    block = [OG_ONLY, FOCAL, EMPHASIS, OG_RT(s,:)', OG_Hit(s,:)', PM_RT(s,:)', PM_Hit(s,:)', PM_miss_OG_hit(s,:)', s, block_id, first_PM_RT(s), PM_miss_OG_RT(s,:)'];
                    data = [data; block];
                    run_ids = [run_ids; run];
                    subject_pool_ids = [subject_pool_ids; subject_pool_id];

                    temp_params = model_params;
                    temp_params(which_params_we_vary_across_subjects) = subject_params(s,:);
                    if debug_mode
                        subject_extra = {sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, responses(s,:)', RTs(s,:)', squeeze(act(s,:,:)), squeeze(acc(s,:,:)), onsets(s,:)', offsets(s,:)', squeeze(nets(s,:,:)), s, block, temp_params};
                        extra = [extra; subject_extra];
                    else
                        extra = [extra; temp_params];
                    end
                end
            end
        end % if exp_id == 1, 3, 4, 5, 6 / elseif exp_id == 2

        % show picture of whole thing (for debugging)
        % NOTE: doesn't work with parfor!! need regular forloop
        %
        if debug_mode
            for s = 1:subject_pool_size
                temp_params = model_params;
                temp_params(which_params_we_vary_across_subjects) = subject_params(s,:);
                fprintf('   curpar(1:4) = %.3f %.3f %.3f %.3f\n', temp_params(1), temp_params(2), temp_params(3), temp_params(4));
            end
            getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                responses, RTs, act, acc, onsets, offsets, nets, ...
                is_target{OG_ONLY + 1}, correct{OG_ONLY + 1}, og_correct{OG_ONLY + 1}, is_inter_task{OG_ONLY + 1}, ...
                true, true);
        end

    end % for OG_ONLY = [0 1]

end % for condition = conditions


% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%
% -----------------------------------------------------------------------------------------%

%
% Post-simulation stuffs
%

% sanity check -- % ensure the ordering by checking things that vary across subjects always (even in
% the OG_ONLY half), e.g. the WM bias
% ...mostly relevant for experiment 4

% Since we split the subject pool into several smaller subject pools (for faster evaluation),
% we now have to merge them. This is a bit tricky --
% remember, because of the way we do things, for each condition,
% the OG_ONLY = 0 and OG_ONLY = 1 entries should correspond to the two
% experiment halves for same sequences of subjects.
% But this order is not guarandeed across the conditions, and not even across the different subject pools
% => we reorder them here and then make sure we didn't screw up the ordering
%
data_pools_joined = [];
extra_pools_joined = [];
for subject_pool_id = 1:n_subject_pools
    data_pools_joined = [data_pools_joined; data(subject_pool_ids(:) == subject_pool_id, :)];
    extra_pools_joined = [extra_pools_joined; extra(subject_pool_ids(:) == subject_pool_id, :)];
end

% sanity check -- because of the way we do things, for each condition,
% the OG_ONLY = 0 and OG_ONLY = 1 entries should correspond to the two
% experiment halves for same sequences of subjects.
% ensure the ordering by checking things that vary across subjects always (even in
% the OG_ONLY half), e.g. the WM bias
% ...mostly relevant for experiment 4
%
if exp_id ~= 5 % no OG_ONLY half in experiment 5
    for cond_id = 1:size(conditions, 1)
        condition = conditions(cond_id, :);
        run = condition(1);
        FOCAL = condition(2);
        EMPHASIS = condition(3);
        TARGETS = condition(4);
        
        wm_bias_params = [5 6 9];
        og_wm_biases = extra(run_ids(:) == run & data(:, 1) == 1 & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS & data(:, 9) == TARGETS, wm_bias_params);
        pm_wm_biases = extra(run_ids(:) == run & data(:, 1) == 0 & data(:, 2) == FOCAL & data(:, 3) == EMPHASIS & data(:, 9) == TARGETS, wm_bias_params);
        assert(sum(sum(og_wm_biases - pm_wm_biases)) == 0);
    end
end

% reshape the data so its a cell array, one cell for each run (only relevant when fitting)
%
data_reshaped = cell(experiment_runs, 1);
extra_reshaped = cell(experiment_runs, 1);
for run = 1:experiment_runs
    data_reshaped{run} = data(run_ids(:) == run, :);
    extra_reshaped{run} = extra(run_ids(:) == run, :);
end
data = data_reshaped;
extra = extra_reshaped;
