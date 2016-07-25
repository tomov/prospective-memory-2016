function [data, extra, run_ids] = EM2005_condition(params, exp_id, condition, subjects_per_condition, blocks_per_condition, trials_per_block, debug_mode, fitting_mode, do_print, pm_blocks_exp1, pm_trials_exp2, pm_trials_exp3, pm_trials_exp6)
% Run a bunch of subjects in a given condition (eg Nonfocal, high emph PM)
%
run = condition(1);
OG_ONLY = condition(2);
FOCAL = condition(3);
EMPHASIS = condition(4);
TARGETS = condition(5);

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

data = []; 
extra = [];
run_ids = [];

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
og_correct_pattern = {'Yes'; 'No'; 'No'; 'Yes'; 'No'; 'Yes'};

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

% generate trial sequence (all blocks concatenated)
stimuli = repmat(og_block, blocks_per_condition, 1);
correct = repmat(og_block_correct, blocks_per_condition, 1);
og_correct = correct;
is_target = zeros(blocks_per_condition * trials_per_block, 1);
is_inter_task = [];

if fitting_mode
    assert(exp_id ~= 5);

    % when fitting, have PM task more often so we can get a
    % good estimate of the PM hit rate with less subjects
    %
    reps = blocks_per_condition * trials_per_block;

    stimuli_pattern = og_stimuli_pattern;
    correct_pattern = og_correct_pattern;
    is_target_pattern = zeros(length(og_stimuli_pattern), 1);
    if ~OG_ONLY
        stimuli_pattern([3 6], :) = pm_targets_pattern;
        correct_pattern([3 6]) = pm_correct_pattern;
        is_target_pattern([3 6]) = 1;
    end

    % repeat patterns
    stimuli = repmat(stimuli_pattern, reps, 1);
    correct = repmat(correct_pattern, reps, 1);
    og_correct = repmat(og_correct_pattern, reps, 1);
    is_target = repmat(is_target_pattern, reps, 1);

    % truncate
    stimuli = stimuli(1:reps, :);
    correct = correct(1:reps, :);
    og_correct = og_correct(1:reps, :);
    is_target = is_target(1:reps, :);

else % if not fitting_mode i.e. regular simulations
    % insert the PM targets
    %
    if ~OG_ONLY
        if debug_mode                        
            % in debug mode (but not fitting mode), every third trial is a PM trial -- this is only for
            % testing; not used in any of E&M's experiments
            %
            for i = 1:length(stimuli)
                if mod(i,40) == 0
                    target_id = mod(i, size(pm_targets_pattern, 1)) + 1;
                    middle = i;
                    stimuli(middle,:) = pm_targets_pattern(target_id, :);
                    correct(middle) = pm_correct_pattern(target_id);
                    og_correct(middle) = pm_og_correct_pattern(target_id);
                    is_target(middle) = 1;
                end
            end
        else % if not debug_mode nor fitting_mode
            if exp_id == 1
                % insert one PM target in each of the PM blocks
                % in experiment 1, there is a target in blocks 1, 3, 6, 7
                for i = 1:length(pm_blocks_exp1)
                    b = pm_blocks_exp1(i);
                    block_start = (b - 1) * trials_per_block + 1;
                    block_end = b * trials_per_block;
                    middle = int32((block_start + block_end) / 2);
                    target_id = mod(i, size(pm_targets_pattern, 1)) + 1;

                    stimuli(middle,:) = pm_targets_pattern(target_id, :);
                    correct(middle) = pm_correct_pattern(target_id);
                    og_correct(middle) = pm_og_correct_pattern(target_id);
                    is_target(middle) = 1;
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
                    stimuli(trial,:) = pm_targets_pattern(target_id, :);
                    correct(trial) = pm_correct_pattern(target_id);
                    og_correct(trial) = pm_og_correct_pattern(target_id);
                    is_target(trial) = 1;                        
                end
            end
        end % if debug_mode / else
    end % if ~OG_ONLY

    if exp_id == 5
        % experiment 5 is special altogether
        % Note that the "inter task" is fake -- it's
        % actually just the OG task without the PM
        % instruction. This is fine b/c there's no priming
        % in our model.
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

        stimuli = repmat(stimuli_pattern, reps, 1);
        correct = repmat(correct_pattern, reps, 1);
        og_correct = repmat(og_correct_pattern, reps, 1);
        is_target = repmat(is_target_pattern, reps, 1);
        is_or_was_target = repmat(is_or_was_target_pattern, reps, 1);
        is_inter_task = repmat(is_inter_task_pattern, reps, 1);

        % truncate
        stimuli = stimuli(1:reps, :);
        correct = correct(1:reps, :);
        og_correct = og_correct(1:reps, :);
        is_target = is_target(1:reps, :);
        is_or_was_target = is_or_was_target(1:reps, :);
        is_inter_task = is_inter_task(1:reps, :);
    else
        is_or_was_target = zeros(); % hacks to make parfor work
    end % if exp_id == 5
end % if fitting_mode / else
assert(length(correct) == length(stimuli));
assert(length(is_target) == length(stimuli));

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

model_params = zeros(1,6);
model_params(7) = 1; % init PM context
model_params(8) = 0; % init Other context
model_params(5) = bias_for_task;
model_params(6) = bias_for_attention;
model_params(9) = bias_for_context;
if exp_id == 6 && EMPHASIS == 0
    % Brewer et al. 2010, low WM capacity
    % note that we use EMPHASIS to denote high/low wm
    % capacity #hacksauce
    %
    model_params(5) = params(26); % bias for task
    model_params(6) = params(27); % bias for attention
    model_params(9) = params(28); % bias for context
end
model_params(10) = gamma;
model_params(11) = noise_sigma_ffwd;
model_params(12) = noise_sigma_wm;
if OG_ONLY
    model_params(1:4) = [1 0 1 0];
else       
    if FOCAL
        if ~EMPHASIS
            % focal, low emphasis
            model_params(1:4) = focal_low_init_wm;
        else
            % focal, high emphasis
            model_params(1:4) = focal_high_init_wm;
        end
    else
        if ~EMPHASIS
            % nonfocal, low emphasis
            model_params(1:4) = nonfocal_low_init_wm;
        else
            % nonfocal, high emphasis
            model_params(1:4) = nonfocal_high_init_wm;
        end
    end
end

% personalize perception & response monitoring for each subject, optionally
% adding cross-subject variability
%
assert(exp_id ~= 5); % TODO WTF FIXME make 5 work (split 1k in buckets), also fix it in the Simulator (see switch_to_Inter_task stuff -- BIG TODO)
% Parameters we add noise to across subjects, in order:
% 1. initial PM Task activation
% 2. initial PM target activation
% 3. WM tasks bias
% 4. WM attention bias
% 5. WM context bias
which_params_we_vary_across_subjects = [2 4 5 6 9];
subject_params = repmat(model_params(which_params_we_vary_across_subjects), subjects_per_condition, 1);
if ~OG_ONLY
    % PM task noise
    %subject_params(:, 1) = subject_params(:, 1) + unifrnd(-init_pm_task_noise_sigma, init_pm_task_noise_sigma, subjects_per_condition, 1);
    subject_params(:, 1) = subject_params(:, 1) + normrnd(0, init_pm_task_noise_sigma, size(subject_params(:, 1)));
    % PM task cannot be > OG task
    bad_ones = subject_params(:, 1) > model_params(1) - 0.07;
    subject_params(bad_ones, 1) = model_params(1) - 0.07;
    % PM target noise
    %subject_params(:, 2) = subject_params(:, 2) + unifrnd(-init_pm_target_noise_sigma, init_pm_target_noise_sigma, subjects_per_condition, 1);
    subject_params(:, 2) = subject_params(:, 2) + normrnd(0, init_pm_target_noise_sigma, size(subject_params(:, 2)));
    % PM target cannot be > OG features
    bad_ones = subject_params(:, 2) > model_params(3) - 0.07;
    subject_params(bad_ones, 2) = model_params(3) - 0.07;
end
% WM bias noise across subjects
% IMPORTANT -- make sure noise term is the same for all 3 biases for a given subject
which_biases = [3 4 5];
subject_params(:, which_biases) = subject_params(:, which_biases) + repmat(normrnd(0, wm_bias_noise_sigma, size(subject_params, 1), 1), 1, length(which_biases));

% initialize simulator (for multiple subjects)
%
sim = Simulator(FOCAL, model_params, subjects_per_condition, subject_params, fitting_mode);

if do_print
    for s=1:subjects_per_condition
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
            sim.instruction({'tortoise', 'dog', 'cat', 'kiwi', 'panda', 'monkey'}, true);
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

% run the actual simulations of all subjects
%
[responses, RTs, act, acc, onsets, offsets, nets] = sim.run(stimuli);

% collect the relevant data
%
if exp_id == 1 || exp_id == 3 || exp_id == 4 || exp_id == 5 || exp_id == 6
    % for experiment 1, each subject = 1 sample
    %
    [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_RT, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
        responses, RTs, act, acc, onsets, offsets, nets, ...
        is_target, correct, og_correct, is_inter_task, ...
        false, do_print);
    if exp_id == 5
        % extra analysis for experiment 5
        %
        it_targets = logical(is_or_was_target) & logical(is_inter_task);
        IT_TAR_RT = mean(RTs(it_targets));
        IT_TAR_SEM = std(RTs(it_targets)) / sqrt(length(RTs(it_targets)));
        it_non_targets = logical(~is_or_was_target) & logical(is_inter_task);
        IT_NONTAR_RT = mean(RTs(it_non_targets));
        IT_NONTAR_SEM = std(RTs(it_non_targets)) / sqrt(length(it_non_targets));
        if do_print, fprintf(' bonus Exp 5: target RT = %.2f (%.2f), nontarget RT = %.2f (%.2f)\n', ...
            IT_TAR_RT, IT_TAR_SEM, IT_NONTAR_RT, IT_NONTAR_SEM); end

        it_tar_resp = responses(it_targets);
        it_tar_correct = correct(it_targets);
        IT_TAR_HIT = sum(strcmp(it_tar_resp, it_tar_correct)) / length(it_tar_correct) * 100;
        if do_print, fprintf('            : accuracy on targets = %.2f\n', IT_TAR_HIT); end

        it_nontar_resp = responses(it_non_targets);
        it_nontar_correct = correct(it_non_targets);
        IT_NONTAR_HIT = sum(strcmp(it_nontar_resp, it_nontar_correct)) / length(it_nontar_correct) * 100;
        if do_print, fprintf('            : accuracy on non-targets = %.2f\n', IT_NONTAR_HIT); end
    end

    for s = 1:subjects_per_condition
        subject = [OG_ONLY, FOCAL, EMPHASIS, OG_RT(s,:)', OG_Hit(s,:)', PM_RT(s,:)', PM_Hit(s,:)', PM_miss_OG_hit(s,:)', TARGETS, first_PM_RT(s), PM_miss_OG_RT(s,:)'];
        if exp_id == 5
            subject = [subject, IT_TAR_RT, IT_NONTAR_RT, IT_TAR_HIT, IT_NONTAR_HIT];
        end
        data = [data; subject];
        run_ids = [run_ids; run];

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
    % for experiment 2, each block = 1 sample (i.e. 4
    % samples per subject)
    %
    for block_id = 1:blocks_per_condition
        block_start = (block_id - 1) * trials_per_block + 1;
        block_end = block_id * trials_per_block;                    
        [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_RT, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
            responses(:, block_start:block_end), RTs(:, block_start:block_end), [], [], [], [], [], ...
            is_target(block_start:block_end), ...
            correct(block_start:block_end), ...
            og_correct(block_start:block_end), ...
            [], ...
            false, do_print);

        % put subject and block id's at the end to make it
        % compatible with the data from experiment 1
        %
        for s = 1:subjects_per_condition
            block = [OG_ONLY, FOCAL, EMPHASIS, OG_RT(s,:)', OG_Hit(s,:)', PM_RT(s,:)', PM_Hit(s,:)', PM_miss_OG_hit(s,:)', s, block_id, first_PM_RT(s), PM_miss_OG_RT(s,:)'];
            data = [data; block];
            run_ids = [run_ids; run];

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
    for s = 1:subjects_per_condition
        temp_params = model_params;
        temp_params(which_params_we_vary_across_subjects) = subject_params(s,:);
        fprintf('   curpar(1:4) = %.3f %.3f %.3f %.3f\n', temp_params(1), temp_params(2), temp_params(3), temp_params(4));
        %if ~OG_ONLY
            getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                responses, RTs, act, acc, onsets, offsets, nets, ...
                is_target, correct, og_correct, is_inter_task, ...
                true, true);
        %end
    end
end
