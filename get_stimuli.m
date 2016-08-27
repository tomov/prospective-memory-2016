function [stimuli, correct, og_correct, is_target, is_or_was_target, is_nontarget, is_inter_task] = get_stimuli(exp_id, trials_per_block, blocks_per_condition, debug_mode)
% Return the trial sequence as a list of stimuli and responses + appropriate flags for each trial
%

% TODO move to constants file
pm_blocks_exp1 = [1 3 6 7];
pm_trials_exp2 = [40 80 120 160];
pm_trials_exp3 = [26 52 78 104];
pm_trials_exp6 = [25 50 75 100];

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
        if mod(i,48) == 0
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
% note this is the same in debug_mode -- experiment 5 is a bit special
%
if exp_id == 5
    % experiment 5 is special altogether
    % Note that the "inter task" is fake -- it's
    % actually just the OG task without the PM
    % instruction. This is fine b/c there's no priming
    % in our model.
    % Ref. Table 3 from E&M 2005 (not exactly the same but close enough)
    %
    stimuli_pattern = [
        {'switch back to OG and PM'}, 1; % do the OG + PM task ("Imagery rating" in E&M)
        {'physics,an animal'}, 1;
        {'crocodile,a subject'}, 1;
        {'crocodile,an animal'}, 1;
        {'tortoise,an animal'}, 1; % PM target
        {'physics,a subject'}, 1;
        {'math,an animal'}, 1;
        {'math,a subject'}, 1;

        {'switch to Inter Task'}, 1; % do the Inter task ("Lexical decision task" in E&M)
        {'physics,an animal'}, 1;          % high 
        {'crocodile,a subject'}, 1;         % low 
        {'crocodile,an animal'}, 1;         % high
        {'crocodile,an animal'}, 1; % ARGH  % low 
        {'physics,a subject'}, 1;           % low
        {'math,an animal'}, 1;              % high
        {'math,a subject'}, 1;               % high

        {'physics,an animal'}, 1;           % high
        {'crocodile,a subject'}, 1;         % low
        {'crocodile,an animal'}, 1;         % high
        {'tortoise,an animal'}, 1; % PM target  % low
        {'physics,a subject'}, 1;               % low
        {'math,an animal'}, 1;               % high
        {'math,a subject'}, 1;               % high
        
        {'switch back to OG and PM'}, 1; % do the OG + PM task ("Imagery rating" in E&M)
        {'physics,an animal'}, 1;
        {'crocodile,a subject'}, 1;
        {'crocodile,an animal'}, 1;
        {'tortoise,an animal'}, 1; % PM target
        {'physics,a subject'}, 1;
        {'math,an animal'}, 1;
        {'math,a subject'}, 1;
    ];
    is_target_pattern = zeros(length(stimuli_pattern), 1);
    is_target_pattern([5 28]) = 1;
    is_or_was_target_pattern = zeros(length(stimuli_pattern), 1);
    is_or_was_target_pattern([5 20 28]) = 1;
    % pick the same number of non-target items for the analysis (in this
    % case, the trials right before the target trial in the inter task).
    % since we don't have priming, the concept of "previously presented items" (as in E&M) is irrelevant here
    is_nontarget_pattern = zeros(length(stimuli_pattern), 1);
    is_nontarget_pattern([13]) = 1; 
    is_inter_task_pattern = [zeros(1 + 7, 1); ones(1 + 7 + 7, 1); zeros(1 + 7, 1)]; % count switches as part of inter task
    og_correct_pattern = { ...
        'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; ...
        'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; ...
                  'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; ...
        'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; };
    correct_pattern = og_correct_pattern;
    correct_pattern([5 28]) = {'PM'};

    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(correct_pattern));
    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(is_target_pattern));
    assert(length(stimuli_pattern) == length(is_or_was_target_pattern));
    assert(length(stimuli_pattern) == length(is_nontarget_pattern));
    assert(length(stimuli_pattern) == length(is_inter_task_pattern));

    % copy & trim 'em
    reps = blocks_per_condition * trials_per_block;

    stimuli{1} = repmat(stimuli_pattern, reps, 1);
    correct{1} = repmat(correct_pattern, reps, 1);
    og_correct{1} = repmat(og_correct_pattern, reps, 1);
    is_target{1} = repmat(is_target_pattern, reps, 1);
    is_or_was_target = repmat(is_or_was_target_pattern, reps, 1);
    is_nontarget = repmat(is_nontarget_pattern, reps, 1);
    is_inter_task{1} = repmat(is_inter_task_pattern, reps, 1);

    % truncate
    stimuli{1} = stimuli{1}(1:reps, :);
    correct{1} = correct{1}(1:reps, :);
    og_correct{1} = og_correct{1}(1:reps, :);
    is_target{1} = is_target{1}(1:reps, :);
    is_or_was_target = is_or_was_target(1:reps, :);
    is_nontarget = is_nontarget(1:reps, :);
    is_inter_task{1} = is_inter_task{1}(1:reps, :);
elseif exp_id == 7
    % this is only for predictions
    %
    stimuli_pattern = [
        {'switch back to OG and PM'}, 1; % do the OG + PM task ("Imagery rating" in E&M)
        {'physics,an animal'}, 1;
        {'crocodile,a subject'}, 1;
        {'crocodile,an animal'}, 1;
        {'tortoise,an animal'}, 1; % PM target
        {'physics,a subject'}, 1;
        {'math,an animal'}, 1;
        {'math,a subject'}, 1;

        {'switch to Inter Task'}, 1; % do the Inter task ("Lexical decision task" in E&M)
        {'physics,an animal'}, 1;          % high 
        {'crocodile,a subject'}, 1;         % low 
        {'crocodile,an animal'}, 1;         % high
        {'crocodile,an animal'}, 1; % ARGH  % low 
        {'physics,a subject'}, 1;           % low
        {'math,an animal'}, 1;              % high
        {'math,a subject'}, 1;               % high

        {'physics,an animal'}, 1;           % high
        {'crocodile,a subject'}, 1;         % low
        {'crocodile,an animal'}, 1;         % high
        {'tortoise,an animal'}, 1; % PM target  % low
        {'physics,a subject'}, 1;               % low
        {'math,an animal'}, 1;               % high
        {'math,a subject'}, 1;               % high
    ];
    is_target_pattern = zeros(length(stimuli_pattern), 1);
    is_target_pattern([5]) = 1;
    is_or_was_target_pattern = zeros(length(stimuli_pattern), 1);
    is_or_was_target_pattern([5 20]) = 1;
    is_nontarget_pattern = zeros(length(stimuli_pattern), 1);
    is_nontarget_pattern([13]) = 1; 
    is_inter_task_pattern = [zeros(1 + 7, 1); ones(1 + 7 + 7, 1)]; % count switches as part of inter task
    og_correct_pattern = { ...
        'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; ...
        'Switch'; 'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'; ...
                  'Yes'; 'No'; 'No'; 'Yes'; 'Yes'; 'No'; 'Yes'};
    correct_pattern = og_correct_pattern;
    correct_pattern([5]) = {'PM'};

    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(correct_pattern));
    assert(length(stimuli_pattern) == length(og_correct_pattern));
    assert(length(stimuli_pattern) == length(is_target_pattern));
    assert(length(stimuli_pattern) == length(is_or_was_target_pattern));
    assert(length(stimuli_pattern) == length(is_nontarget_pattern));
    assert(length(stimuli_pattern) == length(is_inter_task_pattern));

    % copy & trim 'em
    reps = blocks_per_condition * trials_per_block;

    stimuli{1} = [stimuli_pattern; repmat(stimuli_pattern(10:end,:), reps, 1)];
    correct{1} = [correct_pattern; repmat(correct_pattern(10:end), reps, 1)];
    og_correct{1} = [og_correct_pattern; repmat(og_correct_pattern(10:end), reps, 1)];
    is_target{1} = [is_target_pattern; repmat(is_target_pattern(10:end), reps, 1)];
    is_or_was_target = [is_or_was_target_pattern; repmat(is_or_was_target_pattern(10:end), reps, 1)];
    is_nontarget = [is_nontarget_pattern; repmat(is_nontarget_pattern(10:end), reps, 1)];
    is_inter_task{1} = [is_inter_task_pattern; repmat(is_inter_task_pattern(10:end), reps, 1)];

    % truncate
    stimuli{1} = stimuli{1}(1:reps, :);
    correct{1} = correct{1}(1:reps, :);
    og_correct{1} = og_correct{1}(1:reps, :);
    is_target{1} = is_target{1}(1:reps, :);
    is_or_was_target = is_or_was_target(1:reps, :);
    is_nontarget = is_nontarget(1:reps, :);
    is_inter_task{1} = is_inter_task{1}(1:reps, :);
else
    is_or_was_target = zeros(); % hacks to make parfor work
    is_nontarget = zeros(); % hack to make parfor work
end % if exp_id == 5; elseif exp_id == 7
assert(length(correct{1}) == length(stimuli{1}));
assert(length(is_target{1}) == length(stimuli{1}));
