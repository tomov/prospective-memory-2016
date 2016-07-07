function [data, extra] = EM2005( params, exp_id, fitting_mode, debug_mode, do_print)
% run a simulation of the E&M with certain parameters and spit out the data
% for all subjects


% create parallel pool
%
if (strfind(version('-date'), '2013')) % rondo lives in 2013
    if matlabpool('size') == 0
        matlabpool;
    end
    fprintf('num of 2013 parallel workers = %d\n', matlabpool('size'));
else
    poolobj = gcp;
    fprintf('num of parallel workers = %d\n', poolobj.NumWorkers);
end


% parse parameters

if do_print, params
end

focal_low_init_wm = params(1:4);
focal_high_init_wm = params(5:8);
nonfocal_low_init_wm = params(9:12);
nonfocal_high_init_wm = params(13:16);
bias_for_task = params(17);
bias_for_attention = params(18);
bias_for_context = params(19);
param_noise_sigma_1 = params(20);
param_noise_sigma_2 = params(21);
gamma = params(22);

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

data = [];
extra = [];

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
    og_range = 0;
    %focal_range = 1;
    %emphasis_range = 0;
    %target_range = [1,6];
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


for OG_ONLY = og_range
    for FOCAL = focal_range
        for EMPHASIS = emphasis_range
            for TARGETS = target_range

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
                                if mod(i,4) == 0
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
                %
                %{
                idx = randperm(size(stimuli, 1))';
                stimuli = stimuli(idx, :);
                is_target = is_target(idx, :);
                correct = correct(idx, :);
                %}

                % get appropriate parameters depending on the condition
                %
                curpar = zeros(1,6);
                curpar(7) = 1;
                curpar(8) = 0;
                curpar(5) = bias_for_task;
                curpar(6) = bias_for_attention;
                curpar(9) = bias_for_context;
                if exp_id == 6 && EMPHASIS == 0
                    % Brewer et al. 2010, low WM capacity
                    % note that we use EMPHASIS to denote high/low wm
                    % capacity #hacksauce
                    %
                    curpar(5) = params(23); % bias for task
                    curpar(6) = params(24); % bias for attention
                    curpar(9) = params(25); % bias for context
                end
                curpar(10) = gamma;
                if OG_ONLY
                    curpar(1:4) = [1 0 1 0];
                else       
                    if FOCAL
                        if ~EMPHASIS
                            % focal, low emphasis
                            curpar(1:4) = focal_low_init_wm;
                        else
                            % focal, high emphasis
                            curpar(1:4) = focal_high_init_wm;
                        end
                    else
                        if ~EMPHASIS
                            % nonfocal, low emphasis
                            curpar(1:4) = nonfocal_low_init_wm;
                        else
                            % nonfocal, high emphasis
                            curpar(1:4) = nonfocal_high_init_wm;
                        end
                    end
                end

                % simulate subjects in parallel; must be serial in
                % debug_mode (i.e. regular for)
                %
                for subject_id = 1:subjects_per_condition
                    % optionally add cross-subject variability
                    %
                    subjpar = curpar;
                    if ~OG_ONLY
                        subjpar(2) = subjpar(2) + unifrnd(-param_noise_sigma_1, param_noise_sigma_1);
                        subjpar(4) = subjpar(4) + unifrnd(-param_noise_sigma_1, param_noise_sigma_2);
                    end
                    
                    % initialize simulator             
                    %
                    sim = Simulator(FOCAL, subjpar, fitting_mode);
                
                    if do_print, fprintf('\nsubject %d: curpar(2,4) = %.2f %.2f curpar(5,6,9) = %.2f %.2f %.2f\n', subject_id, subjpar(2), subjpar(4), subjpar(5), subjpar(6), subjpar(9)); end
                    
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
                    
                    % run the actual simulation
                    %
                    [responses, RTs, act, acc, onsets, offsets, nets] = sim.trial(stimuli);

                    save('shitface.mat');
                    sdfsdf;

                    % collect the relevant data
                    %
                    if exp_id == 1 || exp_id == 3 || exp_id == 4 || exp_id == 5 || exp_id == 6
                        % for experiment 1, each subject = 1 sample
                        %
                        [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_hit, first_PM_RT] = getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                            responses, RTs, act, acc, onsets, offsets, ...
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
                     
                        subject = [OG_ONLY, FOCAL, EMPHASIS, OG_RT, OG_Hit, PM_RT, PM_Hit, PM_miss_OG_hit, TARGETS, first_PM_RT];
                        if exp_id == 5
                            subject = [subject, IT_TAR_RT, IT_NONTAR_RT, IT_TAR_HIT, IT_NONTAR_HIT];
                        end
                        data = [data; subject];
                        
                        if debug_mode
                            subject_extra = {sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, responses, RTs, act, acc, onsets, offsets, nets, subjpar};
                            extra = [extra; subject_extra];
                        else
                            extra = [extra; subjpar];
                        end
                    elseif exp_id == 2
                        % for experiment 2, each block = 1 sample (i.e. 4
                        % samples per subject)
                        %
                        for block_id = 1:blocks_per_condition
                            block_start = (block_id - 1) * trials_per_block + 1;
                            block_end = block_id * trials_per_block;                    
                            [OG_RT, ~, OG_Hit, PM_RT, ~, PM_Hit, PM_miss_OG_hit, ~] = ...
                                getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                                responses(block_start:block_end), RTs(block_start:block_end), [], [], [], [], ...
                                is_target(block_start:block_end), ...
                                correct(block_start:block_end), ...
                                og_correct(block_start:block_end), ...
                                [], ...
                                false, do_print);

                            % put subject and block id's at the end to make it
                            % compatible with the data from experiment 1
                            %
                            block = [OG_ONLY, FOCAL, EMPHASIS, OG_RT, OG_Hit, PM_RT, PM_Hit, PM_miss_OG_hit, subject_id, block_id];
                            data = [data; block];
                            if debug_mode
                                subject_extra = {sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, responses, RTs, act, acc, onsets, offsets, nets, subject_id, block, subjpar};
                                extra = [extra; subject_extra];
                            else
                                extra = [extra; subjpar];
                            end
                        end
                    end % if exp_id / elseif exp_id
                    
                    % show picture of whole thing (for debugging)
                    % NOTE: doesn't work with parfor!! need regular forloop
                    %
                    if debug_mode
                        fprintf('   curpar(1:4) = %.3f %.3f %.3f %.3f\n', subjpar(1), subjpar(2), subjpar(3), subjpar(4));
                        %if ~OG_ONLY
                            getstats(sim, OG_ONLY, FOCAL, EMPHASIS, TARGETS, ...
                                responses, RTs, act, acc, onsets, offsets, ...
                                is_target, correct, og_correct, is_inter_task, ...
                                true, true);
                        %end
                    end
                end % parfor subject_id
            
                
                
            end % for TARGETS
        end % for EMPHASIS 
    end % for FOCAL
end % for OG_ONLY
