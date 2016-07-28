classdef Simulator < Model
    % Class for simulating trials with different inputs and outputs
    % and potentially different model parameters than the default ones
    %
    
    properties (Access = public)
        wm_capacity = 2;
        attention_factor = 0; % obsolete... for now
        next_available_hippo_id = 1;
        resting_wm;   % snapshot of wm activations at start of the experiment
        last_trial_wm; % snapshot of wm activations at start of previous trial
        net_input;
        net_input_avg;
        accumulators; % the evidence accumulators
        activation;   % activation of all nodes
        wm_act;
        wm_net;
        Nout;
        
        fitting_mode; % are we just fitting parameters? use more efficient setup
    end
    
    methods
        function self = Simulator(FOCAL, model_params, n_subjects, subject_params, fitting_mode)
            self = self@Model(FOCAL, model_params, n_subjects, subject_params);
            self.Nout = size(self.output_ids, 2);
            self.fitting_mode = fitting_mode;
        end
        
        function ids = string_to_ids(self, stimulus)
            units = strsplit(stimulus, ',');
            ids = [];
            for i=1:size(units, 2)
                active_unit = units{i};
                ids = [ids, self.unit_id(active_unit)];
            end
        end
        
        % from Jon's 1990 Stroop paper
        function act = logistic(self, net)
            act = 1 ./ (1 + exp(-net));
        end
        
        function threewayEM(self, stimulus, context, task, opposite_task)
            assert(self.next_available_hippo_id + 1 <= length(self.hippo_ids));
            % stimulus + context => task
            %
            stimulus_id = self.unit_id(stimulus);
            context_id = self.unit_id(context);
            task_id = self.unit_id(task);
            hippo_id = self.hippo_ids(self.next_available_hippo_id);
            self.weights(stimulus_id, hippo_id) = self.STIMULUS_TO_HIPPO;
            self.weights(context_id, hippo_id) = self.CONTEXT_TO_HIPPO;
            self.weights(hippo_id, task_id) = self.HIPPO_TO_TASK;
            self.next_available_hippo_id = self.next_available_hippo_id + 1;

            % context => opposite task (so it gets same amount of baseline excitation)
            % TODO HACK FIXME WTF
            %
            opposite_task_id = self.unit_id(opposite_task);
            hippo_id = self.hippo_ids(self.next_available_hippo_id);
            %self.weights(context_id, hippo_id) = self.CONTEXT_TO_HIPPO;
            %self.weights(hippo_id, opposite_task_id) = self.HIPPO_TO_TASK;
            self.next_available_hippo_id = self.next_available_hippo_id + 1;
        end
        
        function instruction(self, targets, include_in_WM)
            for target = targets
                target_unit = strcat('see', {' '}, target); % get perception unit name
                self.threewayEM(target_unit{1}, 'PM Context', 'PM Task', 'OG Task');
                target_monitor_unit = strcat('Monitor ', {' '}, target);
                target_monitor_id = self.unit_id(target_monitor_unit{1});
                if include_in_WM
                    % make WM unit available
                    self.bias(target_monitor_id) = self.BIAS_FOR_ATTENTION; % DEPRECATED -- see line below
                    self.biases(:, target_monitor_id) = self.subject_bias_for_attention;
                    % and give it some initial activation (according to params)
                    self.init_wm(:, self.wm_ids == target_monitor_id) = self.target_init;
                end
            end
        end
        
        % converting WM activations to activations for the feedforward
        % network
        function adapt_wm_act_to_ffwd_act(self)
            self.activation(:, self.wm_ids) = self.wm_act;
            %self.activation(self.wm_ids) = self.logistic(self.wm_act * 12 - 6);
        end
        
        function [responses, RTs, activation_log, accumulators_log, onsets, offsets, net_log] = run(self, stimuli)
            % initialize activations and outputs
            trial_duration = sum(cat(2, stimuli{:, 2})) * self.CYCLES_PER_SEC;
            self.accumulators = zeros(self.n_subjects, self.Nout);
            self.net_input = zeros(self.n_subjects, self.N);
            self.net_input_avg = zeros(self.n_subjects, self.N);
            activation_log = zeros(self.n_subjects, trial_duration, self.N);
            accumulators_log = zeros(self.n_subjects, trial_duration, self.Nout);
            net_log = zeros(self.n_subjects, trial_duration, self.N);
            self.activation = zeros(self.n_subjects, self.N);
            self.wm_act = zeros(self.n_subjects, size(self.wm_ids, 2));
            self.wm_net = zeros(self.n_subjects, size(self.wm_ids, 2));
            self.resting_wm = self.wm_act;
            self.last_trial_wm = self.wm_act;
            responses = cell(self.n_subjects, size(stimuli, 1));
            RTs = zeros(self.n_subjects, size(stimuli, 1));
            onsets = zeros(self.n_subjects, size(stimuli, 1));
            offsets = zeros(self.n_subjects, size(stimuli, 1));
            cycles = 0;
            switched_to_PM_task = logical(zeros(self.n_subjects, 1)); % hacksauce
            switched_to_Inter_task = logical(zeros(self.n_subjects, 1));
            switched_to_OG_and_PM_from_Inter_task = logical(zeros(self.n_subjects, 1)); % hacksauce
            
            %
            % for each input from the time series (i.e. for each trial)
            %
            for ord=1:size(stimuli, 1)
                % get active input units for given stimulus
                % each stimulus string must be a comma-separated list of names of
                % input units
                %
                stimulus = stimuli{ord, 1};
                if strcmp(stimulus, 'switch to Inter Task')
                    active_ids = [];
                    switched_to_Inter_task(:) = true;
                elseif strcmp(stimulus, 'switch back to OG and PM')
                    active_ids = [];
                    switched_to_OG_and_PM_from_Inter_task(:) = true;
                    for target = {'tortoise'} % TODO dedupe with instruction() % TODO pass target as param
                        target_monitor_unit = strcat('Monitor ', {' '}, target);
                        target_monitor_id = self.unit_id(target_monitor_unit{1});
                        % make WM unit available
                        self.bias(target_monitor_id) = self.BIAS_FOR_ATTENTION; % DEPRECATED -- see line below
                        self.biases(:, target_monitor_id) = self.subject_bias_for_attention;
                        % and give it some initial activation (according to params)
                        self.init_wm(:, self.wm_ids == target_monitor_id) = self.target_init;
                    end
                else
                    active_ids = self.string_to_ids(stimulus);
                end
                timeout = stimuli{ord, 2} * self.CYCLES_PER_SEC;
                
                % reset feedforward part of the network
                % but not today...
                %self.activation(self.ffwd_and_em_ids) = 0;
                %self.net_input_avg(self.ffwd_and_em_ids) = -10;
                
                % reset response, output, and monitoring activations
                self.accumulators = zeros(self.n_subjects, size(self.output_ids, 2));
                
                %
                % simulate response to given stimulus i.e. a trial
                % for all subjects simultaneously (woohoo)
                %
                responded = logical(zeros(self.n_subjects, 1));
                settled = logical(zeros(self.n_subjects, 1));
                settle_cycle = zeros(self.n_subjects, 1);
                %fprintf('ord = %d\n', ord);
                for cycle=1:timeout
                    %fprintf('   cycle = %d\n', cycle);
                    
                    % set input activations
                    %
                    self.activation(:, self.input_ids) = 0;                    
                    self.activation(settled, active_ids) = self.INPUT_ACTIVATION;
                    
                    % hack for testing different activations
                    %self.wm_act = self.init_wm;
                    %self.activation(self.wm_ids) = (self.wm_act + 5) / 10;%self.logistic(self.wm_act);
                   
                    % initialize WM at beginning of block (i.e. first
                    % trial),
                    % or after a PM switch
                    if cycle < self.INSTRUCTION_CYLCES
                        % initial WM
                        who_needs_wm_init = ord == 1 | switched_to_OG_and_PM_from_Inter_task;
                        self.wm_act(who_needs_wm_init, :) = self.init_wm(who_needs_wm_init, :);

                        % switch to inter task
                        % BIG TODO WTF -- see below

                        % switched to PM task
                        self.wm_act(switched_to_PM_task, self.wm_ids == self.unit_id('OG Task')) = self.last_trial_wm(switched_to_PM_task, self.wm_ids == self.unit_id('OG Task'));
                        self.wm_act(switched_to_PM_task, self.wm_ids == self.unit_id('PM Task')) = self.last_trial_wm(switched_to_PM_task, self.wm_ids == self.unit_id('PM Task'));

                        self.adapt_wm_act_to_ffwd_act();
                    end

                    %{
                    if cycle < self.INSTRUCTION_CYLCES && (ord == 1 || switched_to_PM_task || switched_to_Inter_task || switched_to_OG_and_PM_from_Inter_task)
                        if ord == 1 || switched_to_OG_and_PM_from_Inter_task
                            self.wm_act = self.init_wm;
                        elseif switched_to_Inter_task
                            % reset all
                            % TODO parametrize
                            for target = {'tortoise'} % TODO dedupe with instruction() % TODO pass target as param
                                target_monitor_unit = strcat('Monitor ', {' '}, target);
                                target_monitor_id = self.unit_id(target_monitor_unit{1});
                                % make WM unit available
                                self.bias(target_monitor_id) = self.BIAS_FOR_ATTENTION;
                                % and give it some initial activation (according to params)
                                self.init_wm(self.wm_ids == target_monitor_id) = 0;
                            end
                            self.wm_act(self.wm_ids == self.unit_id('OG features')) = self.init_wm(self.wm_ids == self.unit_id('OG features'));
                            self.wm_act(self.wm_ids == self.unit_id('OG Task')) = 1;
                            self.wm_act(self.wm_ids == self.unit_id('PM Task')) = 0;
                            self.wm_act(self.wm_ids == self.unit_id('PM Context')) = 0.2;
                            self.wm_act(self.wm_ids == self.unit_id('Other Context')) = 1;
                        else
                            assert(switched_to_PM_task);
                            % only reset tasks
                            self.wm_act(self.wm_ids == self.unit_id('OG Task')) = self.init_wm(self.wm_ids == self.unit_id('OG Task'));
                            self.wm_act(self.wm_ids == self.unit_id('PM Task')) = self.init_wm(self.wm_ids == self.unit_id('PM Task'));
                        end
                        self.adapt_wm_act_to_ffwd_act();
                    end
                    %}

                    % log activation for plotting
                    %
                    activation_log(:, cycles + cycle, :) = self.activation;
                    accumulators_log(:, cycles + cycle, :) = self.accumulators;
                    net_log(:, cycles + cycle, :) = self.net_input;
                   
                    % see if network has settled
                    %
                    if cycle > self.SETTLE_LEEWAY
                        from = cycles + cycle - self.SETTLE_LEEWAY + 1;
                        to = cycles + cycle - 1;
                        m = activation_log(:, from:to, :) - activation_log(:, from-1:to-1, :);
                        m = abs(mean(m, 3));
                        were_settled = settled;
                        settled = settled | (mean(m, 2) < self.SETTLE_MEAN_EPS & std(m, 0, 2) < self.SETTLE_STD_EPS);
                        just_settled = xor(were_settled, settled); % only update subjects whose activations have settled
                        settle_cycle(just_settled) = cycle;
                        if ord == 1
                            self.resting_wm(just_settled, :) = self.wm_act(just_settled, :);
                        end
                        self.last_trial_wm(just_settled, :) = self.wm_act(just_settled, :);
                        onsets(just_settled, ord) = cycles + cycle;
                    end
                    
                    % calculate net inputs for all units
                    %
                    self.net_input(~responded, self.ffwd_and_em_ids) = self.activation(~responded, :) * self.weights(:, self.ffwd_and_em_ids) ...
                        + self.biases(~responded, self.ffwd_and_em_ids);
                        %+ repmat(self.bias(self.ffwd_and_em_ids), sum(~responded), 1);
                    self.net_input(~responded, self.wm_ids) = self.activation(~responded, self.ffwd_and_em_ids) * self.weights(self.ffwd_and_em_ids, self.wm_ids) ...
                        + self.wm_act(~responded, :) * self.weights(self.wm_ids, self.wm_ids) ...
                        + self.biases(~responded, self.wm_ids);
                        %+ repmat(self.bias(self.wm_ids), sum(~responded), 1);
                    % unless we're fitting (i.e. when doing regular
                    % simulations), add noise to the net inputs
                    %
                    if ~self.fitting_mode
                        self.net_input(~responded, self.ffwd_and_em_ids) = self.net_input(~responded, self.ffwd_and_em_ids) + normrnd(0, self.NOISE_SIGMA_FFWD, size(self.net_input(~responded, self.ffwd_and_em_ids)));
                        self.net_input(~responded, self.wm_ids) = self.net_input(~responded, self.wm_ids) + normrnd(0, self.NOISE_SIGMA_WM, size(self.net_input(~responded, self.wm_ids)));
                    end
                                        
                    % update activation levels for feedforward part of the network
                    % only update activations of those who haven't responded yet
                    self.net_input_avg(~responded, self.ffwd_and_em_ids) = self.TAU * self.net_input(~responded, self.ffwd_and_em_ids) + (1 - self.TAU) * self.net_input_avg(~responded, self.ffwd_and_em_ids);
                    self.activation(~responded, self.ffwd_ids) = self.logistic(self.net_input_avg(~responded, self.ffwd_ids));
                    self.activation(~responded, self.em_ids) = self.logistic(self.EM_GAIN * self.net_input_avg(~responded, self.em_ids)); % the sigmoid gain makes EM activation more step-like (all-or-nothing)
                    em = self.activation(~responded, self.em_ids); % TODO HACK FIXME --> artifically thresholding the sigmoid
                    em(em < 0.01) = 0; % TODO PARAM
                    self.activation(~responded, self.em_ids) = em;
                    
                    % same for WM module
                    % for WM module, activation f'n is linear and
                    % thresholded between 0 and 1
                    self.wm_act(~responded, :) = self.wm_act(~responded, :) + self.STEP_SIZE * self.net_input(~responded, self.wm_ids);
                    self.wm_act(~responded, :) = min(self.wm_act(~responded, :), self.MAXIMUM_ACTIVATION);
                    self.wm_act(~responded, :) = max(self.wm_act(~responded, :), self.MINIMUM_ACTIVATION);
                    self.adapt_wm_act_to_ffwd_act();
                   
                    %
                    % update evidence accumulators (after network has settled)
                    %
                    
                    % take only the output activation units we care about
                    act = self.activation(~responded & settled, self.output_ids);
                    % find output unit with max activation for each subject
                    [act_max, max_idx] = max(act, [], 2);
                    % get their indices
                    max_linear_idx = sub2ind(size(act), 1:length(max_idx), max_idx');
                    % set them temporarily to -inf so we can find the second max activation units
                    act(max_linear_idx) = -inf;
                    second_act_max = max(act, [], 2);
                    % scale act_max to same size as the output units
                    act_to_subtract = repmat(act_max, 1, size(self.output_ids, 2));
                    % restore the actual max activations
                    act(max_linear_idx) = act_max;
                    % see evidence accumulation equations in paper
                    % basically, we subtract the max activation from everybody
                    % except from the max activation units themselves -- from them
                    % we subtract the second max activations
                    act_to_subtract(max_linear_idx) = second_act_max;
                    mu = self.EVIDENCE_ACCUM_ALPHA * (act - act_to_subtract);
                    % then we add noise proportional to that to the evidence accumulators
                    add = normrnd(mu, repmat(self.EVIDENCE_ACCUM_SIGMA, size(mu)));
                    self.accumulators(~responded & settled, :) = self.accumulators(~responded & settled, :) + add;

                    % check if response threshold is met
                    %
                    [v, id] = max(self.accumulators, [], 2);
                    % get the output id for all subjects
                    all_output_id = self.output_ids(id);
                    % see which subjects passed the response threshold (and haven't responded yet)
                    who_just_responded = ~responded & settled & v > self.EVIDENCE_ACCUM_THRESHOLD;
                    % set their output ids
                    output_id(who_just_responded) = all_output_id(who_just_responded);
                    % record the actual response
                    responses(who_just_responded & (switched_to_Inter_task | switched_to_OG_and_PM_from_Inter_task), ord) = {'Switch'};
                    who_just_responded_and_is_not_a_switch = who_just_responded & (~switched_to_Inter_task & ~switched_to_OG_and_PM_from_Inter_task);
                    responses(who_just_responded_and_is_not_a_switch, ord) = self.units(all_output_id(who_just_responded_and_is_not_a_switch));
                    % set their response times
                    RTs(who_just_responded, ord) = cycle - settle_cycle(who_just_responded);
                    % mark that they responded
                    responded(who_just_responded) = true;
                    % and do some bookkeeping
                    offsets(who_just_responded, ord) = cycles + cycle;
                    
                    if sum(~responded) == 0
                        break % everyone has responded
                    end
                end
              
                % whoever didn't respond is a timeout
                responses(~responded, ord) = {'timeout'};
                RTs(~responded, ord) = timeout;

                %{
                % record response and response time
                if switched_to_Inter_task || switched_to_OG_and_PM_from_Inter_task
                    output = 'Switch';
                else
                    output = self.units{output_id};
                end
                offsets = [offsets; cycles + cycle];
                responses = [responses; {output}];
                RTs = [RTs; RT];
                %}
                cycles = cycles + cycle;                
                
                %switched_to_PM_task = (self.wm_act(2) > self.wm_act(1) -
                %0.1); <--- DOESN'T quite work (partial switch)
                switched_to_PM_task = (self.wm_act(:, 2) > self.resting_wm(:, 2) + 0.01);
                switched_to_Inter_task(:) = false;
                switched_to_OG_and_PM_from_Inter_task(:) = false;
            end
            
            fprintf('Total simulation cycles = %d', cycles');
            
            activation_log(:,cycles:end,:) = [];
            accumulators_log(:,cycles:end,:) = [];
            net_log(:,cycles:end,:) = [];
        end
    end
end
