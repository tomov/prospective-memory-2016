classdef Simulator < Model
    % Class for simulating trials with different inputs and outputs
    % and potentially different model parameters than the default ones
    %
    
    properties (Access = public)
        wm_capacity = 2;
        attention_factor = 0; % obsolete... for now
        next_available_hippo_id = 1;
        resting_wm;
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
        function self = Simulator(FOCAL, params, fitting_mode)
            self = self@Model(FOCAL, params);
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
        
        function threewayEM(self, stimulus, context, response)
            assert(self.next_available_hippo_id <= length(self.hippo_ids));
            stimulus_id = self.unit_id(stimulus);
            context_id = self.unit_id(context);
            response_id = self.unit_id(response);
            hippo_id = self.hippo_ids(self.next_available_hippo_id);
            self.weights(stimulus_id, hippo_id) = self.STIMULUS_TO_HIPPO;
            self.weights(context_id, hippo_id) = self.CONTEXT_TO_HIPPO;
            self.weights(hippo_id, response_id) = self.HIPPO_TO_TASK;
            self.next_available_hippo_id = self.next_available_hippo_id + 1;
        end
        
        function instruction(self, targets, include_in_WM)
            for target = targets
                target_unit = strcat('see', {' '}, target); % get perception unit name
                self.threewayEM(target_unit{1}, 'PM Context', 'PM Task');
                target_monitor_unit = strcat('Monitor ', {' '}, target);
                target_monitor_id = self.unit_id(target_monitor_unit{1});
                if include_in_WM
                    % make WM unit available
                    self.bias(target_monitor_id) = self.BIAS_FOR_ATTENTION;
                    % and give it some initial activation (according to params)
                    self.init_wm(self.wm_ids == target_monitor_id) = self.target_init;
                end
            end
        end
        
        function instruction_old(self, perceptions, targets, secs)
            from_ids = self.string_to_ids(perceptions);
            to_ids = self.string_to_ids(targets);
            from = zeros(1, self.N);
            from(from_ids) = self.MAXIMUM_ACTIVATION;
            to = zeros(1, self.N);
            to(to_ids) = self.MAXIMUM_ACTIVATION;
            duration = secs * self.CYCLES_PER_SEC;
            for cycle=1:duration
                % hebbian learning
                delta_w = self.LEARNING_RATE * from' * to;
                self.weights = self.weights + delta_w;
                
                % scale weights to fit model constraints
                sub = self.weights(self.perception_ids, self.task_ids);
                sub = sub * self.PERCEPTION_TO_TASK / max(sub(:));
                self.weights(self.perception_ids, self.task_ids) = sub;                
            end
            % add noise
            % ... TODO a little artificial at the end but whatever
            % also noise sigma is hardcoded and made up
            %self.weights(self.perception_ids, self.task_ids) = self.weights(self.perception_ids, self.task_ids) ...
            %    + normrnd(0, 3, size(self.perception_ids, 2), size(self.task_ids, 2));
        end
        
        % from http://grey.colorado.edu/CompCogNeuro/index.php/CCNBook/Networks/kWTA_Equations
        function kWTA_basic(self, k, ids)
            act = sort(self.net_input(ids), 'descend');
            if size(act, 2) <= k
                return
            end
            q = 0.5;
            threshold = act(k+1) + q*(act(k) - act(k+1));
            self.net_input(ids) = self.net_input(ids) - threshold;
        end

        function kWTA_average(self, k, ids)
            act = sort(self.net_input(ids), 'descend');
            if size(act, 2) <= k
                return
            end
            top = mean(act(1:k));
            bottom = mean(act(k+1:end));
            q = 0.5;
            threshold = bottom + q*(top - bottom);
            self.net_input(ids) = self.net_input(ids) - threshold;
        end
        
        % converting WM activations to activations for the feedforward
        % network
        function adapt_wm_act_to_ffwd_act(self)
            self.activation(:, self.wm_ids) = self.wm_act;
            %self.activation(self.wm_ids) = self.logistic(self.wm_act * 12 - 6);
        end
        
        function [responses, RTs, activation_log, accumulators_log, onsets, offsets, net_log] = run(self, stimuli, n_subjects)
            % initialize activations and outputs
            trial_duration = sum(cat(2, stimuli{:, 2})) * self.CYCLES_PER_SEC;
            self.accumulators = zeros(n_subjects, self.Nout);
            self.net_input = zeros(n_subjects, self.N);
            self.net_input_avg = zeros(n_subjects, self.N);
            activation_log = zeros(n_subjects, trial_duration, self.N);
            accumulators_log = zeros(n_subjects, trial_duration, self.Nout);
            net_log = zeros(n_subjects, trial_duration, self.N);
            self.activation = zeros(n_subjects, self.N);
            self.wm_act = zeros(n_subjects, size(self.wm_ids, 2));
            self.wm_net = zeros(n_subjects, size(self.wm_ids, 2));
            self.resting_wm = self.wm_act;
            responses = cell(n_subjects, size(stimuli, 1));
            RTs = zeros(n_subjects, size(stimuli, 1));
            onsets = zeros(n_subjects, size(stimuli, 1));
            offsets = zeros(n_subjects, size(stimuli, 1));
            cycles = 0;
            switched_to_PM_task = logical(zeros(n_subjects, 1)); % hacksauce
            switched_to_Inter_task = logical(zeros(n_subjects, 1));
            switched_to_OG_and_PM_from_Inter_task = logical(zeros(n_subjects, 1)); % hacksauce
            
            % for each input from the time series
            for ord=1:size(stimuli, 1)
                % get active input units for given stimulus
                % each stimulus string must be a comma-separated list of names of
                % input units
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
                        self.bias(target_monitor_id) = self.BIAS_FOR_ATTENTION;
                        % and give it some initial activation (according to params)
                        self.init_wm(self.wm_ids == target_monitor_id) = self.target_init;
                    end
                else
                    active_ids = self.string_to_ids(stimulus);
                end
                timeout = stimuli{ord, 2} * self.CYCLES_PER_SEC;
                
                % reset feedforward part of the network
                %self.activation(self.ffwd_ids) = 0;
                %self.net_input_avg(self.ffwd_ids) = -10;
                
                % reset response, output, and monitoring activations
                self.accumulators = zeros(n_subjects, size(self.output_ids, 2));
                
                % default output is timeout
                output_id = self.unit_id('timeout') * ones(n_subjects, 1);
                RTs(:, ord) = timeout;
                
                % simulate response to stimulus
                responded = logical(zeros(n_subjects, 1));
                settled = logical(zeros(n_subjects, 1));
                settle_cycle = zeros(n_subjects, 1);
                fprintf('ord = %d\n', ord);
                for cycle=1:timeout
                    %fprintf('   cycle = %d\n', cycle);
                    
                    % set input activations
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
                        self.wm_act(who_needs_wm_init, :) = ones(sum(who_needs_wm_init), 1) * self.init_wm;

                        % switch to inter task
                        % BIG TODO WTF -- see below

                        % switched to PM task
                        self.wm_act(switched_to_PM_task, self.wm_ids == self.unit_id('OG Task')) = self.init_wm(self.wm_ids == self.unit_id('OG Task'));
                        self.wm_act(switched_to_PM_task, self.wm_ids == self.unit_id('PM Task')) = self.init_wm(self.wm_ids == self.unit_id('PM Task'));

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
                    activation_log(:, cycles + cycle, :) = self.activation;
                    accumulators_log(:, cycles + cycle, :) = self.accumulators;
                    net_log(:, cycles + cycle, :) = self.net_input;
                   
                    % see if network has settled
                    if cycle > self.SETTLE_LEEWAY
                        from = cycles + cycle - self.SETTLE_LEEWAY + 1;
                        to = cycles + cycle - 1;
                        m = activation_log(:, from:to, :) - activation_log(:, from-1:to-1, :);
                        m = abs(mean(m, 3));
                        were_settled = settled;
                        settled = settled | (mean(m, 2) < self.SETTLE_MEAN_EPS & std(m, 0, 2) < self.SETTLE_STD_EPS);
                        just_settled = xor(were_settled, settled);
                        settle_cycle(just_settled) = cycle;
                        if ord == 1
                            self.resting_wm(just_settled, :) = self.wm_act(just_settled, :);
                        end
                        onsets(just_settled, ord) = cycles + cycle;
                    end

                    %{
                    if cycle > self.SETTLE_LEEWAY && ~is_settled
                        from = cycles + cycle - self.SETTLE_LEEWAY + 1;
                        to = cycles + cycle - 1;
                        m = activation_log(from:to, :) - activation_log(from-1:to-1, :);
                        m = abs(mean(m, 2));
                        %fprintf('%d -> %.6f, %6f\n', cycle, mean(m), std(m));
                        if mean(m) < self.SETTLE_MEAN_EPS && std(m) < self.SETTLE_STD_EPS
                            % it has settled
                            is_settled = true;
                            settle_cycle = cycle;
                            if ord == 1
                                self.resting_wm = self.wm_act;
                            end
                            % save stimulus onset
                            onsets = [onsets; cycles + cycle];
                        end
                    end
                    %}
                    
                    % calculate net inputs for all units
                    % WTF TODO ~responded insteaf of : ? faster?
                    %
                    self.net_input(~responded, self.ffwd_ids) = self.activation(~responded, :) * self.weights(:, self.ffwd_ids) ...
                        + ones(sum(~responded), 1) * self.bias(self.ffwd_ids);
                    self.net_input(~responded, self.wm_ids) = self.activation(~responded, self.ffwd_ids) * self.weights(self.ffwd_ids, self.wm_ids) ...
                        + self.wm_act(~responded, :) * self.weights(self.wm_ids, self.wm_ids) ...
                        + ones(sum(~responded), 1) * self.bias(self.wm_ids);
                    % unless we're fitting (i.e. when doing regular
                    % simulations), add noise to the net inputs
                    %
                    if ~self.fitting_mode
                        % TODO parametrize the noise
                        self.net_input(~responded, self.ffwd_ids) = self.net_input(~responded, self.ffwd_ids) + normrnd(0, 0.1, size(self.net_input(~responded, self.ffwd_ids)));
                        self.net_input(~responded, self.wm_ids) = self.net_input(~responded, self.wm_ids) + normrnd(0, 0.01, size(self.net_input(~responded, self.wm_ids)));
                    end
                    
                    % on instruction, oscillate around initial WM
                    % activations -- WTF bro
                    %if cycle < self.INSTRUCTION_CYLCES && (ord == 1 || switched_to_PM_task)
                    %    self.net_input(self.wm_ids) = self.net_input(self.wm_ids) + (self.init_wm - self.wm_act);
                    %    self.adapt_wm_act_to_ffwd_act();
                    %end
                    % for debugging...
                    %if cycles + cycle < 10
                    %    fprintf('%d: %.4f %.4f\n', cycles + cycle, self.activation(1), self.net_input(self.unit_id('PM Task')));
                    %end
                                                            
                    %self.bias(self.wm_ids) = self.bias(self.wm_ids) - 0.00001;
                    
                    % update activation levels for feedforward part of the
                    % network
                    self.net_input_avg(~responded, self.ffwd_ids) = self.TAU * self.net_input(~responded, self.ffwd_ids) + (1 - self.TAU) * self.net_input_avg(~responded, self.ffwd_ids);
                    self.activation(~responded, self.ffwd_ids) = self.logistic(self.net_input_avg(~responded, self.ffwd_ids)); % only update activations of those who haven't responded yet
                    
                    % same for WM module
                    % for WM module, activation f'n is linear and
                    % thresholded between 0 and 1
                    self.wm_act(~responded, :) = self.wm_act(~responded, :) + self.STEP_SIZE * self.net_input(~responded, self.wm_ids);
                    self.wm_act(~responded, :) = min(self.wm_act(~responded, :), self.MAXIMUM_ACTIVATION);
                    self.wm_act(~responded, :) = max(self.wm_act(~responded, :), self.MINIMUM_ACTIVATION);
                    self.adapt_wm_act_to_ffwd_act();
                   
                    % ---- update evidence accumulators (after network has settled) -----
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
                    act_to_subtract = act_max * ones(1, size(self.output_ids, 2));
                    % restore the actual max activations
                    act(max_linear_idx) = act_max;
                    % see evidence accumulation equations in paper
                    % basically, we subtract the max activation from everybody
                    % except from the max activation units themselves -- from them
                    % we subtract the second max activations
                    act_to_subtract(max_linear_idx) = second_act_max;
                    mu = self.EVIDENCE_ACCUM_ALPHA * (act - act_to_subtract);
                    % then we add noise proportional to that to the evidence accumulators
                    add = normrnd(mu, ones(size(mu)) * self.EVIDENCE_ACCUM_SIGMA);
                    self.accumulators(~responded & settled, :) = self.accumulators(~responded & settled, :) + add;

                    % check if activation threshold is met
                    %
                    [v, id] = max(self.accumulators, [], 2);
                    % get the output id for all subjects
                    all_output_id = self.output_ids(id);
                    % see which subjects passed the response threshold (and haven't responded yet)
                    who_just_responded = ~responded & settled & v > self.EVIDENCE_ACCUM_THRESHOLD;
                    % set their output ids
                    output_id(who_just_responded) = all_output_id(who_just_responded); % TODO might have to flip
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
                        break % no subjects left over
                    end

                    % single-subject version of the stuff above
                    % update evidence accumulators (after network has
                    % settled)
                    %{
                    if is_settled
                        act_sorted = sort(self.activation(self.output_ids), 'descend');
                        act_max = ones(n_subjects, size(self.output_ids, 2)) .* act_max;
                        act_max(self.activation(self.output_ids) == act_sorted(1)) = act_sorted(2);
                        mu = self.EVIDENCE_ACCUM_ALPHA * (self.activation(self.output_ids) - act_max);
                        add = normrnd(mu, ones(size(mu)) * self.EVIDENCE_ACCUM_SIGMA);
                        self.accumulators = self.accumulators + add;

                        % check if activation threshold is met
                        [v, id] = max(self.accumulators);
                        if ~responded && v > self.EVIDENCE_ACCUM_THRESHOLD
                            % save response and response time
                            output_id = self.output_ids(id);
                            RT = cycle - settle_cycle;
                            responded = true;
                            % a bit hacky, ALSO TODO does not work after
                            % timeout
                            break;
                        end
                    end
                    %}
                end
              
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
                
                %switched_to_PM_task = (self.activation(self.unit_id('PM Task')) > self.activation(self.unit_id('OG Task')));
                % TODO hacky...
                assert(size(self.resting_wm, 2) == size(self.wm_act, 2));
                %switched_to_PM_task = (self.wm_act(2) > self.wm_act(1) - 0.1);
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
