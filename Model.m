classdef Model < handle
    % All constants and predefined variables in the model
    %
    
    properties (Constant = true)
        % PDP model parameters
        
        NOISE_SIGMA = 0.1; % TODO -- ??
        STEP_SIZE = 0.05;
        DECAY = 0.01;
        CYCLES_PER_SEC = 500; % 500 ; CANNOT be lower, e.g. 200...
        SETTLE_MEAN_EPS = 1e-3; % adjust these when you add noise to the model
        SETTLE_STD_EPS = 1e-4; % ...this too
        TAU = 0.1; % rate constant from Jon's paper = 0.1
        INSTRUCTION_CYLCES = 2/Model.TAU;
        RESET_CYCLES = Model.INSTRUCTION_CYLCES;
        SETTLE_LEEWAY = 2*Model.INSTRUCTION_CYLCES;
        EVIDENCE_ACCUM_SIGMA = 0.05;
        EVIDENCE_ACCUM_ALPHA = 0.05;
        EVIDENCE_ACCUM_THRESHOLD = 1;
        
        % activation levels

        MAXIMUM_ACTIVATION = 1;
        MINIMUM_ACTIVATION = 0;
        
        INPUT_ACTIVATION = 1;
        
        LAMBDA = 0; % regularization constant for backprop TODO set to 1 or 0.1 or something
    end

    
    properties (Access = public)
        

        % --- begin connection weights ---
        
        % perception
        
        BIAS_FOR_PERCEPTION = -18;
        PERCEPTION_INHIBITION = 0;
        
        INPUT_TO_PERCEPTION = 15;
        INPUT_TO_PERCEPTION_INHIBITION = 0;
        
        ATTENTION_TO_PERCEPTION = 8;
        ATTENTION_TO_PERCEPTION_INHIBITION = 0;

        % responses
        
        BIAS_FOR_RESPONSES = -12;
        RESPONSE_INHIBITION = -5;
        
        PERCEPTION_TO_RESPONSE = 5;
        PERCEPTION_TO_RESPONSE_INHIBITION = 0;

        TASK_TO_RESPONSE = 8;
        TASK_TO_RESPONSE_INHIBITION = 0;
        
        % outputs
        
        BIAS_FOR_OUTPUTS = 0;
        OUTPUT_INHIBITION = -3;
        
        RESPONSE_TO_OUTPUT = 1;
        RESPONSE_TO_OUTPUT_INHIBITION = 0;
        
        % WM and EM units follow
        %
        
        BIAS_WHEN_OFF = -100;
                
        % task representation (WM)
        
        BIAS_FOR_TASK = 3;
        TASK_INHIBITION = -2;
        TASK_SELF = -2 + 0.0004;
        
        ATTENTION_TO_TASK = -1;
        
        HIPPO_TO_TASK = 10;
        %PERCEPTION_TO_TASK = 1.2;  % EM = speed of task switch --
        %DEPRECATEd; see hippo
        
        % feature attention (WM)
        
        BIAS_FOR_ATTENTION = 3;
        ATTENTION_INHIBITION = -2;
        ATTENTION_SELF = -2 + 0.0004;
        
        TASK_TO_ATTENTION = -1;
        
        % context units
        
        BIAS_FOR_CONTEXT = 3;
        CONTEXT_INHIBITION = -2;
        CONTEXT_SELF = -2 + 0.0004;
        
        TASK_TO_CONTEXT = -1;
        ATTENTION_TO_CONTEXT = -1;
        CONTEXT_TO_TASK_AND_ATTENTION = -1;
        
        % hippocampus
        
        % FUTURE -20
        BIAS_FOR_HIPPO = -23; %-32;  % must be < -10, o/w tasks drift b/c of (super small) input current from hippo
        
        % FUTURE 20, 12
        STIMULUS_TO_HIPPO = 16; % 30
        CONTEXT_TO_HIPPO = 16;  % 20
        
        %OUTPUT_TO_SELF = 0; % makes response->output more like copying rather than integration
        %RESPONSE_TO_SELF = 0;
        
        % --- end of connection weights ---
        
        % EM parameters
        
        LEARNING_RATE = 0.01;
                
        % variables
        
        units
        N
        unit_id
        
        input_units
        perception_units
        response_units
        output_units
        task_units
        attention_units
        context_units
        hippo_units
        
        input_ids      % layer 0 -- input layer: the raw stimulus input
        perception_ids % layer 1 -- feature perception layer: the perceived input = raw input + feature attention
        response_ids   % layer 2 -- associative layer: the response = perceived input + task attention
        output_ids     % layer 3 -- output layer: the motor action
        task_ids       % task representation WM units
        attention_ids  % feature attention WM units
        context_ids    % context WM units
        hippo_ids      % hippocampus units encoding three-way associations from features (layer 1) and context units to task units
        
        wm_ids         % units in the WM part of the model (the top-down attentional modulation)
        ffwd_ids       % units in the feed-forward part of the model (the bottom-up stimulus-response pathway)
        
        connections
        weights        % the weights of the network
        bias
        init_wm
        target_init
    end
    
    methods
        function lateral_inhibition(self, units, weight)
            for i=1:size(units, 2)
                for j=1:size(units, 2)
                    if i ~= j
                        if ~ismember([units(i), units(j)], self.connections(:,1:2), 'rows')
                            self.connections = [self.connections;
                                units(i), units(j), weight];
                            %fprintf('%s -> %s: %d (LI)\n', self.units{units(i)}, self.units{units(j)}, weight);
                        end
                    end
                end
            end
        end
        
        
        function forward_parallel(self, from, to, weight)
            assert(size(from, 2) <= size(to, 2));
            for i=1:size(from, 2)
                self.connections = [self.connections;
                    from(i), to(i), weight];
                %fprintf('%s -> %s: %d\n', self.units{from(i)}, self.units{to(i)}, weight);
            end
        end

        function forward_all_to_all(self, from, to, weight)
            for i=1:size(from, 2)
                for j=1:size(to, 2)
                    if ~ismember([from(i), to(j)], self.connections(:,1:2), 'rows')
                        self.connections = [self.connections;
                            from(i), to(j), weight];
                        %fprintf('%s -> %s: %d (FI)\n', self.units{from(i)}, self.units{to(j)}, weight);
                    end
                end
            end
        end

        function self_excitation(self, units, weight)
            for i=1:size(units, 2)
                self.connections = [self.connections;
                    units(i), units(i), weight];
                %fprintf('%s -> %s: %d (SE)\n', self.units{units(i)}, self.units{units(i)}, weight);
            end
        end
        
        function self = Model(FOCAL, params, have_third_task)
            % specify unit names in each layer
            words = {
                'tortoise', 'physics', 'crocodile', 'math', ... % words
                'dog', 'cat', 'panda', 'kiwi', 'monkey', ... % additional words for Exp 3/4
                };
            categories = {
                'a subject', 'an animal', ... % categories
                };
            syllables = {
                'tor'
                };
            self.input_units = [words categories];
            % PERCEPTION MUST BE IN SAME ORDER AS INPUT UNITS
            self.perception_units = strcat('see', {' '}, self.input_units')';
            self.perception_units = [self.perception_units, strcat('see', {' '}, syllables')';];
            self.response_units = {
                'A Subject', 'An Animal', 'No Match 1', 'No Match 2', ... % OG task
                'PM Response', ...   % PM task
                'Domestic', 'Wild'   % Interleaved task
                };
            self.output_units = {
                'Yes', 'No', 'PM'
                };
            self.task_units = {
                'OG Task', 'PM Task', 'Inter Task'
                };
            self.context_units = {
                'PM Context', 'Other Context'
                };
            self.attention_units = {
                'OG features'
                };
            self.attention_units = [self.attention_units, strcat('Monitor ', {' '}, words')'];
            self.attention_units = [self.attention_units, strcat('Monitor ', {' '}, syllables')'];
            self.hippo_units = {
                'hippo 1', ...
                'hippo 2', ...
                'hippo 3', ...
                'hippo 4', ...
                'hippo 5', ...
                'hippo 6', ...
                'hippo 7', ...
                'hippo 8'
                };
            self.units = [
                self.input_units, ...
                self.perception_units, ...
                self.response_units, ...
                self.output_units, ...
                self.task_units, ...
                self.attention_units, ...
                self.context_units, ...
                self.hippo_units, ...
                {'timeout'}
                ];
            
            % generate indices (for convenience)
            self.N = size(self.units, 2);
            self.unit_id = containers.Map(self.units, 1:self.N);

            self.input_ids = cellfun(@self.unit_id, self.input_units);
            self.perception_ids = cellfun(@self.unit_id, self.perception_units);
            self.response_ids = cellfun(@self.unit_id, self.response_units);
            self.output_ids = cellfun(@self.unit_id, self.output_units);
            self.task_ids = cellfun(@self.unit_id, self.task_units);
            self.attention_ids = cellfun(@self.unit_id, self.attention_units);
            self.context_ids = cellfun(@self.unit_id, self.context_units);
            self.hippo_ids = cellfun(@self.unit_id, self.hippo_units);
            
            self.ffwd_ids = [
                self.input_ids ...
                self.perception_ids ...
                self.response_ids ...
                self.output_ids ...
                self.hippo_ids];
            self.wm_ids = [self.task_ids self.attention_ids self.context_ids];

            % initialize free parameters (based on PM instruction, task, etc)
            self.init_wm = zeros(1, length(self.wm_ids));
            self.init_wm(self.wm_ids == self.unit_id('OG Task')) = params(1);
            self.init_wm(self.wm_ids == self.unit_id('PM Task')) = params(2);
            self.init_wm(self.wm_ids == self.unit_id('OG features')) = params(3);
            self.target_init = params(4);
            self.BIAS_FOR_TASK = params(5);
            self.BIAS_FOR_ATTENTION = params(6);
            self.init_wm(self.wm_ids == self.unit_id('PM Context')) = params(7);
            self.init_wm(self.wm_ids == self.unit_id('Other Context')) = params(8);
            self.BIAS_FOR_CONTEXT = params(9);

            % ---==== specify connections between units ====---
            
            self.connections = [
                % task monitoring to responses
                self.unit_id('OG Task')        , self.unit_id('A Subject')         , self.TASK_TO_RESPONSE;
                self.unit_id('OG Task')        , self.unit_id('An Animal')         , self.TASK_TO_RESPONSE;
                self.unit_id('OG Task')        , self.unit_id('No Match 1')        , self.TASK_TO_RESPONSE;
                self.unit_id('OG Task')        , self.unit_id('No Match 2')        , self.TASK_TO_RESPONSE;
                self.unit_id('PM Task')        , self.unit_id('PM Response')       , self.TASK_TO_RESPONSE;
                self.unit_id('Inter Task')     , self.unit_id('Domestic')          , self.TASK_TO_RESPONSE;
                self.unit_id('Inter Task')     , self.unit_id('Wild')              , self.TASK_TO_RESPONSE;
                
                % perception to response mapping (direct OG pathway)
                %
                
                % -- categories to categories
                self.unit_id('see a subject')                  , self.unit_id('A Subject')          , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see an animal')                  , self.unit_id('An Animal')          , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see a subject')                  , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see an animal')                  , self.unit_id('No Match 2')         , self.PERCEPTION_TO_RESPONSE;

                
                % -- animals to matching categories
                self.unit_id('see physics')                , self.unit_id('A Subject')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see math')                   , self.unit_id('A Subject')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see tortoise')               , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see crocodile')              , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see dog')                    , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see cat')                    , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see panda')                  , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see kiwi')                   , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see monkey')                 , self.unit_id('An Animal')         , self.PERCEPTION_TO_RESPONSE;
                
                % -- default response is No Match
                self.unit_id('see physics')                , self.unit_id('No Match 2')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see math')                   , self.unit_id('No Match 2')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see tortoise')               , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see crocodile')              , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see dog')                    , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see cat')                    , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see panda')                  , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see kiwi')                   , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                self.unit_id('see monkey')                 , self.unit_id('No Match 1')         , self.PERCEPTION_TO_RESPONSE;
                
                % perception to response mapping (Interleaved task)
                % notice all is *2 b/c there's only 1 input per response =>
                % weaker
                %

                self.unit_id('see tortoise')               , self.unit_id('Wild')         , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see crocodile')              , self.unit_id('Wild')         , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see dog')                    , self.unit_id('Domestic')     , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see cat')                    , self.unit_id('Domestic')     , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see panda')                  , self.unit_id('Wild')         , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see kiwi')                   , self.unit_id('Wild')         , self.PERCEPTION_TO_RESPONSE * 2;
                self.unit_id('see monkey')                 , self.unit_id('Wild')         , self.PERCEPTION_TO_RESPONSE * 2;
                
                % raw inputs to perception -- nonfocal PM targets
                self.unit_id('tortoise')               , self.unit_id('see tor')         , self.INPUT_TO_PERCEPTION;
                
                % attention to perception -- nonfocal PM targets
                self.unit_id('Monitor tor')            , self.unit_id('see tor')              , self.ATTENTION_TO_PERCEPTION;
                
                % responses to outputs                
                self.unit_id('A Subject')           , self.unit_id('Yes')            , self.RESPONSE_TO_OUTPUT;
                self.unit_id('An Animal')           , self.unit_id('Yes')            , self.RESPONSE_TO_OUTPUT;
                self.unit_id('No Match 1')          , self.unit_id('No')             , self.RESPONSE_TO_OUTPUT;
                self.unit_id('No Match 2')          , self.unit_id('No')             , self.RESPONSE_TO_OUTPUT;
                self.unit_id('PM Response')         , self.unit_id('PM')             , self.RESPONSE_TO_OUTPUT;
                self.unit_id('Domestic')            , self.unit_id('Yes')              , self.RESPONSE_TO_OUTPUT;
                self.unit_id('Wild')                , self.unit_id('No')              , self.RESPONSE_TO_OUTPUT;
            ];
        
        
            % -- TODO HACK to make the PM task work, you need to put
            % it up to baseline (the winning OG response gets x2 inputs)
            %if FOCAL
            % FUTURE
                self.connections = [self.connections;
                    self.unit_id('see a subject')              , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see an animal')              , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see physics')                , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see math')                   , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see tortoise')               , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see crocodile')              , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see dog')                    , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see cat')                    , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see panda')                  , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see kiwi')                   , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                    self.unit_id('see monkey')                 , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                ];
            %else
            %FUTURE
            
                self.connections = [self.connections;
                   % self.unit_id('see tor')                    , self.unit_id('PM Response')         , self.PERCEPTION_TO_RESPONSE;
                ];
            
            %end
        
            % WM LCA vertical mutual inhibition
            self.forward_all_to_all(self.attention_ids, self.task_ids, self.ATTENTION_TO_TASK);
            self.forward_all_to_all(self.task_ids, self.attention_ids, self.TASK_TO_ATTENTION);
            
            self.forward_all_to_all(self.attention_ids, self.context_ids, self.ATTENTION_TO_CONTEXT);
            self.forward_all_to_all(self.task_ids, self.context_ids, self.TASK_TO_CONTEXT);
            self.forward_all_to_all(self.context_ids, [self.task_ids self.attention_ids], self.CONTEXT_TO_TASK_AND_ATTENTION);
            
            % perception to task representation (indirect PM pathway)
            self.forward_all_to_all(self.perception_ids, self.task_ids, 0); % EM!!!
            
            % OG attention to perception
            from = self.unit_id('OG features');
            to = cellfun(@self.unit_id, strcat('see', {' '}, [words categories]')');
            self.forward_all_to_all(from, to, self.ATTENTION_TO_PERCEPTION);
            
            % attention (target monitoring) to perception
            from = cellfun(@self.unit_id, strcat('Monitor', {' '}, words')');
            to = cellfun(@self.unit_id, strcat('see', {' '}, words')');
            self.forward_parallel(from, to, self.ATTENTION_TO_PERCEPTION);

            % raw inputs to perception (cont'd)
            self.forward_parallel(self.input_ids, self.perception_ids, self.INPUT_TO_PERCEPTION);
            
            % forward inhibitions
            self.forward_all_to_all(self.input_ids, self.perception_ids, self.INPUT_TO_PERCEPTION_INHIBITION);
            self.forward_all_to_all(self.perception_ids, self.response_ids, self.PERCEPTION_TO_RESPONSE_INHIBITION);
            self.forward_all_to_all(self.task_ids, self.response_ids, self.TASK_TO_RESPONSE_INHIBITION);
            self.forward_all_to_all(self.response_ids, self.output_ids, self.RESPONSE_TO_OUTPUT_INHIBITION);

            % lateral inhibitions
            self.lateral_inhibition(self.perception_ids, self.PERCEPTION_INHIBITION);
            self.lateral_inhibition(self.response_ids, self.RESPONSE_INHIBITION);
            self.lateral_inhibition(self.output_ids, self.OUTPUT_INHIBITION);
            self.lateral_inhibition(self.task_ids, self.TASK_INHIBITION);
            self.lateral_inhibition(self.attention_ids, self.ATTENTION_INHIBITION);
            self.lateral_inhibition(self.context_ids, self.CONTEXT_INHIBITION);

            % self excitations
            self.self_excitation(self.task_ids, self.TASK_SELF);
            self.self_excitation(self.attention_ids, self.ATTENTION_SELF);
            self.self_excitation(self.context_ids, self.CONTEXT_SELF);
            
            % generate weight matrix from defined connections
            self.weights = sparse(self.connections(:,1), self.connections(:,2), self.connections(:,3), ...
                self.N, self.N);

            self.bias = zeros(1, self.N);
            self.bias(self.perception_ids) = self.BIAS_FOR_PERCEPTION;
            self.bias(self.response_ids) = self.BIAS_FOR_RESPONSES;
            self.bias(self.output_ids) = self.BIAS_FOR_OUTPUTS;
            self.bias(self.task_ids) = self.BIAS_FOR_TASK;
            self.bias(self.attention_ids) = self.BIAS_WHEN_OFF; % all attention units don't exist by default
                self.bias(self.unit_id('OG features')) = self.BIAS_FOR_ATTENTION; % ...except for OG features
            self.bias(self.context_ids) = self.BIAS_FOR_CONTEXT;
            self.bias(self.hippo_ids) = self.BIAS_FOR_HIPPO;            
            if ~have_third_task
                self.bias(self.unit_id('Inter Task')) = self.BIAS_WHEN_OFF;
            end
        end
        
        function EM = print_EM(self)
            EM = full(self.weights(self.perception_ids, self.task_ids));
            EM = EM';
        end
        
        function show_EM(self)
            EM = self.print_EM();
            m = max(EM(:));
            figure;
            imshow(imresize(1-EM/m, 60, 'box'))
        end
    end
end
