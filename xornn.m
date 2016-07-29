clear ; close all; clc

%load('what-the-fuck.mat');

NOISE_SIGMA = 0.1; % TODO -- ??
STEP_SIZE = 0.05;
DECAY = 0.01;
CYCLES_PER_SEC = 500; % 500 ; CANNOT be lower, e.g. 200...
SETTLE_MEAN_EPS = 1e-3; % adjust these when you add noise to the model
SETTLE_STD_EPS = 1e-4; % 1e-4; % ...this too
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

LAMBDA = 0; %0.001; % regularization constant for backprop TODO set to 1 or 0.1 or something

Nout = 3;
N = 23;

% the actual values

%X = rand(5000, 2) > 0.5;
%y = xor(X(:, 1), X(:, 2)) + 1;

stimuli = X;
correct = zeros(size(X, 1), length(output_ids));
for ord = 1:size(X, 1)
    correct(ord, y(ord)) = 1;
end

input_ids = 1:12;
hidden_ids = 13:20;
output_ids = 21:23;

% training stuff
%


weights = zeros(N, N);
bias = zeros(1, N);

del = zeros(1, N);
weightsGradient = zeros(size(weights));
biasGradient = zeros(size(bias));

epsilon_init = 0.12; % from ML course TODO move to self.

% initialize weights & biases

weights(hidden_ids, output_ids) = rand(length(hidden_ids), length(output_ids)) * 2 * epsilon_init - epsilon_init;
weights(input_ids, hidden_ids) = rand(length(input_ids), length(hidden_ids)) * 2 * epsilon_init - epsilon_init;

bias(output_ids) = rand(1, length(output_ids)) * 2 * epsilon_init - epsilon_init;
bias(hidden_ids) = rand(1, length(hidden_ids)) * 2 * epsilon_init - epsilon_init;

m = length(X);
expected = zeros(m, length(output_ids));
idx = sub2ind(size(expected), 1:size(expected, 1), y');
expected(idx) = 1;


% -------- TRAIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

initial_nn_params = [weights(:); bias(:)];

options = optimset('MaxIter', 1000);

costFn = @(p) costFunction(p, ...
                           N, ...
                           input_ids, ...
                           hidden_ids, ...
                           output_ids, ...
                           X, expected, LAMBDA);

[nn_params, cost] = fmincg(costFn, initial_nn_params, options);

weights = reshape(nn_params(1:N*N), N, N);
bias = reshape(nn_params(N*N+1:end), 1, N);

% --- end of training.................
            
for shit = 3:3

            experiment_duration = length(stimuli) * CYCLES_PER_SEC;
            net_input = zeros(1, N);
            net_input_avg = zeros(1, N);
            activation_log = zeros(experiment_duration, N);
            accumulators_log = zeros(experiment_duration, Nout);
            net_log = zeros(experiment_duration, N);
            activation = zeros(1, N);
            net_input_avg = zeros(1, N);
            responses = [];
            RTs = [];
            onsets = [];
            offsets = [];
            cycles = 0;
            switched_to_PM_task = false;
            
            
            hits = 0;
            misses = 0;
            is_hit = zeros(size(stimuli, 1), 1);
            
            last_trial = size(stimuli, 1);
            % for each input from the time series
            for ord=1:100
                timeout = CYCLES_PER_SEC;
                if shit == 3
                    timeout = 1
                end
                
                % reset feedforward part of the network
                %self.activation(self.ffwd_ids) = 0;
                %self.net_input_avg(self.ffwd_ids) = -10;
                
                % default output is timeout
                %output_id = self.unit_id('timeout');
                %RT = timeout;
                
                % simulate response to stimulus
                responded = false;
                is_settled = false;
                settle_cycles = 0;
                for cycle=1:timeout
                    % set input activations
                    activation(input_ids) = 0;
                    if is_settled
                        activation(input_ids) = stimuli(ord, :);
                    end
                    
                    % log activation for plotting
                    activation_log(cycles + cycle, :) = activation;
                    net_log(cycles + cycle, :) = net_input;
                    
                    % see if network has settled
                    if cycle > SETTLE_LEEWAY && ~is_settled
                        from = cycles + cycle - SETTLE_LEEWAY + 1;
                        to = cycles + cycle - 1;
                        m = activation_log(from:to,:) - activation_log(from-1:to-1,:);
                        m = abs(mean(m, 2));
                        %fprintf('%d -> %.6f, %6f\n', cycle, mean(m), std(m));
                        if mean(m) < SETTLE_MEAN_EPS && std(m) < SETTLE_STD_EPS
                            % it has settled
                            is_settled = true;
                            settle_cycles = cycle;
                            % save stimulus onset
                            onsets = [onsets; cycles + cycle];
                        end
                    end
                    
                    % calculate net inputs for all units
                    net_input = activation * weights ...
                        + bias; % TODO REVERT ME + normrnd(0, 0.1, N);
                                        
                    % update activation levels for feedforward part of the
                    % network
                    net_input_avg = TAU * net_input + (1 - TAU) * net_input_avg;
                    
                    if shit == 1
                        activation = logistic(net_input);
                    elseif shit == 2
                        activation = logistic(net_input_avg);
                    elseif shit == 3
                        activation(input_ids) = X(ord, :);
                        activation(hidden_ids) = logistic(activation(input_ids) * weights(input_ids, hidden_ids) + bias(hidden_ids));
                        activation(output_ids) = logistic(activation(hidden_ids) * weights(hidden_ids, output_ids) + bias(output_ids));
                    end

                    
                    
                    %% SIMULTANEOUS UPDATE!!! #FUCKUP
                    
                    [max_act, max_idx] = max(activation(output_ids), [], 2);
                    output = activation(output_ids) == max_act;
                    
                    
                end % for cycle = 1:timeout
                
                %output = self.units{output_id};
                
                offsets = [offsets; cycles + cycle];
                if ~is_settled 
                    % we never settled => no stimulus onset (put it same as
                    % offset for consistency
                    onsets = [onsets; cycles + cycle];
                end
                %responses = [responses; {output}];
                %RTs = [RTs; RT];
                cycles = cycles + cycle;

                expected = correct(ord, :);

                if output == expected
                    hits = hits + 1;
                    is_hit(ord) = 1;
                else
                    misses = misses + 1;
                end
                fprintf('hits = %d, misses = %d, hit ratio = %f\n', hits, misses, hits / (hits + misses));
            end % for ord in stimuli
            
            activation_log(cycles:end,:) = [];
            accumulators_log(cycles:end,:) = [];
            net_log(cycles:end,:) = [];
            weights
            bias
            

            
            
    act = activation_log;
            
    figure;

    %x_lim = [onsets(last_trial) - 1000 onsets(last_trial)];
    x_lim = [0 10000];
    y_lim = [MINIMUM_ACTIVATION - 0.1 MAXIMUM_ACTIVATION + 0.1];
    onset_plot = onsets;
    offset_plot = offsets;
    
    subplot(3, 1, 1);
    plot(act(:, output_ids));
    legend({'a', 'b'});
    title('Outputs');
    xlim(x_lim);
    ylim(y_lim);
    line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
    line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

    subplot(3, 1, 2);
    plot(act(:, hidden_ids));
    legend({'a', 'b'});
    title('Hiddens');
    xlim(x_lim);
    ylim(y_lim);
    line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
    line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

    subplot(3, 1, 3);
    plot(act(:, input_ids));
    legend({'a', 'b'});
    title('Inputs');
    xlim(x_lim);
    ylim(y_lim);
    line([onset_plot onset_plot],y_lim,'Color',[0.5 0.5 0.5])
    line([offset_plot offset_plot],y_lim, 'LineStyle', '--', 'Color',[0.5 0.5 0.5])

            
end