function [J grad] = costFunction(nn_params, ...
                                   N, ...
                                   input_ids, ...
                                   hidden_ids, ...
                                   output_ids, ... 
                                   X,  ... % aka stimuli
                                   expected,  ... % aka correct (almost)
                                   lambda)

assert(isempty(intersect(input_ids, hidden_ids)));
assert(isempty(intersect(hidden_ids, output_ids)));
assert(isempty(intersect(input_ids, output_ids)));
assert(length(input_ids) + length(hidden_ids) + length(output_ids) == N);
                               
weights = reshape(nn_params(1:N*N), N, N);
bias = reshape(nn_params(N*N+1:end), 1, N);

m = size(X, 1);
num_labels = size(output_ids, 2);

% expected output



% feedforward 

activation = zeros(m, N);

% ---- instant propagation
activation(:, input_ids) = X;
net_input_to_hidden = activation(:, input_ids) * weights(input_ids, hidden_ids) + ones(m, 1) * bias(hidden_ids);
activation(:, hidden_ids) = logistic(net_input_to_hidden);
activation(:, output_ids) = logistic(activation(:, hidden_ids) * weights(hidden_ids, output_ids) + ones(m, 1) * bias(output_ids));

% ---- cascading propagation
net_input_avg = zeros(m, N);
for cycle=1:Model.CYCLES_PER_SEC / 2
    activation(:, input_ids) = X;
    net_input = activation * weights + ones(m, 1) * bias;
    net_input_avg = Model.TAU * net_input + (1 - Model.TAU) * net_input_avg;
    activation = logistic(net_input_avg);
    
    %activation(hidden_ids)
    %activation(hidden_ids)
end


% cost

Theta1 = [bias(hidden_ids)' weights(input_ids, hidden_ids)'];
Theta2 = [bias(output_ids)' weights(hidden_ids, output_ids)'];

J = 1.0/m * sum(sum(-expected .* log(activation(:, output_ids)) - (1 - expected) .* log(1 - activation(:, output_ids))));
                    
% regularize

TT1 = (Theta1 .* Theta1);
TT2 = (Theta2 .* Theta2);
reg = 0.5 * lambda / m * (sum(sum(TT1(:,2:end))) + sum(sum(TT2(:,2:end))));
J = J + reg;

% calculate the gradients

del = zeros(m, N);
del(:, output_ids) = activation(:, output_ids) - expected;
del(:, hidden_ids) = (del(:, output_ids) * weights(hidden_ids, output_ids)') .* logisticGradient(net_input_to_hidden);
%self.del(self.perception_ids) = (self.del(self.response_ids) * self.weights(self.perception_ids, self.response_ids)') .* self.logisticGradient(self.net_input_avg(self.perception_ids)); % TODO same here

weightsGradient = zeros(size(weights));
weightsGradient(hidden_ids, output_ids) = 1.0 / m * activation(:, hidden_ids)' * del(:, output_ids);
weightsGradient(input_ids, hidden_ids) = 1.0 / m * activation(:, input_ids)' * del(:, hidden_ids);
%self.weightsGradient(self.task_ids, self.response_ids) = 1.0 / m * self.wm_act(ismember(self.wm_ids, self.task_ids))' * self.del(self.response_ids);
%self.weightsGradient(self.input_ids, self.perception_ids) = 1.0 / m * self.activation(self.input_ids)' * self.del(self.perception_ids);

biasGradient = zeros(size(bias));
biasGradient(output_ids) = 1.0 / m * sum(del(:, output_ids));
biasGradient(hidden_ids) = 1.0 / m * sum(del(:, hidden_ids));
%self.biasGradient(self.perception_ids) = 1.0 / m * self.del(self.perception_ids);

% gradient regularization

weightsGradient(hidden_ids, output_ids) = lambda / m * weights(hidden_ids, output_ids) + weightsGradient(hidden_ids, output_ids);
weightsGradient(input_ids, hidden_ids) = lambda / m * weights(input_ids, hidden_ids) + weightsGradient(input_ids, hidden_ids);
%self.weightsGradient(self.task_ids, self.response_ids) = self.lambda / m * self.weights(self.task_ids, self.response_ids) + self.weightsGradient(self.task_ids, self.response_ids);
%self.weightsGradient(self.input_ids, self.perception_ids) = self.lambda / m * self.weights(self.input_ids, self.perception_ids) + self.weightsGradient(self.input_ids, self.perception_ids);

biasGradient(output_ids) = lambda / m * bias(output_ids) + biasGradient(output_ids);
biasGradient(hidden_ids) = lambda / m * bias(hidden_ids) + biasGradient(hidden_ids);
%self.biasGradient(self.perception_ids) = self.lambda / m * self.bias(self.perception_ids) + self.biasGradient(self.perception_ids);

% flatten to nn_params

grad = [weightsGradient(:); biasGradient(:)];
