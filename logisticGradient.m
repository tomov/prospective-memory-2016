function actGradient = logisticGradient(net)
    actGradient = logistic(net) .* (1 - logistic(net));
end