function act = logistic(net)
    act = 1 ./ (1 + exp(-net));
end