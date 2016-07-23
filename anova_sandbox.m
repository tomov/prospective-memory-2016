n = 100;
y = rand(n, 1) + 0.1;
x = rand(n, 1);
z = [ones(n, 1); 2 * ones(n, 1)];
barweb([mean(x) mean(y)], [std(x) / sqrt(n) std(y) / sqrt(n)], 1, {});
[p, table] = anovan([x; y], z);
F = table{2, 6};
fprintf('p = %.5f, F = %.4f\n', p, F);