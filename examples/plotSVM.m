mat = readmatrix("data.txt");
X = mat(1:end-1, 1:2);
lab = mat(1:end-1, 3);
idx = (lab == 1);
scatter(X(idx, 1), X(idx, 2), 'r');
hold on
idx = (lab == -1);
scatter(X(idx, 1), X(idx, 2), 'b');

w = mat(end, :);
f = @(x) -w(1)/w(2) * x - w(3)/w(2);
Xmin = min(X(:,1));
Xmax = max(X(:,1));
fplot(f, [Xmin, Xmax], 'LineWidth', 1.5);