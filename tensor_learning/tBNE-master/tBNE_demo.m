clear
clc

addpath(genpath('./tensor_toolbox'));
addpath(genpath('./FOptM'));
rng(5489, 'twister');

m = 10;
n = 10;
k = 10; % rank for tensor
[X, Z, Y] = tBNE_data(m, n, k); % generate the tensor, guidance and label

[T, W] = tBNE_fun(X, Z, Y, k);

[~, y1] = max(Y, [], 2);
[~, y2] = max(T{3} * W, [], 2);
fprintf('accuracy %3.2e\n', sum(y1 == y2) / n);