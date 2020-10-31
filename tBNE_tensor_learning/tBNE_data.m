function [X, Z, Y] = tBNE_data(m, n, k)
    B = randn(m, k);
    S = randn(n, k);
    A = {B, B, S};
    X = ktensor(A);

    Z = randn(n, 4);

    Y = zeros(n, 2);
    l = ceil(n / 2);
    Y(1 : l, 1) = 1;
    Y(l + 1 : end, 2) = 1;

    X = tensor(X);
end