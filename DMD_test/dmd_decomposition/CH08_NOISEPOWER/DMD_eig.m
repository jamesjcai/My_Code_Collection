function lambda = DMD_eig(Yin, flavor)

% build data matrices
Y = Yin(:, 1:end-1);
Yp = Yin(:, 2:end);

r = size(Yin, 2);

switch lower(flavor),
    case 'dmd',
        [U, S, V] = svd(Y, 'econ');
        Atilde = U' * Yp * V / S;
        
    case 'fbdmd',
        % compute forward DMD
        [U, S, V] = svd(Y, 'econ');
        f_Atilde = U' * Yp * V / S;

        % compute backward DMD
        [U, S, V] = svd(Yp, 'econ');
        b_Atilde = U' * Y * V / S;

        % estimate forward/backward DMD
        Atilde = (f_Atilde * inv(b_Atilde)) ^ 0.5;
        
    case 'tlsdmd',
        % stack
        Z = [Y; Yp];

        % tls DMD
        [U, ~, ~] = svd(Z, 'econ');
        U11 = U(1:r, 1:r);
        U21 = U(r+1:end, 1:r);

        Atilde = U21 * inv(U11);

    otherwise,
        error('WTF, do not understand which algorithm to execute.');
        Atilde = [];
end;

[~, D] = eig(Atilde);
lambda = diag(D); % eigenvalues
