function Y = SPCA_Algorithm(xMat, beta, k, n)
    obj1 = 0;
    obj2 = 0;
    thresh = 1e-50;
    V = eye(n);     % Identity matrix (500, 500)
    vMat = mat(V);
    for m = 1:10
        Z = -(xMat' * xMat) + beta * vMat;  % Same as Python
        [Z_eigVals, Z_eigVects] = eig(Z);
        eigValIndice = sort(Z_eigVals);
        n_eigValIndice = eigValIndice(1:k);
        n_Z_eigVect = Z_eigVects(:, n_eigValIndice);
        Q = double(n_Z_eigVect);  % Convert to double precision
        q = norm(Q, 2, 1);  % Norm along first dimension
        qq = 1.0 ./ (q * 2);
        VV = diag(qq);
        vMat = mat(VV);
        qMat = mat(Q);
        Y = xMat * qMat;
        % obj1 calculation commented out as in Python
        obj1 = norm(qMat);
        if m > 1
            diff = obj2 - obj1;
            if diff < thresh
                break
            end
        end
        obj2 = obj1;
    end
end

function Y = SPCA_cal_projections(X_data, beta1, k_d)
    % nclass = 2;
    % nclass = B_data.shape[1];
    n = length(X_data);  % Number of rows in X_data
    Y = SPCA_Algorithm(X_data', beta1, k_d, n);
end
