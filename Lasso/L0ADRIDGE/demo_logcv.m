% rng('default') % for reproducibility
n = 100;
p = 1000;

S =100;

Et = zeros(S, 1);
optlam =zeros(S,1);
Bt = zeros(p+1, S);

for k =1:S,
    X = randn(n,p);             % generate a random design matrix
    X = [ones(size(X,1),1) X];
    b = zeros(p+1,1);           % true signalfirst ten predictors are 5
    b([1, 2, 6, 11]) = [0,  0.5, 0.5, -0.4]';
    % first 5 predictors are 1
    %b([16, 21, 26]) = [-1, -1, -1];               % next 5 predictors are -1
    inner = X*b + randn(n,1)/5;                % linear parts
    prob = 1./(1+exp(-inner));

    %y = double(prob > 0.5);
    y = prob;

    L=100;
    cv = 5;
    lammax = max([log(n), log(p)])+2;
    lammin = 0.0001;
    llamax = log(lammax);
    llamin = log(lammin);
    llam = llamin:((llamax -llamin)/L):llamax;
    lam = exp(llam);
    epn = 1e-4;

    [outopt, E, B] = l0glmpath(X, y, lam,  cv, epn);
    Et(k) = outopt.e;
    optlam(k) = outopt.lam;
    Bt(:,k) = outopt.b;
    k
end

save logcvout1000 optlam Bt Et;


