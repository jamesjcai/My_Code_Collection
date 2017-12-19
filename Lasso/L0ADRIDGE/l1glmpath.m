function [outopt, E, B] = l1glmpath(X, y, cv, model, penalty, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-validation with generalized linear model
% Written by Zhenqiu Liu
% Cedars-Sinai Medical Center
%  12/18/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if nargin < 6,
     L = 50;
 end

if nargin < 5,
   penalty = 'scad';           % set penalty to lasso

end

if nargin < 4,
   model = 'logistic';         % set model to logistic
end

if strcmp(penalty, 'scad')
    penparam = 3.7;
else
    penparam = 1;
end

penidx = true(size(X,2),1);  % leave intercept unpenalized
lambdastart = 0;            % find the maximum tuning parameter to start
for j=1:size(X,2)
    if (penidx(j))
    lambdastart = max(lambdastart, ...
        glm_maxlambda(X(:,j),y,model,'penalty',penalty,'penparam',penparam));
    end
end

lammin = 0.0001;
llamax = log(lambdastart);
llamin = log(lammin);
llam = llamin:((llamax -llamin)/L):llamax;
lam = exp(llam);
p =length(lam);

[n, m] = size(X);
cvidx = crossvalind('Kfold', n, cv);
y = y(:);
%cvidx = crossvalid('Kfold', y, cv); % when y is binary
mx  = max(cvidx);
m = size(X,2);

E = zeros(p, 1);

B = zeros(m, p);



for k =1:p, 
    Re = [];
    for i =1:mx,
        te = (cvidx ==i);
        tr = ~te;
        Xr = X(tr,:);
        yr = y(tr);
        Xt = X(te,:);
        yt = y(te);
         w = ...               % sparse regression
            glm_sparsereg(Xr,yr,lam(k),model,'penidx',penidx,'penalty',penalty,...
            'penparam',penparam);
      
         Et = logistestt(Xt, yt, w);
         Re = [Re; Et];
    end
    rmse = sqrt(mean(Re.^2));
    E(k) = rmse;
    bt = ...               % sparse regression
            glm_sparsereg(Xr,yr,lam(k),model,'penidx',penidx,'penalty',penalty,...
            'penparam',penparam);
    B(:, k) = bt;      

end
[e, j] = min(E);
 outopt.e = e
 outopt.b =B(:, j);
 outopt.lam = lam(j);
end


