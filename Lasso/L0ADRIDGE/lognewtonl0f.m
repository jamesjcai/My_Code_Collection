function w = lognewtonl0f(X, y, lambda, itn, epn)
%l0 iterative logistic regression 
% Zhenqiu Liu
% 01/25/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n, m] = size(X);
        

if nargin < 5,
   epn =1e-6;
end

if nargin < 4,
    itn = 1000;
end

if nargin < 3,
    lambda = 1; % AIC criteria
end

w = zeros(m,1);

i=1;
Xt= X;
while i < itn,
    old_w = w; % s1 is nx1
    xw = X*w; 
    s1 = 1./(1+exp(-xw));
    a = s1.*(1-s1);
    A = diag(a);
    z = A*xw + (y-s1);
   
    if n > m,
       w = (Xt'*A*X + lambda*eye(m))\(Xt'*z);
    else
         w = Xt'*((A*X*Xt' + lambda*eye(n))\ z);
    end
    % w = w./max(abs(w));
    w2 = w.*w;
    w2 = w2';
    Xt = repmat(w2,n,1).*X;
    if norm(w - old_w) < epn,
       break;
    end
    i= i+1;    
end

w(abs(w) < 1e-4) =0;

end

