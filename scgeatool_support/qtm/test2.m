R0  = FastPairMI(data, h);

R = R0(1:end-1,1:end-1)/(K-1);
J = R0(end,1:end-1);

tic;
fun = @(alpha)howmany(alpha,R,J) - K;
alphasol = fzero(fun,[0 1]);
% alphasol = 0.3;
[~,xsol] = howmany(alphasol,R,J);
g(find(xsol.BestX))
toc;

tic;
B = lasso(X',y');
g(B(:,55)~=0)
toc;

intersect(g(find(xsol.BestX)), g(B(:,55)~=0))


function [n,result] = howmany(alpha,R,J)
    Q = qubo((1-alpha)*R - alpha*diag(J));
    result = solve(Q);
    n = sum(result.BestX);
end


