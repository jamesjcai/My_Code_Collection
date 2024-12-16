load D1_1f_G2M_5000g.mat
%y = sce.list_cell_attributes{2};
% X = sc_transform(sce.X);
% [~, X, g] = sc_splinefit(sce.X, sce.g);
%X = sc_norm(X);
X = full(sc_transform(sce.X, "type", ...
    "PearsonResiduals"));

idx = find(contains(sce.list_cell_attributes(1:2:end), 'monocle3_pseudotime'));
if isempty(idx), returen; end
y = sce.list_cell_attributes{idx*2};


% h = 1.2;
K = 15;

data = [X; y.'];

% tic
% R0  = FastPairMIpar(data, h);
% disp('R0par')
% toc
% 
% tic
% R0x  = FastPairMI(data, h);
% disp('R0x')
% toc


addpath('MI')
tic
if exist('D1_1f_R0.mat', 'file')
    load('D1_1f_R0.mat','R0')
else
    R0  = pairmi(data);
end
toc

R = R0(1:end-1,1:end-1)/(K-1);
J = R0(end,1:end-1);

tic;
fun = @(alpha)howmany(alpha,R,J) - K;
alphasol = fzero(fun,[0 1]);
% alphasol = 0.3;
[~,xsol] = howmany(alphasol,R,J);
g = sce.g;

g(find(xsol.BestX))
toc;

tic;
B = lasso(full(X'),y');
g(B(:,80)~=0)
toc;

intersect(g(find(xsol.BestX)), g(B(:,80)~=0))


function [n,result] = howmany(alpha,R,J)
    Q = qubo((1-alpha)*R - alpha*diag(J));
    result = solve(Q);
    n = sum(result.BestX);
end


