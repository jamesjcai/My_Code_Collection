load('fisheriris.mat')
X = meas(1:50,:);
X=zscore(X);

[rd1]=mahal_robust(X,X,0.75);
[rd2]=mahal_robust(X,X,0.75,[1 1 1 1]);
assert(norm(rd1-rd2)<1e-7)

[rd3]=mahal_robust(X(:,1:3),X(:,1:3),0.75);
[rd4]=mahal_robust(X,X,0.75,[1 1 1 0]);
[rd3 rd4]
assert(norm(rd3-rd4)<1e-7)
