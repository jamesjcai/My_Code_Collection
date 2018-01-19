%%
% GRM Z
load ../s1_mcode_make_grm/grm_bin
n=3925;
Z = triu(ones(n),1);
Z(~~Z)=M.off;
Z=Z+Z';
Z(1:(n+1):end)=M.diag;

%%

load test_inputs X y
X(isnan(X))=1;


%%
ys=zscore(y);
Y=ys*ys';  % outer product
lX=s_extract_triu(Z);
ly=s_extract_triu(Y);
f=fitlm(lX,ly)
% [ones(size(lX)) lX]\ly


%%
Y2=s_square_diff(ys);
lX=s_extract_triu(Z);
ly=s_extract_triu(Y2);
f=fitlm(lX,ly)
b=[ones(size(lX)) lX]\ly;
h2=-b(2)/b(1)





