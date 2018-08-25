S=[10 1; 9 1.5; 2 8; 1 9; 4 5];
A=[0.8 0.9 0.2];
A=[A;1-A];

%[Z,W]=fastICA(S*A,2)
[W1,H1]=nnmf(S*A,2)
addpath('nmfv1_4\')

%rng(1)
[W2,H2]=nmfrule(S*A,2)
%rng(1)
[W3,H3]=supersimple_nnmf(S*A,2)

addpath('nmf-toolbox\');
[W4, H4] = nmf(S*A,2)
