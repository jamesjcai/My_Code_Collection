function [P3]=interprates2(P,wei)
%
[P2]=interprates(P);
index2=(wei<20000);
P2(index2)=0;
P3=interprates(P2);







