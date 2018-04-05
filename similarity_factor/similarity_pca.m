%dvenugopalarao%
function [TPCA TPCA1 PC1 PC2]=similarity_pca(X1,X2)
%% //////////////////////////////////////////////////////////////
% This function determines ordinary and weighted similarity between two data
% sets as X1 and X2. The programme also returns PCs of X1 and X2. TPCA and
% TPCA1 coressond to ordinary and weighted similarity 
% References:
% Ashish Singhal,Dale E. Seborg,Matching Patterns from Historical Data
% Using PCA and Distance Similarity Factors,Proceedings of the American Control Conference
% Arlington, VA June 25-27, 2001.
% Authors: 
% Seshu Kumar Damarla, Ph.D Student,Department of Chemical Engineering, National Institute of Technology Rourkela, India
% Prof.Madhusree Kundu, Associate Professor,Department of Chemical
% Engineering, National Institute of Technology Rourkela, India
%% /////////////////////////////////////////////////////////////
pi=22/7;
TPCA=[];
TPCA1=[];
PHI=[];

phi=[];
SPCA=[];
SPCA1=[];

size(X1) ;
%Find PCs of X1 and X1
[pc1,score1,latent1,tsquare] = pca(X1);
[pc2,score2,latent2,tsquare] = pca(X2);
pc1;
pc2;
PC1=pc1(:,1:2);
PC2=pc2(:,1:2);
% Find eigen values of X1 and X2

e1=latent1;
e2=latent2;
E1=sum(e1);
E2=sum(e2);
e11=(e1/E1)*100;
e12=(e2/E2)*100;
kn=size(e1);
k=0;
S=0;
for i=1:kn(1,1)
   S=S+e11(i);
   if S>95
       break
   end
   
end
k=i;
if k<2
    k=2;
else
    k=i;
end

% Consider 1st two eigen values corresponding two 1st two principle eigen  vectors
e1=e1(1:k);
e2=e2(1:k);


n1=size(PC1);
n2=size(PC2);


T=[];
for i=1:n1(1,2)
    th=[];
for j=1:n2(1,2)
a=PC1(:,i);
b=PC2(:,j);
b=b';
d=0;

for k=1:n1(1,1)
    
   d=d+a(k)*b(k) ;
end

 d1=sum(a.^2);
d2=sum(b.^2);
th=[th;(d/(d1*d2))];

end
T=[T th];
end
T;
T1=T.^2;

% Find ordinary similarity between X1 and X2
Spca=0;
for i=1:n1(1,2)
    
    for j=1:n2(1,2)
        
       Spca=Spca+(T(i,j) )^2;
        
    end
end
Spca=Spca/n1(1,2);
% Find weighted similarity between X1 and X2
Spca1=0;
for i=1:n1(1,2)
    
    for j=1:n2(1,2)
        
       Spca1=Spca1+e1(i)*e2(j)*(T(i,j))^2 ;
    end
    
end

l=0;
for k=1:n1(1,2)
    
   l=l+e1(k)*e2(k); 
end
Spca1=Spca1/l;


SPCA=[SPCA Spca ];
SPCA1=[SPCA1 Spca1];

TPCA=[TPCA;SPCA];
TPCA1=[TPCA1;SPCA1];
