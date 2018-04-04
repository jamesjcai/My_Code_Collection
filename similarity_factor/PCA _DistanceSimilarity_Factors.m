load hald
%   ingredients      13x4               416  double
X1=ingredients; %%treat as Historical data
X2=ingredients+rand(size(ingredients,1),size(ingredients,2)); %%%treat as  Snapshot data


[coeff1,score1,latent1,tsquared1,explained1] = pca(X1);
s=0
for i=1:1;length(explained1)
    s=s+explained1(i);
    if s>0.95
        break;
    end
end
i=i+1;
% Principal component of Snapshot
L=score1(:,1:i) 

[coeff2,score2,latent2,tsquared2,explained2] = pca(X2);
s=0
for i=1:1;length(explained2)
    s=s+explained2(i);
    if s>0.95
        break;
    end
end
i=i+1;

%%%principal component of Historical
M=score2(:,1:i)

% number of component
index=i;

%  Equation 1
result1=trace(L'*M*M'*L)/index  

ss=0
for i=1:1:index
    for j=1:1:index
        a=L(:,i);
        b=M(:,j);
%  According to the paper, Equation 2 should be calculated as 
%         ss=ss+(dot(a,b)/(norm(a)*norm(b))).^2; 
% infact Equation 2 should be calculated as follow,then it is equal to the result of Equation 1
       ss=ss+(dot(a,b)).^2;
    end
end
ss=ss/index %%%result of Equation 2

xh=mean(X2)
xs=mean(X1)

%%Calculation of seudo-inverse of Snapshot
A=X1;
 [U,S,V] = svd(A); % A = U*S*V'
% Step2: 将S中的非零元素求倒
 T=S;
 T(find(S~=0)) = 1./S(find(S~=0));
% Step3: 求invA
svdInvA = V * T' * U';
xfc=T(1:min(size(X1)),1:min(size(X1)))

% Equation 3
fi=sqrt((xh-xs)*xfc*(xh-xs)')

syms x;
e4=sqrt(2/pi)*int(exp(-x.^2/2),x,fi,inf)
result=vpa(e4,4) %%%Equation 4

