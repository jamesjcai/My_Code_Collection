function S=mat2spss(M)

S=[];

[n,m]= size(M);

for k=1:m
    x=[ones(n,1)*k,M(:,k)];
    S=[S;x];
end
    

