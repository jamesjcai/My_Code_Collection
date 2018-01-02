% genotype X
s_readX;

%%
% GRM Z
load grm_bin
n=3925;
Z = triu(ones(n),1);
Z(~~Z)=M.off;
Z=Z+Z';
Z(1:(n+1):end)=M.diag;

%%
f=0.5*nansum(X)./sum(~isnan(X));
GRM=zeros(n);
for i=1:100
    for j=i:100
        GRM(i,j)=0.5*nanmean((X(i,:)-2*f).*(X(j,:)-2*f)./(f.*(1-f)));        
    end
end 
GRM=GRM+triu(GRM,1)';


%%
% GRM Z
load grm_bin_01_04
n=3925;
Z = triu(ones(n),1);
Z(~~Z)=M.off;
Z=Z+Z';
Z(1:(n+1):end)=M.diag;

%%
s_readX;
f=0.5*nansum(X)./sum(~isnan(X));
idx=f>0.1 & f<0.4;
X=X(:,idx);
f=f(idx);

GRM=zeros(n);
for i=1:100
    for j=i:100
        GRM(i,j)=0.5*nanmean((X(i,:)-2*f).*(X(j,:)-2*f)./(f.*(1-f)));        
    end
end 
GRM=GRM+triu(GRM,1)';


