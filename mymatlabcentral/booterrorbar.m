function [resc,resl,resu,P]=booterrorbar(D,idxedD)

if nargin<2
    idxedD=0;
end
    
if ~(idxedD)    
if iscell(D)
    Dx=D;
else
    for k=1:size(D,2)
        Dx{k}=D(:,k);
    end
end
else
    idx=D(:,1);
    Data=D(:,2);
    for k=1:max(idx)
    Dx{k}=Data(idx==k,:);
    end    
end
    



n=length(Dx);

resc=zeros(1,n);
resu=zeros(1,n);
resl=zeros(1,n);


for k=1:n
    Dxa=Dx{k};
    len=length(Dxa);
    replis=bootstrp(1000,@median,Dxa);
    REP{k}=replis;
    resc(k)=2*median(Dxa)-median(replis);
    
%    resc(k)=2*mean(Dxa)-mean(replis);
%    resl(k)=resc(k)-1.96*(std(replis)/sqrt(len));
%    resu(k)=resc(k)+1.96*(std(replis)/sqrt(len));
    
    [res]=bootci(1000,@median,Dxa);
    resl(k)=res(1);
    resu(k)=res(2);    
end

bar(1:n,resc)
hold on
h=errorbar(1:n,resc,resl,resu,'r');
set(h,'linestyle','none');
colormap summer
hold off



P=ones(n);

for i=1:n-1
for j=i+1:n
    X=REP{i};
    Y=REP{j};    
    P(i,j)=signrank(X,Y);
end
end


