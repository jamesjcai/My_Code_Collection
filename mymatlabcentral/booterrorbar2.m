function [resc,resl,resu,Pmwu,Pkst]=booterrorbar2(ka,lsg)
%
% for example, booterrorbar2(ka,lsg)

ulsg=unique(lsg);

n=length(ulsg);

resc=zeros(1,n);
resu=zeros(1,n);
resl=zeros(1,n);


bootn=100;

for k=1:n
    
    kaa=ka(lsg==ulsg(k));
    len=length(kaa);
    
    replis=bootstrp(bootn,@nanmedian,kaa);
    REP{k}=replis;
    resc(k)=2*nanmedian(kaa)-nanmedian(replis);
    
%    resc(k)=2*median(Dxa)-median(replis);
%    resl(k)=resc(k)-1.96*(std(replis)/sqrt(len));
%    resu(k)=resc(k)+1.96*(std(replis)/sqrt(len));
    
    [res]=bootci(bootn,@nanmedian,kaa);
    if isnan(res(1))
        [res]=bootci(bootn,@nanmedian,kaa(kaa<Inf));
    end
    if isnan(res(1))
        brep=bootstrp(bootn,@nanmedian,kaa(isfinite(kaa)));
        res=repli2ci(brep(isfinite(brep)));        
    end
    resl(k)=res(1);
    resu(k)=res(2);    
end

if nargout<1
    %figure;
    bar(1:n,resc)
    hold on
    %h=errorbar(1:n,resc,resl,resu,'r.');
    h=errorbar(1:n,resc,resc-resl,resu-resc,'r');
    set(h,'linestyle','none');
    %colormap summer
    colormap(vivid(2))
    hold off
end

if nargout>3
Pmwu=ones(n)*-1;
Pkst=ones(n)*-1;
for i=1:n-1
for j=i+1:n
    X=REP{i};
    Y=REP{j};    
    Pmwu(i,j)=ranksum(X,Y);   %  Mann-Whitney U
    [h,p]=kstest2(X,Y);
    Pkst(i,j)=p;
end
end
end


