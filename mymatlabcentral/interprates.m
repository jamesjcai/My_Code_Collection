function [P2]=interprates(P)

[n,m]=size(P);
P2=zeros(n,m);
for (k=1:m)
r=P(:,k)';
idx=find(r==0);
x=1:n;
rx=[r;x];
rx(:,idx)=[];
inr=interp1(rx(2,:),rx(1,:),idx);
r(idx)=inr;
    
   P2(:,k)=r';
end






