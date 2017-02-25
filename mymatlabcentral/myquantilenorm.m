function [x]=myquantilenorm(x)
%x=[2 4 4 ; 5 4 14; 4 6 8; 3 5 8; 3 3 9];
[xs,ix]=sort(x);
t=mean(xs,2);
[~,iy]=sort(ix);
for k=1:size(x,2)
    x(:,k)=t(iy(:,k));
end
