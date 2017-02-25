function [bins]=buildequalbins(x)

x=sort(x);
totaln=length(x);
stepn=round(totaln/100);

k=1:stepn:totaln;
bs=x(k);

bins=[];
for (k=1:length(bs)-1),
     bins=[bins; [bs(k) bs(k+1)]];
end

