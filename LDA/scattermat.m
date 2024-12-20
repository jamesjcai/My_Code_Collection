function [B, W]=scattermat(data,Y) 
%FUNCTION THAT CALCULATES SCATTER MATRIX: 
% B:BETWEEN CLASS SCATTER MATRIX 
% W:WITHIN CLASS SCATTER MATRIX 
%

[total_length, l]=size(data); %CALCULATE SIZE OF DATA 
clases=unique(Y); %GET VECTOR OF CLASSES 
tot_clases=length(clases); %HOW MANY CLASSES 
B=zeros(l,l); %INIT B AND W 
W=zeros(l,l); 
overallmean=mean(data); %MEAN OVER ALL DATA 
for i=1:tot_clases 
clasei = find(Y==clases(i)); %GET DATA FOR EACH CLASS 
xi=data(clasei,:); 

mci=mean(xi); %MEAN PER CLASS 
xi=xi-repmat(mci,length(clasei),1); %Xi-MeanXi 
% W=W+xi'*xi; %CALCULATE W 
W=W+((length(clasei)./total_length)*(xi'*xi)); 
% B=B+length(clasei)*(mci-overallmean)'*(mci-overallmean); %CALCULATE B 
B=B+((length(clasei)./total_length)*(overallmean-mci)'*(overallmean-mci)); 
end 
end

