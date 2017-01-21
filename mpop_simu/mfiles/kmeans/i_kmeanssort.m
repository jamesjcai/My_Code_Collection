function [hap,d]=i_kmeanssort(hap,k)
if nargin<2
    k=2;
end

try
    
    
    %[idx]=kmeans(hap,k);
    
    
       X=hap; CUTOFF=2;
       Y = pdist(X,'euclid'); 
       Z = linkage(Y,'single'); 
       idx = cluster(Z,'maxclust',CUTOFF);
    
    
    [~,i]=sort(idx);
    hap=hap(i,:);
    d=sum(idx==1);
catch
    hap=sortrows(hap);
    d=0;
end