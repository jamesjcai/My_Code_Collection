load('GPL6244_Expr_intersect');
N=15507;
ff = @(x)fitnessf(x,C,P);
if max(size(gcp)) == 0 % parallel pool needed
    parpool % create the parallel pool
end
%options = gaoptimset('PlotFcns',@gaplotbestf);
gaoptions = gaoptimset('Display','iter','UseParallel',true,...
            'StallGenLimit',150);
for k=16:30
    tic;
[x, fval] = ga(ff,5,[],[],[],[],ones(1,5),...
          N*ones(1,5),[],1:5,gaoptions);
save(sprintf('solution_%d_%d',k,round(-fval)),'x','fval');
toc;
end
if max(size(gcp)) > 0 % parallel pool exists
    delete(gcp) % delete the pool
end
      

function [f]=fitnessf(idx,C,P)
    if length(unique(idx))<5   
        f=0;
    else
     d1=C(:,idx);
     d2=P(:,idx);    
     [d_obs]=mahal_robust(d1,d2);
     f=-sum(d_obs);
    end