        load('data/hvg.mat','x');
        n = length(x(:,1));        
        mass = 1/n*ones(1,n);
        y0 = []; cut_indices0 = [];
        
        lambda1 = .001;
        lambda2 = 5;
        
        rho = 1.0;
        tol = 10^-2;
        max_m = 200;
        max_avg_turn = 30;
        normalize_data = 1;
        pause_bool = 0;
        
        tic;
        [yfinal,cut_indices,I,iters] = mppc(y0,cut_indices0,x,mass,lambda1,lambda2,tol,rho,...
            max_m, max_avg_turn,normalize_data,pause_bool);
        toc;
        z=yfinal;

        
        a=x;

lgu=a(:,1); dropr=a(:,2); lgcv=a(:,3);
xyz=a';

s = cumsum([0;sqrt(diff(lgu(:)).^2 + diff(dropr(:)).^2 ...
    + diff(lgcv(:)).^2)]);

pp1 = splinefit(s,xyz,15,0.75);
b = ppval(pp1,s);
       figure; 
       scatter3(a(:,1),a(:,2),a(:,3),'.'); 

       hold on
       scatter3(b(1,:)',b(2,:)',b(3,:)');
       scatter3(z(:,1),z(:,2),z(:,3)); 
