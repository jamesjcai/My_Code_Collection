n = 500;
m = 10000;
Fb = [];
Fa =[];
Bl =[];
Bs =[];
Bm =[];
E1 =[];
E2 =[];
E3 =[];
E4 =[];
E5 =[];

for i =1:100,
    Z = randn(n, m);
    mu = exp(Z(:,[5 10 15])*[0.5 ;0.5; 0.4] + 1);
   %mu = exp(X(:,[5 10 15])*[0.5; 0.5; 0.4] );
   y = poissrnd(mu);
   b = zeros(m+1,1);
   b([1, 6, 11,  16]) =[1, 0.5, 0.5, 0.4];
    f = poisnewtonl0f(Z, y, log(n));
    p1 = glmval(f,Z,'log');
    e1 = y-p1;
    sq1 = sqrt((e1'*e1)/n);
    E1 =[E1; sq1];
    
    Fb =[Fb, f];
    f2 = poisnewtonl0f(Z, y, 2);
    p2 = glmval(f2,Z,'log');
    e2 = y-p2;
    sq2 = sqrt((e2'*e2)/n);
    E2 =[E2; sq2];
    Fa =[Fa, f2];
    
    
    X =[ones(n,1) Z];
    model = 'loglinear';        % set model to logistic
    penidx = true(size(X,2),1);  % leave intercept unpenalized
       
    penparam = 1;
    lambdastart = 0;          %b = zeros(p+1,1);
      
    penalty = 'enet';           % set penalty to lasso
  % find the maximum tuning parameter to start
  for j=1:size(X,2)
     if (penidx(j))
        lambdastart = max(lambdastart, ...
        glm_maxlambda(X(:,j),y,model,'penalty',penalty,'penparam',penparam));
     end
  end
 
   lambda = 0.5*lambdastart;   % tuning parameter value
    betahat1 = ...               % sparse regression
    glm_sparsereg(X,y,lambda,model,'penidx',penidx,'penalty',penalty, ...
    'penparam',penparam);
  
    p3 = glmval(betahat1,Z,'log');
    e = y- p3;
    sq3 = sqrt((e'*e)/n);
    E3 =[E3; sq3];
    Bl =[Bl, betahat1];
      
    penparam = 3.7;
    lambdastart = 0;          %
      
    penalty = 'scad';           % set penalty to lasso
  % find the maximum tuning parameter to start
  for j=1:size(X,2)
     if (penidx(j))
        lambdastart = max(lambdastart, ...
        glm_maxlambda(X(:,j),y,model,'penalty',penalty,'penparam',penparam));
     end
   end
  


    lambda = 0.5*lambdastart;   % tuning parameter value
    betahat2 = ...               % sparse regression
    glm_sparsereg(X,y,lambda,model,'penidx',penidx,'penalty',penalty, ...
    'penparam',penparam);

    p = glmval(betahat2,Z,'log');
    e = y- p;
    sq4 = sqrt((e'*e)/n);
    E4 =[E4; sq4];
    Bs = [Bs, betahat2];

    penparam = 1;
    lambdastart = 0;          %b = zeros(p+1,1);
      
    penalty = 'mcp';           % set penalty to lasso
  % find the maximum tuning parameter to start
  for j=1:size(X,2)
     if (penidx(j))
        lambdastart = max(lambdastart, ...
        glm_maxlambda(X(:,j),y,model,'penalty',penalty,'penparam',penparam));
     end
  end
 

  lambda = 0.5*lambdastart;   % tuning parameter value
  betahat3 = ...               % sparse regression
    glm_sparsereg(X,y,lambda,model,'penidx',penidx,'penalty',penalty, ...
    'penparam',penparam);
    p = glmval(betahat3,Z,'log');
    e = y- p;
    sq5 = sqrt((e'*e)/n);
    E5 =[E5; sq5];
  Bm =[Bm, betahat3];
      
 i
end
 Etot = [E1, E2, E3, E4, E5];
 save poisout110000 Etot Bl Bs Bm b Fb Fa;