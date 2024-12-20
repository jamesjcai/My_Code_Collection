n=100;
y=[]; x=[]; g=[];
for k=1:5
    t=rand(n,1);
    x=[x; t];
    y=[y; (k+1)*t*3.0+(1.0+0.5*(k-1))+randn(n,1)*0.99];
    g=[g;k*ones(n,1)];    
end


%% treat catogrical varible G as continous  -- 3 free variables
fitlm([x g],y,'y~1+x1*x2')
tbl=table(y,x,g);
fitlm(tbl,'y~1+x*g')


%% set G as categorical variable -- # of free variables is 9
g=categorical(g);
tbl=table(y,x,g);
fitlm(tbl,'y~1+x*g')

%%
lme=fitlme(tbl,'y~1+x+(x|g)','fitmethod','REML')

%%
X=[ones(size(x)),x];
Z=[ones(size(x)),x];
lme1=fitlmematrix(X,y,Z,g);

%%
FixedEff=fixedEffects(lme)
RndEff=randomEffects(lme)

%%
figure;
ax1=scatter(x,y,'k.');
%hold on
%ax1=gscatter(x,y,g);
% refline(m,b) adds a reference line with slope m and intercept b to the current axes.
% refline(FixedEff(2)+RndEff(2),FixedEff(1)+RndEff(1))
% refline(FixedEff(2)+RndEff(4),FixedEff(1)+RndEff(3))
% refline(FixedEff(2)+RndEff(6),FixedEff(1)+RndEff(5))
% refline(FixedEff(2)+RndEff(8),FixedEff(1)+RndEff(7))
% refline(FixedEff(2)+RndEff(10),FixedEff(1)+RndEff(9))


%%

z1=lme.designMatrix('Random');
x1=lme.designMatrix('Fixed');
res=y-x*(x\y);

fun = @(b)parameterfun0(b,y,z1,x);
x0 = [12, 1, 0.1, 0.5, 0.5, 0.1, 0.5, 0.5, 0.1, 0.5, 0.5, 0.1];
[xe,fval]=fminunc(fun,x0);

[lme.fixedEffects;
lme.randomEffects]'

[x1\y; xe']'


function d = parameterfun(b,res,Z)
    d=sum((res-(Z(:,[1,3,5,7,9])*b([1 3 5 7 9])'+Z(:,[2, 4, 6 8 10])*b([2 4 6 8 10])')).^2);
end

function d = parameterfun0(b0,y,Z0,X)
    
    Z=Z0(3:end);
    b=b0(3:end);
    yp=X(:,1)*b0(1)+X(:,2)*b0(2)+(Z(:,[1,3,5,7,9])*b([1 3 5 7 9])'+Z(:,[2, 4, 6 8 10])*b([2 4 6 8 10])');
    d=sum((y-yp).^2);
end