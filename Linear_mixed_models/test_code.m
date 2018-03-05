% example from mixed-workshop_ulf.ppt

y=[21 19 20 22 14 15 13 16 14 17 15 17 12 11 12 14 16 20 18 19 14 14 14 12]';
brand=[1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2 2 2]';
site =[1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3]';

% see https://www.mathworks.com/help/stats/prepare-data-for-linear-mixed-effects-models.html

%mdl=fitlm([brand site],y,'y~x1*x2');
% brand=brand-1;
brand=categorical(brand);
site=categorical(site);

tbl = table(brand,site,y);
mdl = fitlm(tbl,'y~brand*site');

lme = fitlme(tbl,'y~brand+(brand|site)','FitMethod','ML');

%%


n=length(y);
r=length(unique(site));
[~,~,idx]=unique(site);
Z=zeros(n,r);

for k=1:n, Z(k,idx(k))=1; end
z1=sparse(Z);
z2=lme.designMatrix('Random');
% assert(isequal(z1,z2))


%x1=x2fx(brandx);
x1=lme.designMatrix('Fixed');
%assert(isequal(x1,x2))

% see https://www.mathworks.com/help/stats/relationship-between-formula-and-design-matrix-.html

%betaF=zeros(size(z1,2),1);   % lme.fixedEffects [22.9167 -4.7500]
%betaR=zeros(size(g1,2),1);   % lme.randomEffects [1.5119 -1.5857 0.0738]
%d=abs(y-(z1*betaF+g1*betaR));

fun0 = @(b)parameterfun0(b,y,x1,z1);
res=y-x1*(x1\y);
fun = @(b)parameterfun(b,res,z2);
% fun = @(b)sum((res-Z*b').^2);

x0 = [0.1, 0.5, 0.5, 0.1, 0.5, 0.5];
%options = optimset('LargeScale', 'on', 'Display', 'iter-detailed', ...
%    'TolX', 0.00001, 'TolFun', 0.001, 'GradObj', 'off', 'DerivativeCheck', 'off');
% [x,fval]=fminsearch(fun,x0);
[x,fval]=fminunc(fun,x0);


[lme.fixedEffects;
lme.randomEffects]'

[x1\y; x']'






%%
lme2 = fitlme(tbl,'y~brand+(1|site)+(brand-1|site)');
X=lme2.designMatrix('Fixed');
Z1=lme.designMatrix('Random');
Z2=lme2.designMatrix('Random');
figure;
subplot(1,2,1)
spy(Z2,'.k')
subplot(1,2,2)
imagesc(Z2)
colorbar
colormap bone(3)



%%
function d = parameterfun(b,res,Z)
    d=sum((res-(Z(:,[1,3,5])*b([1 3 5])'+Z(:,[2, 4, 6])*b([2 4 6])')).^2);
end


function d = parameterfun1(b,res,Z)    
    d=sum((res-Z*b').^2);
end


function d = parameterfun0(b,y,X,Z)
    betaF=X\y;
    betaR=b(1:3)';
    d=sum((y-(X*betaF+Z*betaR)).^2);
end
