n=100;

y=[];
x=[];
g=[];
for k=1:5
    t=rand(n,1);
    x=[x; t];
    y=[y; (k+1)*t*3.0+(1.0+0.5*(k-1))+randn(n,1)*0.99];
    g=[g;k*ones(n,1)];    
end

%%
fitlm([x g],y,'y~1+x1*x2')
tbl=table(y,x,g);
fitlm(tbl,'y~1+x*g')

%%
g=categorical(g);
tbl=table(y,x,g);
fitlm(tbl,'y~1+x*g')

%%
lme=fitlme(tbl,'y~1+x+(x|g)')

%%
X=[ones(size(x)),x];
Z=[ones(size(x)),x];
lme1=fitlmematrix(X,y,Z,g);

%%
FixedEff=fixedEffects(lme);
RndEff=randomEffects(lme);

%%
figure;
ax1=gscatter(x,y,g);
refline(FixedEff(2)+RndEff(2),FixedEff(1)+RndEff(1))
refline(FixedEff(2)+RndEff(4),FixedEff(1)+RndEff(3))
refline(FixedEff(2)+RndEff(6),FixedEff(1)+RndEff(5))
refline(FixedEff(2)+RndEff(8),FixedEff(1)+RndEff(7))
refline(FixedEff(2)+RndEff(10),FixedEff(1)+RndEff(9))





