load phenotype
i=~isnan(TOD);
load expr_b47_b11 data genelist

data_b47=data(:,1:210);
data_b11=data(:,211:end);


data_b47i=data_b47(:,i);
data_b11i=data_b11(:,i);
TODi=TOD(i);
Agei=Age(i);
P1Ii=P1I(i);
pHi=pH(i);
Racei=Race_isWhite(i);
Sexi=Sex_isMale(i);
RINi=RIN(i);


%%

[y,idx]=ismember('PER1',genelist);

figure;
subplot(2,1,1)
scatter(TODi,data_b11i(idx,:))
xline([0 12])
box on

subplot(2,1,2)
scatter(TODi,data_b47i(idx,:))
xline([0 12])
box on

%%

figure;
x=TODi;
y=data_b47i(idx,:)';
 plot(x,y,'o')
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,4.0],...
               'Upper',[5,5,13],...
               'StartPoint',[0.5 0.5 6.0]);         
ft = fittype('a*sin((pi/12)*x+b)+c',...
       'independent','x','options',fo);
[f, goodness, output] = fit(x,y,ft);
hold on
plot(f,x,y)
 ax=-10:0.1:20;
 plot(ax,f.a*sin((pi/12)*ax+f.b)+f.c)
xlim([-5 19])
goodness.sse




