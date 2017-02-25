function addregsemilogy(Dn,Pneu,showqudratic,linestyle)
%addreg(Dn,Pneu,showqudratic,linestyle)

if nargin<3
    showqudratic=1;
end
if nargin<4
    linestyle='r-';
end

Dn=Dn(~isnan(Pneu));
Pneu=Pneu(~isnan(Pneu));
Pneu=Pneu(~isnan(Dn));
Dn=Dn(~isnan(Dn));


%plot(Dn, Pneu,'.')
[orixy]=axis;

stepl=(max(Dn)-min(Dn))./100;
x=[min(Dn):stepl:max(Dn)];
[a]=regress1d(Dn,Pneu);
y=a(1).*x+a(2);
hold on;

%line(x,y)
semilogy(x,y,linestyle,'LineWidth',3)
%semilogy(x,y,linestyle,'LineWidth',2)
xlim([orixy(1) orixy(2)])
ylim([orixy(3) orixy(4)])

%{
try
[a]=L1LinearRegression(Dn, Pneu);
catch
[a]=L1LinearRegression(Dn', Pneu');
end
y=a(1)*x+a(2);
semilogy(x,y,'b-')
%}

%[r,p]=corr(Dn,Pneu,'type','spearman');
%fprintf('Spearman rho=%r (P=%f)\n',r,p);

if showqudratic
% Calculate fit parameters
p=polyfit(Dn,Pneu,2);
x2=linspace(min(Dn),max(Dn));
y2=polyval(p,x2);
%semilogy(x2,y2,'--','color',[0.25098     0.50196     0.50196]);
semilogy(x2,y2,'g--','LineWidth',2);
end

hold off;

