function addregloglog(Dn,Pneu,showqudratic,linestyle)
%addreg(Dn,Pneu,showqudratic,linestyle)

if nargin<3
    showqudratic=0;
end
if nargin<4
    linestyle='g--';
end

Dn=Dn(~isnan(Pneu));
Pneu=Pneu(~isnan(Pneu));
Pneu=Pneu(~isnan(Dn));
Dn=Dn(~isnan(Dn));


%plot(Dn, Pneu,'.')
[orixy]=axis;


X=log10(Dn); Y=log10(Pneu);

stepl=(max(X)-min(X))./100;
x1=[min(X):stepl:max(X)];
[a]=regress1d(X,Y);
y1=a(1).*x1+a(2);

% log(P(k))=log(alpha) - gamma*log(k)

gammax=a(1);
logalpha=a(2);

step2=(max(Dn)-min(Pneu))./100;
x2=[min(Dn):step2:max(Pneu)];
y2=(10^logalpha).*x2.^gammax;

% P(k)=alpha*k^-gamma

hold on;

fprintf('P(k)=%.2f*k^%.2f\n', 10^logalpha, gammax);

%line(x,y)
loglog(x2,y2,linestyle,'LineWidth',2)
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

