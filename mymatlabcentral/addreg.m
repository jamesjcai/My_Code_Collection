function addreg(Dn,Pneu,showqudratic,linestyle,lineweight)
%addreg(Dn,Pneu,showqudratic,linestyle)

if nargin<5
    lineweight=2;
end

if nargin<4
    linestyle='r-';
end

if nargin<3
    showqudratic=false;
end

Dn=Dn(~isnan(Pneu));
Pneu=Pneu(~isnan(Pneu));
Pneu=Pneu(~isnan(Dn));
Dn=Dn(~isnan(Dn));


%plot(Dn, Pneu,'.')
[orixy]=axis;

x=[min(Dn),max(Dn)];
[a]=regress1d(Dn, Pneu);
y=a(1)*x+a(2);
hold on;

plot(x,y,linestyle,'LineWidth',lineweight,'LineSmoothing','on')
xlim([orixy(1) orixy(2)])
ylim([orixy(3) orixy(4)])

%{
try
[a]=L1LinearRegression(Dn, Pneu);
catch
[a]=L1LinearRegression(Dn', Pneu');
end
y=a(1)*x+a(2);
plot(x,y,'b-')
%}

%[r,p]=corr(Dn,Pneu,'type','spearman');
%fprintf('Spearman rho=%r (P=%f)\n',r,p);

if showqudratic
    % Calculate fit parameters
    p=polyfit(Dn,Pneu,2);
    x2=linspace(min(Dn),max(Dn));
    y2=polyval(p,x2);
    %plot(x2,y2,'--','color',[0.25098     0.50196     0.50196]);
    %plot(x2,y2,'g--','LineWidth',2);
    %plot(x2,y2,linestyle,'LineWidth',lineweight,'linesmoothing','on');
    plot(x2,y2,linestyle,'LineWidth',lineweight);
end
hold off;

