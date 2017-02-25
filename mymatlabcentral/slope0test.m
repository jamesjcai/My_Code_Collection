function [slope,p]=slope0test(x,y)

xtmp=[x ones(length(x),1)]; %input matrix for regress function
ytmp=y;

%regression coefficients
[p,pINT,R,Rint] = regress(ytmp,xtmp);

xtmp(:,2)=[]; %delete column 2
%save coefficients value
m(1)=p(1); q(1)=p(2);

n=length(xtmp); 
xm=mean(xtmp); xsd=std(xtmp);

%regression standard error (RSE)
RSE=realsqrt((sum(R.^2))/(n-2)); %RSE

%standard error of regression coefficients
%cv=tinv(0.975,n-2); %Student's critical value
cv=tinv(0.9995,n-2); %Student's critical value
m(2)=(pINT(3)-p(1))/cv; %slope standard error
%m=[m pINT(1,:)]; %add slope 95% C.I.
%q(2)=(pINT(4)-p(2))/cv; %intercept standard error
%q=[q pINT(2,:)]; %add intercept 95% C.I.

[rp,pr]=corrcoef(xtmp,ytmp);
slope=m(1);
p=pr(2);

    %test on slope
    t=abs(m(1)/m(2)); %Student's t
    disp('Student''s t-test on slope=0')
    fprintf('t = %0.4f    Critical Value = %0.4f     p = %0.4f\n',t,cv,pr(2))
    if t>cv
        disp('Test passed: slope ~= 0')
    else
        disp('Test not passed: slope = 0')
    end

   