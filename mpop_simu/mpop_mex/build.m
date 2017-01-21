%mex mpop_mex.cpp haplotype.cpp population.cpp C:/gsl-1.13-windows-binaries/gsl/lib/cblas_d.lib C:/gsl-1.13-windows-binaries/gsl/lib/gsl_d.lib -IC:/gsl-1.13-windows-binaries/gsl/include -IC:/boost_1_43_0/

if 1<0
if strcmp(mexext,'mexw64')
    mex mpop_mex.cpp haplotype.cpp population.cpp C:/gsl-1.13-windows-binaries/gsl/lib_x64/cblas.lib C:/gsl-1.13-windows-binaries/gsl/lib_x64/gsl.lib -IC:/gsl-1.13-windows-binaries/gsl/include -IC:/boost_1_44_0/
elseif strcmp(mexext,'mexw64')
    mex mpop_mex.cpp haplotype.cpp population.cpp C:/gsl-1.13-windows-binaries/gsl/lib_w32/cblas.lib C:/gsl-1.13-windows-binaries/gsl/lib_w32/gsl.lib -IC:/gsl-1.13-windows-binaries/gsl/include -IC:/boost_1_44_0/    
end
end



%%
tic
N=2000;
[a,b]=ms_mex(N,4,0);

a=logical(a);
c.N=N;
c.m=0.001;
c.r=0.001;
c.s=0.01;
c.h=0.5;
c.g=200;
c.e=ceil(99999*rand);
[x,y]=mpop_mex(a,b,c);

fprintf('\nmpop -N %d -m %f -r %f -S -s %f -h %f -g %d -i msin -o msout -e %d\n\n',...
        c.N,c.m,c.r,c.s,c.h,c.g,c.e);

toc
%return;
%%

%{
[a,b]=ms_mex(100,0,20);
a=logical(a);
c.N=100;
c.s=0.0;
c.h=0.5;
c.r=0.0;
c.m=0.0;
c.g=1;
c.e=ceil(9999999*rand);
[x,y]=mpop(a,b,c);
%}