function histline(x,binum,islogy,islogx,linestyle)
if nargin<5, linestyle='rs-'; end
if nargin<4, islogx=false; end
if nargin<3, islogy=false; end
if nargin<2, binum=20; end

[a,b]=hist(x,binum);
if islogy && islogx
    loglog(b,a,linestyle)
    addregloglog(b,a);    
elseif islogy
    semilogy(b,a,linestyle)
elseif islogx
     semilogx(b,a,linestyle)
else
    plot(b,a,linestyle)
end


