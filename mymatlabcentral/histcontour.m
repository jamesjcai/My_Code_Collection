function [xi,f]=histcontour(x,nbin,drawcurve,staircolor,ksdenstyle,linewid)

if nargin<6, linewid=2; end
if nargin<5, ksdenstyle='r-'; end
if nargin<4, staircolor=[0.50196 0.50196 0.50196]; end
if nargin<3, drawcurve=true; end
if nargin<2, nbin=50; end

if isempty(staircolor)
    staircolor=[0.50196 0.50196 0.50196];
end
x(isnan(x))=[];
n=length(x);

    [count,xout] = hist(x,nbin);  
    dx = xout(2)-xout(1);
    f = count/n/dx;  % Probability density    
    if ~strcmpi(staircolor,'w')
        h=stairs([xout-dx/2 xout(end)+dx/2],[f f(end)]);
        set(h,'color',staircolor)
    end

xlabel('X'), ylabel('Probability Density')
% title('Observed maximum of n exponentially distributed variables')

if drawcurve
    [f,xi] = ksdensity(x);
    %hold all
    %plot(xi,f,ksdenstyle,'linewidth',linewid)
    plot(xi,f,'-','linewidth',linewid,'color',[0.50196 0.50196 0.50196])    
    %hold off
end
end


%x = linspace(0,max(xout),500);
%for k = 1:length(n)
%    line(x,gevpdf(x,0,1,log(n(k))),'color','k')
    % Matlab <7: use   n(k)*exp(-x).*exp(-n(k)*exp(-x))
%end



