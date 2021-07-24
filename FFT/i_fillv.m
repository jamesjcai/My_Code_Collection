function [v,x]=i_fillv(v,x)
if ~isempty(x)
    n=round(length(v)/2);
    v1=v(1:n);
    v2=v(n+1:end);
    [m,idx]=maxk(x,2);
    v1(round(length(v1)/2))=m(1);
    v2(round(length(v2)/2))=m(2);
    x(idx)=[];
    v=[v1,v2];
    [v1,x]=i_fillv(v1,x);
    [v2,x]=i_fillv(v2,x);
end
end