function d=markuniqueentries(x)
[a,b]=unique(x,'first');
d=zeros(1,length(x));
d(b)=1;

