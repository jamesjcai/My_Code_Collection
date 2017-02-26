function [h,hraw]=genomicfaywuh(chrid,startn,endn,pos,freq)

%chrid=10;
%startn=8000000;
%endn=startn+400000;
if (nargin<5)
  [pos,freq]=textread(sprintf('Y:/perlegen/ancestralfreq/amaskedsancAfricanPerlegenchr%d.tab',chrid),'%d%f');
  pos=pos(freq<1 & freq>0);
  freq=freq(freq<1 & freq>0);
end
h=nan;
hraw=nan;

freq=freq(pos>=startn & pos<=endn);
if isempty(freq)
    return;
end

q=freq;
[Sn]=length(q);

if (Sn>0)
	p=1-q;   % derived freq
	smpln=48;

	theh=sum(2.*p.*p)*(smpln/(smpln-1));
	thepi=sum(2.*p.*q)*(smpln/(smpln-1));

	hraw=thepi-theh;

	nx=1:(smpln-1);
	an=sum(1./nx);
	bn=sum(1./(nx.^2));
	bn2=sum(1./([1:smpln].^2));
	t1=Sn./an;
	t2=Sn*(Sn-1)./(an.^2+bn);

	n=smpln;
	hvar=t1*(n-2)/(6*(n-1))+t2*((18*n^2)*(3*n+2)*bn2-(88*n^3+9*n^2-13*n+6))/(9*n*(n-1)^2);
    
	h=hraw./sqrt(hvar);


%nsam=smpln
%nx=1:(nsam-1);
%th=2*sum((nx.^2).*(thissfs))./(nsam*(nsam-1));
%tp=2*sum((nx.*(nsam-nx)).*(thissfs))./(nsam*(nsam-1));  %thetapi
%hraw=tp-th;


end
