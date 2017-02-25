function [d]=genomictajimad(chrid,startn,endn,pos,freq)

%chrid=10;
%startn=8000000;
%endn=startn+400000;
if (nargin<5)
  [pos,freq]=textread(sprintf('Y:/perlegen/ancestralfreq/amaskedsancAfricanPerlegenchr%d.tab',chrid),'%d%f');
  pos=pos(freq<1 & freq>0);
  freq=freq(freq<1 & freq>0);
end
d=nan;

freq=freq(pos>=startn & pos<=endn);
if isempty(freq)
    return;
end

q=freq;
[Sn]=length(q);

if (Sn>0)
	p=1-q;   % derived freq
	smpln=48;

	%theh=sum(2.*p.*p)*(smpln/(smpln-1));
	thepi=sum(2.*p.*q)*(smpln/(smpln-1));

    d=tajima89d(smpln, Sn, thepi);
    
	

end
