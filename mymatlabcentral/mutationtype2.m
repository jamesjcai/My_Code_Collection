function [typen,anccodon,codonpos]=mutationtype2(nuc,exons,strand,pos,ancnuc,mutnuc)

%{
exons=[869937      870043;
      870300      870389;
      870761      870896;
      871416      871529;
      871645      871788;
      873374      873475;
      873733      873846;
      876370      876481;
      877243      877382;
      877655      877843;
      878418      878531;
      879025      879135;
      879247      879325;
      881166      881256;
      881338      881458;
      882137      882268;
      882342      882516;
      884172      884324;
      884458      884483];

  pos=884324;
  pos=884458+10
  %pos=  870761;
  %}
  

%typen=zeros(size(pos));
maxexon=size(exons,1);
% at idx-th exon
idx=inladder2(exons,pos);
if ~(idx>0)
    %error('xxx')
    typen=0;
    anccodon=0;
    codonpos=0;
    return;
end
% all cds length
cdslen=exons(:,2)-exons(:,1)+1;


if (strand>0)
    if idx>1, leadinglen=sum(cdslen(1:idx-1)); else leadinglen=0; end
    poslen=leadinglen+pos-exons(idx,1)+1;
else    
    if idx<maxexon, leadinglen=sum(cdslen(idx+1:end)); else leadinglen=0; end
    poslen=leadinglen+exons(idx,2)-pos+1;    
end
nthcodon=floor(poslen/3)+(mod(poslen,3)>0);

anccodon=nuc([nthcodon*3-2,nthcodon*3-1,nthcodon*3]);
codonpos=mod(poslen,3);
if codonpos==0
    codonpos=3;
end

codonpos
%anccodon(codonpos)
%ancnuc
%mutnuc
%pause

if(anccodon(codonpos)~=ancnuc && anccodon(codonpos)~=mutnuc)
   % error('something wrong')
   %    anccodon
   %    codonpos
   %    ancnuc
   %    mutnuc
   %    pos
   %    exons
   %    strand   
   typen=0;
   anccodon=0;
else
    mutcodon=anccodon;
	mutcodon(codonpos)=mutnuc;
    anccodon(codonpos)=ancnuc;

	if (translateseq(mutcodon)==translateseq(anccodon))
	    typen=1;    
	else
	    typen=2;   %nonsyn
	end
end


