function [seq,bas] = chrnuc_chimp(chrid,startn,endn)
if (nargin<3), endn=startn; end

if (isa(chrid,'memmapfile'))
    chrmm=chrid;
else
    %mmFilename=sprintf('Y:/chimpaln/mmFilenamechr%d',chrid);
    mmFilename=sprintf('C:/biodata/Alignments/chimpalnPanTro2/download/seq/mmf/alnPanTroIIchr%d.mm',chrid);
    chrmm = memmapfile(mmFilename, 'format', 'uint8');
end

if isscalar(startn)
seq=chrmm.Data(startn:endn);
else
seq=chrmm.Data(startn);
end

clear chrmm
seq=seq';
if (nargout>1),
  bas=int2nt(seq);
end

