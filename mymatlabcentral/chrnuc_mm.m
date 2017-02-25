function [seq,bas] = chrnuc_mm(chrid,startn,endn)

if (isa(chrid,'memmapfile'))
    chrmm=chrid;
else    
    mmFilename=['C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.',int2str(chrid),'.mm'];
    chrmm = memmapfile(mmFilename, 'format', 'uint8');
end
seq=chrmm.Data(startn:endn);
clear chrmm
seq=seq';
if (nargout>1),
  bas=int2nt(seq);
end

