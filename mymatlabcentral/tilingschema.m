function [out]=tilingschema(tiling)

if nargin<1, tiling=1000000; end
%chromset={'1','10','11','12','13','14','15','16','17','18','19','2','20','21','22','3','4','5','6','7','8','9','23','24'};
chromset=[1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,23,24];

out=[];

switch (tiling)
    case (1000000)        
        id=321;
    case (500000)
        id=3409;         
end

for (k=1:length(chromset)),
	chrn=chromset(k);
	chrl=chrlen(chrn);
	x=1:tiling:chrl;
	for (p=1:length(x)-1),
        out=cat(1,out,[id, chrn, x(p), x(p+1)-1]);
		%fprintf('%d\t%s\t%d\t%d\n',id, chrn, x(p), x(p+1)-1) 
		id=id+1;
    end
    out=cat(1,out,[id, chrn, x(end), chrl]);
	%fprintf('%d\t%s\t%d\t%d\n',id, chrn, x(end), chrl)
	id=id+1;
end


