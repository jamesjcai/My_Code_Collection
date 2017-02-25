function [ihs,pos]=getdownloadedihs(popid,chrid)

switch (upper(popid))
    case ('CEU')
	popid='ceu';
    case ('YRI')
	popid='yri';
    case ('JNC')
	popid='asn';
end

[a,pos,c,ihs]=textread(sprintf('y:/iHs_data/res/%s.ch%d',popid,chrid),...
    '%s%d%s%f','delimiter','\t','emptyvalue',NaN);
pos=pos(ihs~=-99);
ihs=ihs(ihs~=-99);
