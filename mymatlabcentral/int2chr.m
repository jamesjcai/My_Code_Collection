function [chrn] = int2chr(chrn)
%FUNCNAME - 

if ~(isstr(chrn)), 
chrn = round(real(chrn));
if ~(chrn>=1 & chrn<=24), error('Wrong Chromosome Number'); end
	if (chrn==23),
		chrn='X';
	elseif(chrn==24),
		chrn='Y';
	else
		chrn=int2str(chrn);
	end
else
	chromset={'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',...
		  '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'};
	chrn=upper(chrn);
	if ~(ismember(chrn,chromset)), error('Wrong Chromosome Number'); end
	if (strcmp(chrn,'23')),
		chrn='X';
	elseif (strcmp(chrn,'24')),
		chrn='Y';
	end
end
