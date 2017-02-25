function region=regioncoord(regtype,within)
% intergenic intron inalu outalu

if (nargin<2)
	regtype='gene';
	within=0;
end



conn = database('snpdata', '', '');
%ping(conn)

T=tilingschema;
%T=T([1:5],:);

region=[];
oldchr=0;
for (k=1:length(T)),
      regionid=T(k,1);
      chrid=T(k,2);
      if (chrid~=oldchr),
	disp(sprintf('At chromosome %d\n',chrid));
	oldchr=chrid;
      end

startn=T(k,3); endn=T(k,4);


switch lower(regtype)
    case {'gene'}
	sqldb='gene';
    case {'alu'}
	sqldb='alu';
end

sql=sprintf('SELECT * FROM feacoord f where type=''%s'' and chr_id=%d and start >=%d and end <=%d;',sqldb,chrid,startn,endn);
curs = exec(conn, sql);
setdbprefs('DataReturnFormat','structure')
curs = fetch(curs);
D=curs.Data;
close(curs)


if (within)
	if (isstruct(D))
		intergencoord=[D.start,D.xEnd];
		region=[region;[ones(size(intergencoord,1),1)*regionid, ones(size(intergencoord,1),1)*chrid, intergencoord]];
	end
else
	if (isstruct(D))
		[intergencoord]=excludregions([D.start,D.xEnd],startn,endn);
		region=[region;[ones(size(intergencoord,1),1)*regionid,ones(size(intergencoord,1),1)*chrid, intergencoord]];
	else
		region=[region; [regionid,chrid,startn,endn]];
	end
end



end
close(conn)
