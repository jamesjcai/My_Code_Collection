function region=regioncoord(chrid,regtype,within)
% intergenic intron inalu outalu

if (nargin<3)
    chrid=1;
	regtype='gene';
	within=0;
end

conn = database('snpdata', '', '');
%ping(conn)

T=tilingschema;
T=T(find(T(:,2)==chrid),:);
%T=T([1:5],:);

region=[];

for (k=1:length(T)),
    regionid=T(k,1);
    startn=T(k,3); endn=T(k,4);
    switch lower(regtype)
        case {'gene'}
    	sqldb='gene';
        case {'alu'}
    	sqldb='alu';
    end
sql=sprintf('SELECT distinct chr_id,start,end FROM feacoord f where type=''%s'' and chr_id=%d and start >=%d and end <=%d;',sqldb,chrid,startn,endn);
curs = exec(conn, sql);
setdbprefs('DataReturnFormat','structure')
curs = fetch(curs);
D=curs.Data;
close(curs)

if (within)
	if (isstruct(D))
		intergencoord=[D.start,D.xEnd];
		region=[region;[ones(size(intergencoord,1),1)*regionid, ones(size(intergencoord,1),1)*chrid, intergencoord]];
    else
   		region=[region; [regionid,chrid,startn,endn]];
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
