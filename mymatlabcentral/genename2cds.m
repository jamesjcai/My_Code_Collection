function [transid,transcord,D] = genename2cds(genename)

%genename='ABCE1';
%genename='A4GALT';

if strcmp(upper(genename),'APEX')
    genename='APEX1';
elseif strcmp(upper(genename),'CKN1')
    genename='ERCC8';    
elseif strcmp(upper(genename),'DNCH1')
    genename='DYNC1H1';
elseif strcmp(upper(genename),'HCNP')
    genename='XAB2';
elseif strcmp(upper(genename),'NBS1')
    genename='NBN';
elseif strcmp(upper(genename),'REQ')
    genename='DPF2';
elseif strcmp(upper(genename),'DIA1')
    genename='CYB5R3';
elseif strcmp(upper(genename),'G22P1')
    genename='XRCC6';
elseif strcmp(upper(genename),'GAPD')
    genename='GAPDH';
elseif strcmp(upper(genename),'SEI1')
    genename='SERTAD1';
% sei1 --> SERTAD1
% G22P1 --> XRCC6    
% GAPD --> GAPDH
%DIA1 -->CYB5R3
 %DNCH1 - >DYNC1H1    
%HCNP --> XAB2
%NBS1 --> NBN
%REQ --> DPF2
%SEI1 --> CDK4  (wrong)
elseif strcmp(upper(genename),'CRF')
    genename='C1QL1';
elseif strcmp(upper(genename),'DO')
    genename='ART4';
elseif strcmp(upper(genename),'FSBP')
    genename='RAD54B';
elseif strcmp(upper(genename),'SCYA2')
    genename='CCL2';    
elseif strcmp(upper(genename),'SMP1')
    genename='TMEM50A';
%seattlesnps
%CRF --> C1QL1
%DO --> ART4
%FSBP  --> RAD54B
%SCYA2 --> CCL2
% smp1 -- >TMEM50A
end


conn = database('refgene', '', '');
%sql=sprintf('SELECT RefGene_prime.* FROM RefGene_prime WHERE RefGene_prime.name2=''%s'';',genename);
%curs = exec(conn, sql);      
%setdbprefs('DataReturnFormat','structure')
%curs = fetch(curs);
%D=curs.Data;
%close(curs)

%if (~isstruct(D))
	sql=sprintf('SELECT RefGene.* FROM RefGene WHERE RefGene.name2=''%s'';',genename);
	curs = exec(conn, sql);      
	setdbprefs('DataReturnFormat','structure')
	curs = fetch(curs);
	D=curs.Data;
	close(curs)
%end
close(conn)




if (isstruct(D))
  a=D.exonStarts{1};
  b=D.exonEnds{1}; 
  s1=D.cdsStart(1);
  e1=D.cdsEnd(1);
  transid=D.name{1};
  transcord=i_makecoord(a,b,s1,e1);  
else
  transid='';
  transcord=[];
end



function [data]=i_makecoord(a,b,s1,e1)
%%%%%%
[sx] = strread(a,'%d','delimiter',',');
[sy] = strread(b,'%d','delimiter',',');

data=[sx, sy];
data=sortrows(data);
id1=inladder2(data,s1);
id2=inladder2(data,e1);
%id1=sum(cdspos1>=e(:,1));
%id2=sum(cdspos2>=e(:,1));
data=data([id1:id2],:);
data(1)=s1;
data(end)=e1;

data(:,1)=data(:,1)+1;    %%%% <----------- important


if (mod(sum(data(:,2)-data(:,1)+1),3)>0)
	data(end)=e1-mod(sum(data(:,2)-data(:,1)+1),3);
end

