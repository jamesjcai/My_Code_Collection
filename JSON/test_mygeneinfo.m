s1=urlread('http://mygene.info/v3/query?q=ATMIN');
s2=urlread('http://myvariant.info/v1/query?q=rs58991260');
% d2=JSON.parse(s2)
d2=parse_json(s1);


addpath('jsonlab')
dat=loadjson(s1);



s3=urlread('http://epigenomesportal.ca/cgi-bin/api/getDataHub.py?session=1051');
s4=urlread('ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/json/locus_groups/protein-coding_gene.json');

%%
tic
dat4=loadjson(s4);
toc
%Elapsed time is 10002.092182 seconds.

tic
d4=parse_json(s4);
toc
%Elapsed time is 36.767753 seconds.

%%
% tic
% [data4] = parse_json1(s4);
% toc