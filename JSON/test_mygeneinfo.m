s1=urlread('http://mygene.info/v3/query?q=ATMIN');
s2=urlread('http://myvariant.info/v1/query?q=rs58991260');
% d2=JSON.parse(s2)
d2=parse_json(s1);


addpath('jsonlab')
dat=loadjson(s1);