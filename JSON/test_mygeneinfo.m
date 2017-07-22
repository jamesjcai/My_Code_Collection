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

%%

s3=urlread('http://epigenomesportal.ca/cgi-bin/api/getDataHub.py?session=4123');
d=parse_json(s3);
a=fieldnames(d.samples);
%%
for k=1:length(a)    
    b=getfield(d.samples,a{k});
    %if ~strcmp(a{k},b.SAMPLE_ID)
    %       a{k}
    %end
    try
    fprintf('%s\t%s\t%s\t%s\t%s\n',a{k},b.donor_id,b.donor_age,...
        b.donor_sex,b.donor_ethnicity);
    catch
    end    
    try
    fprintf('%s\t%s\t%s\t%s\t%s\n',a{k},b.DONOR_ID,b.DONOR_AGE,...
        b.DONOR_SEX,b.DONOR_ETHNICITY);
    catch
    end

end
