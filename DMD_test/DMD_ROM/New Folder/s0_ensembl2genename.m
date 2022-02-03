a=tempname;
b=websave(a,'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt');
T=readtable(b,'FileType','text');
string([T.ensembl_gene_id, T.symbol, T.name])
%%
keySet=T.ensembl_gene_id;
valueSet=T.symbol;
M = containers.Map(string(keySet),string(valueSet));
%%
input=string({'ENSG00000121410', 'ENSG00000121411'});

%%
ix=isKey(M, cellstr(input));
values(M,cellstr(input(ix)));

