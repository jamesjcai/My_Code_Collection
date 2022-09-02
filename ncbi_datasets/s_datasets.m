[a,b]=system('datasets.exe summary gene symbol Atmin brca1 --taxon human');
x=jsondecode(b);
x.genes(2).gene.description


