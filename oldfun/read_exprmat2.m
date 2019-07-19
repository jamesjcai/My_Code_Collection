function [X0,X1,genelist]=read_exprmat2(filename0,filename1,genecolnum0,genecolnum1)
[X0,genelist0]=read_exprmat(filename0,genecolnum0);
[X1,genelist1]=read_exprmat(filename1,genecolnum1);
[X0,genelist0]=filter_genes(X0,genelist0);
[X1,genelist1]=filter_genes(X1,genelist1);

[X0,X1,genelist]=exprmat_intersect(X0,X1,genelist0,genelist1);

[X0]=filter_samples(X0);
[X1]=filter_samples(X1);






    