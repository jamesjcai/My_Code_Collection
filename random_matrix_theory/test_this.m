load hg_pancreas_1459_cells.mat

%%
% [X]=sc_filterc(X);
% [X,genelist]=sc_filterg(X,genelist,7);

[X,genelist]=sc_selectg(X,genelist,2,1);
[X]=sc_norm(X,'type','deseq');

% web('https://rabadan.c2b2.columbia.edu/html/randomly/tutorial.html')
% MarchenkoPasturLaw;

[T]=sc_hvg(X,genelist,true,true);