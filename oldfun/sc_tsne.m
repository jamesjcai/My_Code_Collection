function [s]=sc_tsne(X,genelist)
s=tsne(X');

%figure; scatter(s(:,1),s(:,2)); xlabel('tSNE\_1'); ylabel('tSNE\_2');
%title('tSNE (Perplexity = 30)');

a=uifig_scatter;
a.X=X;
a.genelistx=genelist;
s=tsne(X');
a.score=s;
% scatter(a.UIAxes,s(:,1),s(:,2));
gscatter(a.UIAxes,s(:,1),s(:,2),X(10,:));

