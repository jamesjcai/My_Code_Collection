% X    age, g - gene list, d - norm expr
s_data_common_skin;

%%
load \\165.91.29.77\disk2t\manuscript_workspace\JGuan_backup\PD\downloaded_data\GSE6613\SSMD\GO_Geneset.mat
N=length(DB.GeneSet);
pv=ones(N,1);
fv=ones(N,1);
rv=ones(N,1);
fid=fopen('res_mdmr_go_skin.txt','w');
for k=1:N
    k
    [yes,idx]=ismember(DB.GeneSet{k},g);    
    if sum(yes)>1 && sum(yes)<=50
        g1=DB.GeneSet{k}(yes);
        Y=d(idx(yes),:)';
        [pv(k),fv(k),rv(k)]=MDMR_test(X,Y);
        s=sprintf('%s,',g1{:});
        fprintf(fid,'%d\t%f\t%f\t%f\t%s\t%s\n',...
            k,pv(k),fv(k),rv(k),...
            cell2str(DB.GeneSetName{k}),s);
    end
end
fclose(fid);