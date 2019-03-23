load example10xdata.mat
[X,genelist]=sc_selectg(X,genelist);
% [X]=sc_norm(X,'type','deseq');
T=sc_hvg(X,genelist,false,false,true);
[~,i]=sort(T.fitratio,'descend');
i=i(1:750);
genelist=genelist(i);
X=X(i,:);
X=X(:,1:300);

% csvwrite('inputdata_example.csv',X);

%     object <- sc3_calc_dists(object)
%     object <- sc3_calc_transfs(object)
%     object <- sc3_kmeans(object, ks)
%     object <- sc3_calc_consens(object)
    
X=log2(X+1);
Dis=squareform(pdist(X'));
%Dis=1-corr(X,'type','s');
%Dis=1-corr(X,'type','p');

A=exp(-Dis./max(Dis(:)));   % adjacency matrix
xD=diag(sum(A).^-0.5);  % D=diag(sum(A)); % d(i) the degree of node i
xA=xD*A*xD;             % normalized adjacenty matrix
L=eye(size(A,1))-xA;    % also L=xD*(D-A)*xD 

% see https://people.orie.cornell.edu/dpw/orie6334/lecture7.pdf
% see https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian_2

[V,D]=eig(L);
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);
% [Vs,Ds]=eigs(L);

% For all six datasets, we found that the best clusterings were achieved when d was between 4% and 7% of the number of cells, N (Fig. 1c, Supplementary Fig. 3a and Online Methods). The robustness of the 4–7% range was supported by a simulation experiment
n=size(A,1);
drge=round(n.*[0.04 0.07]);
drange=drge(1):drge(2);

%%
clust=zeros(size(Vs,1),6);
for i=1:6
    clust(:,i) = kmeans(Vs,i,'emptyaction','singleton','replicate',5);
end
va=evalclusters(Vs,clust,'CalinskiHarabasz');
optimk=va.OptimalK;


MM=zeros(size(Vs,1));
cls=[];
for j=1:length(drange)
    idx=kmeans(Vs(:,1:drange(j)),3);
    cls=[cls; idx'];
    % [~,i]=sort(idx);
    % figure;
    % imagesc(X(:,i));
    % vline(sum(idx==1),'-y')
    % vline(sum(idx<3),'-y')

    % evalclusters
    M=zeros(size(Vs,1));
    for k=1:3
        a=find(idx==k);
        b=combnk(a,2);
        M(sub2ind(size(M),b(:,1),b(:,2)))=1;
        % https://blogs.mathworks.com/steve/2008/02/08/linear-indexing/
    end
    M=M+M';
    MM=MM+M;
end
MM=MM./10;





% https://books.google.com/books?id=RUO6BQAAQBAJ&pg=PA305&lpg=PA305&dq=Chu+et+al.+%5B7%5D+used+consensus+similarity+matrix+clustering+methods&source=bl&ots=StYL2MXyRS&sig=ACfU3U3HRbcclgiwsXufHJug3QfaH0QLlw&hl=en&sa=X&ved=2ahUKEwink8zb_4vhAhUSRqwKHcIoBq0Q6AEwAHoECAMQAQ#v=onepage&q=Chu%20et%20al.%20%5B7%5D%20used%20consensus%20similarity%20matrix%20clustering%20methods&f=false
% http://www.strehl.com/soft.html
%%

addpath('ClusterPack-V2.0\')
cl = clusterensemble(cls,3);

%%
figure;
h=clusion(MM,cl);
set(h,'color','g');

idx=cl;
[~,i]=sort(idx);

figure;
imagesc(MM(i,i));

figure;
imagesc(X(:,idx));
vline(sum(idx==1),'-y')
vline(sum(idx<3),'-y')
     
     