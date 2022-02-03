
% [X,g,c]=sc_readtsvfile('z:\Cailab\scTenifoldTime\Data\Cardiomyocytes\Cardiomyocytes.csv');
% T=readtable('Z:/Cailab/scTenifoldTime/Data/Cardiomyocytes/phenodata.csv');

[yes,idx]=ismember(c,T.X);
X=X(:,yes);
c=c(yes);
T=T(idx(yes),:);

%%
timepoint_levels = ["e14", "e18", "p0", "p1", "p4", "p8", "p11", "p14", "p15", "p18", "p22", "p28", "p35", "p56", "p84"];

valueset = [1:length(timepoint_levels)];
catnames = timepoint_levels;
%%
[~,timepoints]=ismember(string(T.Timepoint),catnames);
timepointc = categorical(timepoints,...
    valueset,catnames,'Ordinal',true);

[timepoints,idx]=sort(timepoints);
sce=SingleCellExperiment(X(:,idx),g,[],timepoints);

%%

load ../sce_data_smaller.mat
[t,idx]=sort(sce.list_cell_attributes{4});
X=sce.X(:,idx);
g=sce.g;

%[~,X,g]=sc_hvg(X,g,true,false);
%X=X(1:min([size(X,1),200]),:);
X=norm_libsize(X);

nsubsmpl=50;
csubsmpl=500;
[n,m]=size(X);
r=round(m/nsubsmpl);    % r = step length
winsize=max([r,csubsmpl]);
startptx=1:r:m;
while startptx(end)+winsize>m && r>1
    r=r-1;
    winsize=max([r,csubsmpl]);
    startptx=1:r:m;
    startptx=startptx(1:nsubsmpl);
end


gsettxt='TPM4;MYH7;CYC1;SLC8A1;COX8A;COX6A1;UQCRC1;TPM1;TNNT2;ATP2A2;ACTC1;TPM3;UQCRC2;UQCR11;MYL2;COX7A2;COX6B1;UQCRQ;COX5B;UQCR10;COX7A1;COX7C;ATP1B1;COX7B;UQCRFS1;COX6C;ATP1A1;COX6A2;RYR2;UQCRH;UQCRB;TNNI3;COX7A2L;COX4I1;MYL3;CACNA1C;TNNC1;COX5A';
[gset,idxa,idxb]=intersect(upper(g),string(strsplit(gsettxt,';')));

D=[];
    for k=1:nsubsmpl
        Xrep=X(idxa,startptx(k):startptx(k)+winsize-1);
        %Xrep=Xrep.';
        D=[D, Xrep(:)];
    end
    ISize=size(Xrep);
X1 = D(:,1:(end-1));
X2 = D(:,2:end);


tic;
NumberDMDModes = 30;
[ EigenVector , EigenValue ] = DMD( X1 , X2 , NumberDMDModes );
toc;




f = figure;
f.WindowState = 'maximized';
k = 1;

for i=1:3
    for j=1:10
        subplot(3,10,k);
        I = reshape( EigenVector(:,k) , ISize );
        imagesc( abs(I) );
        
        D = EigenValue(k,k);
        title( sprintf('EV = %0.2f + i %0.2f', real(D) , imag(D)) );
        
        k = k+1;
    end
end

sgtitle('First 30 DMD with eigen values');

[~,idx]=sort(sum(abs(I),2),'descend');
g(idx)
