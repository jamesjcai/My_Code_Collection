
D=[];
for k=1:20
    a=XM(:,:,k);
    D=[D a(:)];
end
clear XM
    ISize=size(a);
X1 = D(:,1:(end-1));
X2 = D(:,2:end);

%%
load networks/sce_data_smallest.mat
g=sce.g;
a=string(ls('networks/A*.mat'));
a=strtrim(natsortfiles(a));
D=[];
cd networks
for k=1:length(a)
    k
    load(a(k));
    A=A(1:500,1:500);
    A=ten.e_filtadjc(A,0.5,false);
    D=[D A(:)];
end
cd ..
ISize=size(A);
X1 = D(:,1:(end-1));
X2 = D(:,2:end);

% g=sce.g;
%%

tic;
NumberDMDModes = 19;
[ EigenVector , EigenValue ] = DMD( X1 , X2 , NumberDMDModes );
toc;

ISize=[500,500];


f = figure;
f.WindowState = 'maximized';
k = 1;

for i=1:3
    for j=1:6
        subplot(3,10,k);
        I = reshape( EigenVector(:,k) , ISize );
        imagesc( abs(I) );
        
        D = EigenValue(k,k);
        title( sprintf('EV = %0.2f + i %0.2f', real(D) , imag(D)) );
        
        k = k+1;
    end
end

sgtitle('First 30 DMD with eigen values');
I = reshape( EigenVector(:,1) , ISize );
[~,idx]=sort(sum(abs(I),2),'descend');
g(idx(1:10))
