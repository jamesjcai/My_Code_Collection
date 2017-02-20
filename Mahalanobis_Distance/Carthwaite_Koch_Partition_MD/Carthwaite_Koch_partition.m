function [d,f]=Carthwaite_Koch_partition(Y,X)

m=mean(X,1);
SampVar=cov(X);
S=SampVar^-1;
C=Y-m;

% Start that part of the Garthwaite-Koch partition that is the same for all items 
Dmat = diag(1./sqrt(diag(SampVar)));
DSD=Dmat*SampVar*Dmat;

[eigvec,eigval]=eig(DSD);
InvRootEig=1./sqrt(diag(eigval));

InvDSDhalf=eigvec*diag(InvRootEig)*eigvec';

d=zeros(size(Y,1),1);
for k=1:size(Y,1)
    df=C(k,:);
	d(k)=df*S*df';
	W=InvDSDhalf*Dmat*df';
	GKcontrib = diag(W*W');
	f=GKcontrib./sum(GKcontrib);
	% Correlation = diag(inv(InvDSDhalf));
end
