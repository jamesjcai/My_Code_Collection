R=[cosd(75) -sind(75); sind(75) cosd(75)]
a=10*randn(1000,1);
b=randn(1000,1);
G=R*[a b]';
G=G';
plot(a,b,'.'); hold on; plot(G(:,1),G(:,2),'.')
%%
Xs=[a b];
Xh=G;
pca_sim_factor(Xh,Xs)
