figure
[U,S,V] = svd(XDAT,0);
stairs(cumsum(diag(S))/sum(diag(S)));
figure
for i=1:10
    subplot(2,5,i)
    ximg = reshape(U(:,i),n,n);
    surf(X,Y,Z,real(ximg),'EdgeColor','none');    
    view(10,32)
    colormap jet;
    axis equal
    axis off
    drawnow
end
set(gcf,'Position',[100 100 1440 450])