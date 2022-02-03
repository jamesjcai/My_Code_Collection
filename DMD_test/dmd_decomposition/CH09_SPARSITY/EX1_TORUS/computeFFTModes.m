figure

for i=1:2 %real vs imag
    for k=1:K % all sparse vals
        subplot(2,K,(i-1)*K+k)
        x2tilde = 0*xtilde;
        x2tilde(I(k),J(k)) = (sqrt(-1))^(i-1); % = 1 for i=1 and =i for i=2
        x2 = ifft2(x2tilde);
        surf(X,Y,Z,real(x2),'EdgeColor','none');
        MODESFFT(i,k,:) = x2(:);
        view(10,32)
        colormap jet;
        axis off
        axis equal 
        drawnow
    end
end
set(gcf,'Position',[100 100 1440 450])