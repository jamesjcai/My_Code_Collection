function Visualize3(X, Y, Z, E, F, G, k)

    if k==4
        k1=2; k2=2;
    elseif k==2
        k1=1; k2=2;
    end
    ndim=size(E,1);
    figure;
    
    X=X'; Y=Y'; Z=Z'; E=E'; F=F'; G=G';
    
    
    %Show the original manifolds
    subplot(k1,k2,1);
    plot3(X(:,1), X(:,2), X(:,3),'r-', 'LineWidth',2);
    hold on;
    plot3(Y(:,1), Y(:,2), Y(:,3),'b-', 'LineWidth',2);
    plot3(Z(:,1), Z(:,2), Z(:,3),'g-', 'LineWidth',2);
    title({['(A) Comparison of Manifold 1, 2 and 3 (Before Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
  
    % 3D
    subplot(k1,k2,2);
    if ndim>=3
    plot3(E(:,1), E(:,2), E(:,3),'r-', 'LineWidth',2);
    hold on;
    plot3(F(:,1), F(:,2), F(:,3),'b-', 'LineWidth',2);
    plot3(G(:,1), G(:,2), G(:,3),'g-', 'LineWidth',2);
    title({['(B) Comparison of Manifold 1, 2 and 3 (After 3D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
    end
    
    if k==4
    % 2D
    subplot(k1,k2,3);
    if ndim>=2
    plot(E(:,1), E(:,2), 'r-', 'LineWidth',2);
    hold on;
    plot(F(:,1), F(:,2), 'b-', 'LineWidth',2);
    plot(G(:,1), G(:,2), 'g-', 'LineWidth',2);
    title({['(C) Comparison of Manifold 1, 2 and 3 (After 2D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
    end
    
    % 1D
    subplot(k1,k2,4); 
    if ndim>=1;
    plot(E(:,1), 'r-', 'LineWidth',2);
    hold on;
    plot(F(:,1), 'b-', 'LineWidth',2);
    plot(G(:,1), 'g-', 'LineWidth',2);
    title({['(D) Comparison of Manifold 1, 2 and 3 (After 1D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
    end
    end

end