function Visualize2(X, Y, E, F, k)

    if k==4
        k1=2; k2=2;
    elseif k==2
        k1=1; k2=2;
    end
    figure;
    
    X=X'; Y=Y'; E=E'; F=F';
    
    
    %Show the original manifolds
    subplot(k1,k2,1);
    plot3(X(:,1), X(:,2), X(:,3),'r-', 'LineWidth',2);
    hold on;
    plot3(Y(:,1), Y(:,2), Y(:,3),'b-', 'LineWidth',2);
    title({['(A) Comparison of Manifold 1 and 2 (Before Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
  
    % 3D
    subplot(k1,k2,2);
    plot3(E(:,1), E(:,2), E(:,3),'r-', 'LineWidth',2);
    hold on;
    plot3(F(:,1), F(:,2), F(:,3),'b-', 'LineWidth',2);
    title({['(B) Comparison of Manifold 1 and 2 (After 3D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
    
    if k==4
    % 2D
    subplot(k1,k2,3);
    plot(E(:,1), E(:,2), 'r-', 'LineWidth',2);
    hold on;
    plot(F(:,1), F(:,2), 'b-', 'LineWidth',2);
    title({['(C) Comparison of Manifold 1 and 2 (After 2D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
  
    % 1D
    subplot(k1,k2,4);    
    plot(E(:,1), 'r-', 'LineWidth',2);
    hold on;
    plot(F(:,1), 'b-', 'LineWidth',2);
    title({['(D) Comparison of Manifold 1 and 2 (After 1D Alignment)']});
    xlabel('X', 'FontSize',14, 'FontWeight', 'bold');
    ylabel('Y', 'FontSize',14, 'FontWeight', 'bold');
    zlabel('Z', 'FontSize',14, 'FontWeight', 'bold');
    end


end