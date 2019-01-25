function [q2best] = Find_Rotation_and_Seed_unique(X1,X2)
% Resample the curves to have N points
      N =70;
    X1 = ReSampleCurve(X1,N);
    X2 = ReSampleCurve(X2,N);

%Center curves, not really needed but good for display purposes
    X1 = X1 - repmat(mean(X1')',1,size(X1,2));
    X2 = X2 - repmat(mean(X2')',1,size(X2,2));   


% Form the q function for representing curves and find best rotation
    [q1] = curve_to_q(X1);
    [q2] = curve_to_q(X2);
    A = q1*q2';lam=0;
    Disp_geodesic_between_the_curves = 1;
    Disp_registration_between_curves = 1;
    Disp_optimal_reparameterization = 1;
    [U,S,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        if (size(X1,1)==2)
            Ot = U*([V(:,1) -V(:,2)])';
        else
            Ot = U*([V(:,1) V(:,2) -V(:,3)])';
        end
    end
    X2 = Ot*X2;
    q2 = Ot*q2;

% Applying optimal re-parameterization to the second curve
    [gam] = DynamicProgrammingQ(q1/sqrt(InnerProd_Q(q1,q1)),q2/sqrt(InnerProd_Q(q2,q2)),lam,0);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    X2n = Group_Action_by_Gamma_Coord(X2,gamI);
    q2n = curve_to_q(X2n);
    if Disp_optimal_reparameterization
        figure(100); clf;
        plot(linspace(0,1,N),gamI,'LineWidth',2)
    end
    
% Find optimal rotation
    A = q1*q2n';
    [U,S,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        if (size(X1,1)==2)
            Ot = U*([V(:,1) -V(:,2)])';
        else
            Ot = U*([V(:,1) V(:,2) -V(:,3)])';
        end
    end
    X2n = Ot*X2n;
    q2best = Ot*q2n;