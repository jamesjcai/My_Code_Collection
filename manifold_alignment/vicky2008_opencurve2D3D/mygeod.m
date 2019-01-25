function [dist,X2n,q2n,X1,q1,gam]=mygeod(X1,X2)
%input: two curves X1 and X2 as 2xn or 3xn (2D vs. 3D curve)
%output: dist=distance, also will display when you run the program; X2n:
%optimally registered curve X2; q2n: same as X2n except in q-fun form; X1:
%normalized curve X1; q1: same as X1 except in q-fun form;

% Load some parameters, no need to change this
    lam = 0;

% What displays you want to see, set to 1 if you want to see figures
    Disp_geodesic_between_the_curves = 1;
    Disp_registration_between_curves = 1;
    Disp_optimal_reparameterization = 1;
    
% Resample the curves to have N points
    N = 200;
    X1 = ReSampleCurve(X1,N);
    X2 = ReSampleCurve(X2,N);

%Center curves, not really needed but good for display purposes
    X1 = X1 - repmat(mean(X1')',1,size(X1,2));
    X2 = X2 - repmat(mean(X2')',1,size(X2,2));    
    
% Form the q function for representing curves and find best rotation
    [q1] = curve_to_q(X1);
    [q2] = curve_to_q(X2);
    A = q1*q2';
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
%      [gam] = DynamicProgrammingQ(q1/sqrt(InnerProd_Q(q1,q1)),q2/sqrt(InnerProd_Q(q2,q2)),lam,0);
     [gam] = DynamicProgrammingQ(q1,q2,lam,0);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    X2n = Group_Action_by_Gamma_Coord(X2,gamI);
    q2n = curve_to_q(X2n);
    if Disp_optimal_reparameterization
        figure(100); clf;
        plot(linspace(0,1,N),gamI,'LineWidth',2)
    end
    save gam;
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
    q2n = Ot*q2n;

% Forming geodesic between the registered curves
N = size(X1,2);
dist = acos(sum(sum(q1.*q2n))/N);
sprintf('The distance between the two curves is %0.3f',dist)
    if(Disp_geodesic_between_the_curves)
         if (size(X1,1)==2)
             for t=1:7
                 s = dist*(t-1)/6;
                 PsiQ(:,:,t) = (sin(dist - s)*q1 + sin(s)*q2n)/sin(dist);
                 PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
             end
             figure(4); clf; axis equal; hold on;
             for t=1:7
             	z = plot(0.5*t + PsiX(1,:,t), PsiX(2,:,t),'r-'); 
                 set(z,'LineWidth',[2],'color',[(t-1)/6 (t-1)/12 0]);
             end
             axis off;       
         else
            for t=1:7
                s=(t-1)/6;
%                 s = dist*(t-1)/6;
                PsiQ(:,:,t) = (1-s)*q1 + s*q2n;
%                 PsiQ(:,:,t) = (sin(dist - s)*q1 + sin(s)*q2n)/sin(dist);
                PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
            end
            figure(4); clf;
             axis equal;
            hold on;
            for t=1:7
                z = plot3(0.5*t + PsiX(1,:,t), PsiX(2,:,t), PsiX(3,:,t),'r-'); 
                set(z,'LineWidth',[2],'color',[(t-1)/6 (t-1)/12 0]);
            end
             axis equal;
               axis off;
        end
    end
    
% Displaying the correspondence
if(Disp_registration_between_curves)
    X2n=PsiX(:,:,end);
    X1=PsiX(:,:,1);
    if (size(X1,1)==2)
        figure(3); clf;
        z = plot(X1(1,:), X1(2,:),'r');
        set(z,'LineWidth',[2]);
         axis off;
        hold on;
        z = plot(X2n(1,:), 0.15+X2n(2,:),'b-+');
        set(z,'LineWidth',[3]);
        N = size(X1,2);
        for j=1:N/15
            i = j*15;
            plot([X1(1,i) X2n(1,i)],[X1(2,i) 0.15+X2n(2,i)], 'k');
        end
    else
        figure(3); clf;
        z = plot3(X1(1,:), X1(2,:), X1(3,:),'r');
        set(z,'LineWidth',[2]);
        axis off;
        hold on;
        z = plot3(X2n(1,:), 0.15+X2n(2,:), X2n(3,:),'b-+');
        set(z,'LineWidth',[3]);
        N = size(X1,2);
        for j=1:N/15
            i = j*15;
            plot3([X1(1,i) X2n(1,i)],[X1(2,i) 0.15+X2n(2,i)], [X1(3,i) X2n(3,i)], 'k');
        end
    end
end
% %X1=Data(:,:,1);X2=Data(:,:,10);[dist,X2n,q2n,X1,q1]=mygeod(X1,X2)
% % figure;X1=Data(:,:,1);X22=Data(:,:,2);X3=Data(:,:,3);X4=Data(:,:,4);X5=Data(:,:,5);
% % X6=Data(:,:,6);X7=Data(:,:,7);X8=Data(:,:,8);X9=Data(:,:,9);X2=Data(:,:,10);
% % 
% % subplot(1,10,1);plot(X1(1,:),X1(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,2);plot(X22(1,:),X22(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,3);plot(X3(1,:),X3(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,4);plot(X4(1,:),X4(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,5);plot(X5(1,:),X5(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,6);plot(X6(1,:),X6(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,7);plot(X7(1,:),X7(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,8);plot(X8(1,:),X8(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,9);plot(X9(1,:),X9(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% % subplot(1,10,10);plot(X2(1,:),X2(2,:),'b','LineWidth',2);
% % axis off;axis equal;
% 
% 
% 
% %A=C(:,:,1:5:50);Data1=A(:,1:70,:);
% 
% %A=C(:,:,1:8:8*10);Data1=A(:,1:70,:);
% % for i=1:10
% %     Data(:,:,i)=[Data1(:,:,i);1:i:70*i];
% % end
% %X1=Data(:,:,1);X2=Data(:,:,10);[dist,X2n,q2n,X1,q1]=mygeod(X1,X2)
% % figure;X1=Data(:,:,1);X2=Data(:,:,2);X3=Data(:,:,3);X4=Data(:,:,4);X5=Data(:,:,5); X6=Data(:,:,6);X7=Data(:,:,7);X8=Data(:,:,8);X9=Data(:,:,9);X10=Data(:,:,10);
%    
% %       subplot(1,10,1);plot3(X1(1,:),X1(2,:),X1(3,:),'b','LineWidth',2); axis off;axis equal;
% %      subplot(1,10,2);plot3(X2(1,:),X2(2,:),X2(3,:),'b','LineWidth',2);axis off;axis equal;
% %       subplot(1,10,3);plot3(X3(1,:),X3(2,:),X3(3,:),'b','LineWidth',2); axis off;axis equal;
% %       subplot(1,10,4);plot3(X4(1,:),X4(2,:),X4(3,:),'b','LineWidth',2); axis off;axis equal;
% %       subplot(1,10,5);plot3(X5(1,:),X5(2,:),X5(3,:),'b','LineWidth',2); axis off;axis equal;
% %      subplot(1,10,6);plot3(X6(1,:),X6(2,:),X6(3,:),'b','LineWidth',2); axis off;axis equal;
% %       subplot(1,10,7);plot3(X7(1,:),X7(2,:),X7(3,:),'b','LineWidth',2); axis off;axis equal;
% %      subplot(1,10,8);plot3(X8(1,:),X8(2,:),X8(3,:),'b','LineWidth',2); axis off;axis equal;
% %       subplot(1,10,9);plot3(X9(1,:),X9(2,:),X9(3,:),'b','LineWidth',2); axis off;axis equal;
% %       subplot(1,10,10);plot3(X10(1,:),X10(2,:),X10(3,:),'b','LineWidth',2);axis off;axis equal;
