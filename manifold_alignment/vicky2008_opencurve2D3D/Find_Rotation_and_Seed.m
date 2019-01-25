function [q2best,Rbest] = Find_Rotation_and_Seed(q1,q2,reparamFlag)

[n,T] = size(q1);

scl = 4;
minE = 1000;
for ctr = 0:floor((T)/scl)
    q2n = ShiftF(q2,scl*ctr);
    [q2n,R] = Find_Best_Rotation(q1,q2n);
%     function [q2new,R] = Find_Best_Rotation(q1,q2)
% % Assumes the starting points are fixed
% 
% [n,T] = size(q1);
% A = q1*q2';
% [U,S,V] = svd(A);
% if det(A) > 0
%     S = eye(n);
% else
%     S = eye(n);
%     S(:,end) = -S(:,end);
% end
% R = U*S*V';
% q2new = R*q2;
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1,q2n,0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            p2n = q_to_curve(q2n);
            p2new = Group_Action_by_Gamma_Coord(p2n,gamI);
            q2new = curve_to_q(p2new);
%             q2new = ProjectC(q2new);
        else
            q2new = q2n;
        end
        
    else
        q2new  = q2n;
    end
    Ec = acos(InnerProd_Q(q1,q2new));
    if Ec < minE
        Rbest=R;
        q2best  = q2new;
        minE = Ec;
    end
end

return;