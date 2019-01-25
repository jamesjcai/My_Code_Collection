% This function projects the tangent vector f in the
% Tangent space of {\cal C} at q
function fnew = Project_Tangent(f,q)

[n,T] = size(q);
% Project w in T_q({\cal B}), ie the unit sphere
w = f - InnerProd_Q(f,q)*q;
e = eye(n);
% Form the basis for the Normal space of {\cal A}
g = Form_Basis_Normal_A(q);

% Refer to the function Gram_Schmidt for the parameters
Evorth = Gram_Schmidt(g,'InnerProd_Q');
% Unpack Evorth structure
for i = 1:n
    Ev(:,:,i) = Evorth{i};
end

for i = 1:n
    fnew = w - InnerProd_Q(w,Ev(:,:,i))*Ev(:,:,i);
end

return;