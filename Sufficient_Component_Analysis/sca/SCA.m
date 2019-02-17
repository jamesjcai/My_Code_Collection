%Sufficient Component Analysis
function [W_opt, MIh,F,L,Rx] = SCA(X,Y,r)

y_type = 0;
iniflag = 0;

d = size(X,1);
n = size(X,2);

%Subtract mean
Y = Y - mean(Y,2)*ones(1,size(Y,2));
Y = Y + randn(size(Y))*eps;
X = X - mean(X,2)*ones(1,size(X,2));

%Initialization
W = eye(d);
W_opt = W;
MIh_old = -realmax;

if iniflag == 1
    maxite = 1;
else
    maxite = 200;
end

for ite = 1:maxite
    [MIht, alphah,Ky,Kx,U, sigma_chosen]=LSMI(W_opt*W_opt'*X,Y,y_type);
    
    sigma = sigma_chosen;
    MIh(ite) = MIht;
    
    %Computing D matrix
    Kxy = (Kx>0).*Ky;
    Kxy = Kxy.*(alphah*ones(1,n));
    Kxys = Kxy/2/sigma^2;
    vol = sum(sum(Kxy));
    
    tmp = sum(Kxys,1);
    Tmp = ones(d,1)*tmp;
    Rx1 = (X.*Tmp)*X';
    
    tmp = sum(Kxys,2)';
    Tmp = ones(d,1)*tmp;
    Rx2 = (U.*Tmp)*U';
    
    Tmp = U*Kxys;
    Rx12 = X*Tmp';
    
    Tmp = X*Kxys';
    Rx21 = U*Tmp';
    
    Rx = vol*eye(d)/r/n + (- Rx1 - Rx2 + Rx12 + Rx12')/n;
    %Rx = (- Rx1 - Rx2 + Rx12 + Rx12')/n;
    [F L] = eig(Rx);
    %[F L] = eig(Rx12+Rx12', Rx1 + Rx2);
    dd = diag(L);
    [val,ind] = sort(-dd);
    F = F(:,ind);
    L = L(ind,ind);
    
    [Wt L] = qr(F);
   
    if ite == 1
        W = Wt(:,1:r);
        W_opt = W;
        MIh_old = MIh(ite);
    else 
        Wt = Wt(:,1:r);
          
        if MIh(ite) - MIh_old < 0 | abs(MIh(ite) - MIh_old) < 0.000001
            return
        else
            MIh_old = MIh(ite);
            W_opt = Wt;%F(:,1:r);
        end
    end
    MIh
end


