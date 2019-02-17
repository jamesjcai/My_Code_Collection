clear all;
close all;

dataset=1;
n =1000;

seed = 1;

rand('state',seed);
randn('state',seed);

switch dataset
    case 1 
        X = randn(4,n);
        E = randn(1,n);
        sg = 0.1;

        X(1,:) = X(1,:);
        X1 = X(1,:);
        X2 = X(2,:);
        Y = (X1.^2 + X2)./(0.5 + (X2 + 1.5).^2) + (1 + X2).^2+ sg*E;
        W0 = [1 0 0 0; 0 1 0 0]';
        reduce_dim=2;
    case 2 
        E = randn(1,n);
        sg = 0.5;
        
        reduce_dim = 1;
        X=(rand(4,n)*2-1)*1;
        Y=X(2,:) + sg*E;
        W0 = [0 1 0 0]';
    case 3
        E = randn(1,n);
        sg = 0.1;
        
        reduce_dim = 1;
        X=(rand(5,n)*2-1)*1;
        Y=1/2*(X(2,:)).^3.*E;
        W0 = [0 1 0 0 0]';
    case 4 
        m = 30;
        E = randn(1,n);
        sg = 0.1;
        
        reduce_dim = 1;
        X=(randn(m,n));
        
        Y = (X(3,:)).^2 + sg*E;
        W0 = [0 0 1 zeros(1,m-3)]';
end

d = size(X,1);

y_type = 0;
iniflag = 0;
%
t = cputime;
[Wf, MIh,F,L] = SCA(X,Y,reduce_dim);
timescadm(seed) = cputime - t;
err2(seed) = norm(W0*W0' - Wf*Wf','fro')/sqrt(2*reduce_dim)

plot(MIh)
