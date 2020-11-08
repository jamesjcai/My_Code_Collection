function pc_projection=pc_project_multidim(X,Xinit,kernel_sigma,targetdim)

point=Xinit;
N=size(X,2);
dim=size(X,1);
threshold=1e-3;
for ind=1:size(point,2)
    flag=0;
    [p gra H SI] = pgh(point(:,ind),X',N,dim,kernel_sigma);
    point_pc1=point(:,ind);
    [V D]= eig(SI);
    %[dummy index]=max(abs(diag(D)));
    [dummymax indexsorted]=sort(diag(abs(D)),'descend');
    ConstrainedSpace = V(:,indexsorted(1:end-targetdim));
    for tim=1:size(ConstrainedSpace,2),
        direction = ConstrainedSpace(:,tim);
        direction = direction * sign(direction'*gra);
        ConstrainedSpace(:,tim)=direction;
    end
    if abs(gra'*H*gra/(norm(gra'*H)*norm(gra'*H)))<0.01,
        flag=1;
        pc(ind,:)= point_pc1;
    end
        if ~flag,
            for a=1:20,
                G = kernel_matrix(point_pc1,X',dim,1,N,kernel_sigma);
                num1 = sum(repmat(G,dim,1).*X,2);
                den1 = sum(G,2);
                if(direction'*gra<0),keyboard;end
                [p1 gra1 H SI] = pgh(point_pc1,X',N,dim,kernel_sigma);
                point_pc1_old=point_pc1;
                for c=1:size(ConstrainedSpace,2),
                    direction=ConstrainedSpace(:,c);
                    point_pc1 = point_pc1 + direction * (direction'*(num1/den1-point_pc1));
                end

                if abs(point_pc1_old-point_pc1)<threshold,break;end
            end
            pc(ind,:)= point_pc1;
        end
end
pc_projection=pc;

%Xdeflated=X-repmat(direction,1,N).*repmat((direction'*X),dim,1);

function K=kernel_matrix(data1,data2,n,N1,N2,kernel_sigma)
for j=1:N1,
    C=1./((2*pi)^(n/2)*kernel_sigma.^n);
    dx=(repmat(data1(:,j),1,N2)-data2')./repmat(kernel_sigma,n,N2);
    K(j,:)=C.*exp(-0.5*sum(dx.^2,1));
end