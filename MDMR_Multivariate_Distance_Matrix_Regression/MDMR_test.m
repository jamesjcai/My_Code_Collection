function [Pperm,Fstat,var_explained,D]=MDMR_test(X,Y)

% X varibles being tested for regression, row - samples, column - age, etc...
% Y gene expression matrix, row - samples, column - genes

	n=size(Y,1);
    
	X=[ones(n,1),X];
    H=X*inv(X'*X)*X';

    %mdl=fitglm(X,Y(:,1));
    %H=mdl.Diagnostics.HatMatrix;
    
	[D]=i_dmatrix_eu(Y);
	A=-0.5*(D.^2);
	
    
    I=eye(n);
    S=I-(1/n).*(ones(n,1)*ones(n,1)');  % Centering matrix
	G=S*A*S;

	Fstat=trace(H*G*H)/trace((I-H)*G*(I-H));

	% an estimation of the proportion of variation within the 
	% matrix that is explained by a particular set of 
	% M predictor variables can be calculated by dividing tr(HGH) 
	% (i.e., the sum of the diagonal elements of a matrix) by tr(G).
	if nargout>2
	      var_explained=trace(H*G*H)/trace(G);
	end


	M=1000;
	Fr=zeros(M,1);
	[n,m]=size(G);

	parfor k=1:M
	    Gr=G(randperm(n),randperm(m));    
	    Fr(k)=trace(H*Gr*H)/trace((I-H)*Gr*(I-H));
	end
	Pperm=sum(Fr>=Fstat)./M;
end




function [D]=i_dmatrix_md(Y)

        mcd=mcdcov(Y,'plots',0);
        
        n=size(Y,1);
        D=zeros(n);
        for i=1:n-1
            for j=i+1:n
            % MDc(i)=(c1(i,:)-mcd.center)*inv(mcd.cov)*(c1(i,:)-mcd.center)';
            D(i,j)=(Y(i,:)-Y(j,:))*inv(mcd.cov)*(Y(i,:)-Y(j,:))';
            end
        end
        D=D+D';
        D=sqrt(D);
end


function [D]=i_dmatrix_md2(Y)

        mcd=mcdcov(Y,'plots',0);
        
        n=size(Y,1);
        D=zeros(n);
        for i=1:n-1
            for j=i+1:n
            MDi=(Y(i,:)-mcd.center)*inv(mcd.cov)*(Y(i,:)-mcd.center)';
            MDj=(Y(j,:)-mcd.center)*inv(mcd.cov)*(Y(j,:)-mcd.center)';
            D(i,j)=abs(sqrt(MDi)-sqrt(MDj));
            end
        end
        D=D+D';
end



function [D]=i_dmatrix_co(Y)
        n=size(Y,1);
        D=zeros(n);
        for i=1:n-1
            for j=i+1:n            
            r=abs(corr(Y(i,:)',Y(j,:)'));
            D(i,j)=sqrt(2*(1-r));
            end
        end
        D=D+D';
end


function [D]=i_dmatrix_eu(Y)
        n=size(Y,1);
        D=zeros(n);
        for i=1:n-1
            for j=i+1:n 
                D(i,j) = norm(Y(i,:) - Y(j,:));
                % D(i,j) = sqrt((Y(,:)-Y(2,:))*(Y(1,:)-Y(2,:))');
                % D(i,j) = sqrt(sum((Y(1,:) - Y(2,:)).^2));
            end
        end
        D=D+D';
end

