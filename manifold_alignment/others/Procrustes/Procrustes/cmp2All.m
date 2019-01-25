function cmp2All

    %total number of instances.
    N=215;
    %how many corresponding pairs to use= p*N 
    p=0.1; 
    %sample rate.
    step=floor(1/p);
        
    %Retrieve data
    [X1, X2, X3]=getData;


 
    %Procrustes Alignment
    j=0; clear X1train X2train;
    for i=1:step:N;
        j=j+1;
        X1train(:,j)=X1(:,i);
        X2train(:,j)=X2(:,i);
    end
    [Q, k, X1_mean, X2_mean] = Procrustes (X1train, X2train);
    visualize2(X1, X2, X1-repmat(X1_mean',1, N), (k*(X2-repmat(X2_mean', 1, N))'*Q)', 2);
   
end