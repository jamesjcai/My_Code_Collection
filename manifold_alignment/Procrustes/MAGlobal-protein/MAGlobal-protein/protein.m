% A simple example to show how to do manifold alignment preserving pairwise
% distance
function protein
    %total number of instances.
    N=215;
    %how many corresponding pairs to use= p*N 
    p=0.10; 
    %sample rate
    step=floor(1/p);

    %Retrieve data
    [X1, X2, X3]=getData;
    epsilon=1e-8;

    %Building the pairwise distance matrices
    D1=L2_distance(X1, X1, 0);
    D2=L2_distance(X2, X2, 0);
    D3=L2_distance(X3, X3, 0);
 
    %Comment the following 3 lines to do alignment with the original
    %distance matrices
    D1=dijk(D1); 
    D2=dijk(D2);  
    D3=dijk(D3);

    %Create correspondence matrices.
    W12=sparse(N, N);
    W13=sparse(N, N);
    W23=sparse(N, N);
    for i=1:step:N;
        W12(i,i)=1;
        W13(i,i)=1;
        W23(i,i)=1;
    end
    
    
    %Feature-level Manifold alignment preserving global geometry     
    [f1, f2, f3]=  dmapGeneralThree(X1, X2, X3, D1, D2, D3, W12, W13, W23, epsilon);
    visualize3(X1, X2, X3, f1(:,1:3)'*X1, f2(:,1:3)'*X2, f3(:,1:3)'*X3, 4);
    set(gcf,'Name', 'Feature-Level Manifold Alignment preserving Global Geometry'); pause(1); 
    
    %Instance-level Manifold alignment preserving global geometry
    [g1, g2, g3]=  dmapGeneralThreeInstance(X1, X2, X3, D1, D2, D3, W12, W13, W23, epsilon);
    visualize3(X1, X2, X3, g1(:,1:3)', g2(:,1:3)', g3(:,1:3)', 4);
    set(gcf,'Name', 'Instance-Level Manifold Alignment preserving Global Geometry'); pause(1); 
end