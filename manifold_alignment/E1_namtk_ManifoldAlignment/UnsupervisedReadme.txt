~~~~~~~~~~~~~~~~~~ Functions used in Alignment ~~~~~~~~~~~~~~~~~

%Feature-level Manifold Projections Two domain
wmapGeneralTwo.m

%Unsupervised Alignment
%The code is the same as feature-level manifold projections.
%The only difference is how to create the correspondence matrix.
generateWeight3.m:  %This one is used to generate the weight matrix.
%It calls computeOptimalMatch.m and decompose3.m 

%knn Search method
knnsearch.m


%Create All connected Graph from the given data
createAllConnectedGraph.m

%Create knn graph from the given data
createKnnGraph.m

%Compute Euclidean distance
L2_distance.m




~~~~~~~~~~~~~~~~~~~Protein Domain ~~~~~~~~~~~~~~~

%Read protein data into memory
'1G7O-1.pdb.3d.real', '1G7O-10.pdb.3d.real', '1G7O-21.pdb.3d.real'==>
!getData.m


%align the input proteins.

%Feature-level manifold projections, 
%Unsupervised Feature-level manifold projections.
cmp2All.m



% Visualization of the alignment results of cmp2All.m
Visualize2.m

