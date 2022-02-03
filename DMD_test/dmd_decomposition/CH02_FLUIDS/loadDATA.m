nx = 199;  % number of grid points in y-direction
ny = 449;  % number of grid points in x-direction

% create space for 150 snapshots
VORTALL = zeros(nx*ny,150); % vorticity

% extract data from 150 snapshots files
for count=1:150
    num = count*10; % load every 10th file    
    % load file
    fname = ['ibpm',num2str(num,'%05d'),'.plt'];       
    [X,Y,U,V,VORT] = loadIBPM(fname,nx,ny);
    VORTALL(:,count) = reshape(VORT,nx*ny,1); 
end