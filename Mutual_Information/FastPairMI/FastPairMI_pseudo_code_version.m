function MI  = FastPairMI_pseudo_code_version(data,h)
% e.g. : MI  = FastPairMI_pseudo_code_version(data,h)
% data      : the input data, rows correspond to genes
%             columns correspond to arrays (samples) 
% h         : the std of the Gaussian kernel for density estimation 
%             if the data is normalized to 0-mean 1-var, h is recommended to be 0.3

MI = zeros(size(data,1));
h_square = h^2;

M = size(data,2); % M samples
N = size(data,1); % N genes

for i=1:M                                                                       % 1
    Dist=[];                                                                    % 2
    SumMargin = zeros(N,1);                                                     % 3  lines 3 and 4 of the pseudo code in the published manuscript should be switched
    for j=1:M                                                                   % 4
        for k=1:N                                                               % 5
            Dist(k,j) = KernelDistance(data(k,j),data(k,i),h_square);           % 6
            SumMargin(k) = SumMargin(k) + Dist(k,j);                            % 7
        end                                                                     % 8 
    end                                                                         % 9
    %                                                                           % 10
    for k=1:N                                                                   % 11  This line and the following one are a little bit different from the psuedo-code
        for l = k:N                                                             % 12  The psuedo-code does not compute the MI between a gene and it self. With the modification here, this code does that
            SumJoint = 0;                                                       % 13
            for j=1:M                                                           % 14
                SumJoint = SumJoint + Dist(k,j)*Dist(l,j);                      % 15
            end                                                                 % 16
            MI(k,l) = MI(k,l) + log(SumJoint/SumMargin(k)/SumMargin(l));        % 17
        end                                                                     % 18
    end                                                                         % 19
end                                                                             % 20

 
MI = MI + triu(MI,1)';  % the above only fills half of the MI matrix, this line fill the symmetric half
MI = MI/M+log(M);       % this line deals with the two M's in equation (4) in the manuscript

return



%%%% the function below corresponds to the kernel distance in line 6 of the
%%%% pseudo-code
function [d] = KernelDistance(a,b,h_square)                                   
d = exp(-1/2/h_square*((a-b)^2));
