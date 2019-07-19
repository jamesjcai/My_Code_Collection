%--------------------------------------------------------------------------
% mider: Mutual Information Distance & Entropy Reduction
%    
%   [Output] = mider(x,options) computes several entropic measures from the
%   data 'x' using the provided 'options', and stores the results in a 
%   structure named 'Output'.
% 
%   x       = p*n array of input data (p data points; n variables)
%   options = structure with the following fields:
%             q          = entropic parameter
%             fraction   = fraction of bins with at least 5 points
%             taumax     = maximum time lag considered
%             ert_crit   = number of entropy reduction rounds
%             threshold  = minimum entropy reduction for a connection 
%             MItype     = type of MI used to calculate distance map
%             plotMI     = plot MI arrays or not
%   Output  = structure with the following fields:
%             Y          = coordinates of the points from multidim. scaling
%             con_array  = n*n array of connections between variables
%             MI         = n*n*(nlags+1) array of mutual information
%             MIl        = MI normalized as in Linfoot (1957)
%             MIm        = MI normalized as in Michaels et al (1998) 
%             MIs        = MI normalized as in Studholme et al (1999) 
%             H1         = n-vector of entropies
%             H2         = n*n*(nlags+1) array of joint entropy, 2 variables
%             H3         = n*n*n array of joint entropy, 3 variables
%             H4         = n*n*n*n array of joint entropy, 4 variables
%             MI3        = n*n*n array of three-way mutual information
%             dist       = n*n*(nlags+1) array of distance between variables
%             taumin     = n*n array of the time lags that minimize dist
%             cond_entr2 = n*n array of conditional entropies, H(Xm|Xn)
%             cond_entr3 = n*n*n array of cond. entropies of 3 variables
%             cond_entr4 = n*n*n*n array of cond. entropies of 4 variables
%             adaptThres = adaptive threshold value
%             T          = n*n array of transfer entropies
%
%	Dependencies: mider.m calls estimateH2.m, estimateH3.m, estimateH4.m
%   mider.m is called by runMIDER.m
% 
%   Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
%   Created: 07/10/2011, last modified: 25/11/2014
%--------------------------------------------------------------------------
% Copyright (C) 2014  Alejandro Fernández Villaverde
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function [Output] = mider(x,options)

fprintf(1,'\n (... using the Matlab/Octave implementation)');

%--------------------------------------------------------------------------
% Read input and options:

ntotal     = size(x,2); 
q          = options.q; 
fraction   = options.fraction; 
taumax     = options.taumax;
ert_crit   = options.ert_crit; 
threshold  = options.threshold;  
MItype     = options.MItype;  


%--------------------------------------------------------------------------
% Initialize arrays:

MI  = zeros(ntotal,ntotal,taumax+1);      % mutual information
H2  = zeros(ntotal,ntotal,taumax+1);      % joint entropy, 2 variables


%--------------------------------------------------------------------------
% Joint entropy (H2) and mutual information (MI):

tic 
for tau = 0:taumax
    for i=1:ntotal
        for j=1:ntotal
            pb = 5;
            [MI(i,j,tau+1),fracn,H2(i,j,tau+1)] = ...
                estimateH2(x(1:end-tau,i),x(tau+1:end,j),pb,q); 
            while fracn < fraction; 
                pb = pb+1;
                [MI(i,j,tau+1),fracn,H2(i,j,tau+1)] = ...
                    estimateH2(x(1:end-tau,i),x(tau+1:end,j),pb,q);  
            end
        end
    end
end


%--------------------------------------------------------------------------
% Entropy of each variable (H1):

H1 = diag(H2(:,:,1));


%--------------------------------------------------------------------------
% Normalizations of the mutual information:

MIl = zeros(ntotal,ntotal,taumax+1); % Linfoot (1957)
MIm = zeros(ntotal,ntotal,taumax+1); % Michaels et al (1998)
MIs = zeros(ntotal,ntotal,taumax+1); % Studholme et al (1999)

for i=1:ntotal
    for j=1:ntotal
        MIl(i,j,:) = sqrt( 1 - (exp(1)).^(-2*MI(i,j,:)) );
        MIm(i,j,:) = MI(i,j,:) / max(H1(i),H1(j));
        MIs(i,j,:) = ( H1(i) + H1(j) ) ./ H2(i,j,:) ;
    end
end


%--------------------------------------------------------------------------
% EMC distance:

if strcmp(MItype,'MI'),          EMC_dist_tau = exp(-MI);  end
if strcmp(MItype,'MIlinfoot'),   EMC_dist_tau = exp(-MIl); end
if strcmp(MItype,'MImichaels'),  EMC_dist_tau = exp(-MIm); end
if strcmp(MItype,'MIstudholme'), EMC_dist_tau = exp(-MIs); end

[dist,taumin] = min(EMC_dist_tau,[],3); 


%--------------------------------------------------------------------------
% Conditional entropies, H(Xm|Xn):

cond_entr2     = zeros(ntotal,ntotal);
cond_entr2_A   = zeros(ntotal,ntotal); 
cond_entr2_B   = zeros(ntotal,ntotal); 
entropy_error    = zeros(ntotal,ntotal);
rel_entr_reduc_2 = zeros(ntotal,ntotal);
for m=1:ntotal
    for n=1:ntotal
        cond_entr2_A(m,n) = H1(m) - MI(m,n,1); % 1st way of calculating H(X|Y) 
        cond_entr2_B(m,n) = H2(n,m,1) - H1(n); % 2nd way of calculating H(X|Y)
        cond_entr2(m,n)   = cond_entr2_A(m,n); % choice of H(X|Y)
        entropy_error(m,n)= abs(cond_entr2_A(m,n) - cond_entr2_B(m,n));
        rel_entr_reduc_2(m,n) = MI(m,n,1)/H1(m);
    end
end

fprintf(1,'\n H1(Xm), H2(Xm,Xn), H(Xm|Xn), MI(Xm,Xn), '); 
fprintf(1,'and distance calculated in %d seconds.',toc);   

  
% -------------------------------------------------------------------------
% ERT, 1st round:

if ert_crit >= 1 

tic

% For each species, determine which of the remaining species yields the
% smallest conditional entropy (that is the most likely connection): 
[~, spec_ind_2] = min(cond_entr2 + 1e3*eye(ntotal),[],2); % (add 1e3 to avoid minima)
ERT_array_2              = [(1:ntotal)', spec_ind_2];

ERT_time = toc;


%--------------------------------------------------------------------------
% ERT, 2nd round:

if ert_crit >= 2
    
% Obtain H3 first:
tic 
H3  = zeros(ntotal,ntotal,ntotal);
for i=1:ntotal
    for j=1:ntotal
        for k=1:ntotal
            pb = 5;
            [H3(i,j,k),fracn] = estimateH3(x(:,i),x(:,j),x(:,k),pb,q);
            while fracn < fraction;
                pb = pb+1;
                [H3(i,j,k),fracn] = estimateH3(x(:,i),x(:,j),x(:,k),pb,q);
            end
        end
    end
end

% Now MI3 (not a part of ERT) can be obtained from H1, H2, and H3:
MI3 = zeros(ntotal,ntotal,ntotal);
for i=1:ntotal
    for j=1:ntotal
        for k=1:ntotal
            MI3(i,j,k) = -H1(i)-H1(j)+2*H2(i,j)+H2(i,k)+H2(j,k)-2*H3(i,j,k);
        end
    end
end
fprintf(1,'\n H3, MI3 calculated in %d seconds.',toc); 
ERT_time = ERT_time + toc;

% Conditional entropy of triplets of variables
% (i.e. entropy of one variable conditional to the entropy of a pair): 
% H(Xm|(Xn,Xp)) = cond_entr3(m,n,p)
tic
cond_entr3   = zeros(ntotal,ntotal,ntotal); 
rel_entr_reduc_3 = zeros(ntotal,ntotal,ntotal);
for m=1:ntotal
    for n=1:ntotal
        for p=1:ntotal
            cond_entr3(m,n,p) = H3(n,p,m) - H2(n,p);
            rel_entr_reduc_3(m,n,p) =...
                ( cond_entr2(m,n) - cond_entr3(m,n,p) ) / H1(m);
        end
    end
end

% Find, for each variable, which of the remaining pairs of variables yields 
% the smallest conditional entropy (i.e. the 2 most likely connections):
diag_ones_3 = zeros(ntotal,ntotal,ntotal);
for i=1:ntotal
    diag_ones_3(i,i,:) = 1;
    diag_ones_3(i,:,i) = 1;
    diag_ones_3(:,i,i) = 1;
end
cond_entr3m = cond_entr3 + 1e3*diag_ones_3; % (add 1e3 to avoid minima)
min_diff_3 = zeros(ntotal,1);
spec_ind_3 = zeros(ntotal,1);
for i=1:ntotal
    [min_diff_3(i), spec_ind_3(i)] = ...
        min(cond_entr3m(ERT_array_2(i,1),ERT_array_2(i,2),:));
end
ERT_array_3 = [(1:ntotal)', spec_ind_2, spec_ind_3];

ERT_time = ERT_time + toc;


%--------------------------------------------------------------------------
% ERT, 3rd round:

if ert_crit >= 3
    
% Obtain H4 first:
tic 
H4  = zeros(ntotal,ntotal,ntotal,ntotal); 
for i=1:ntotal
    for j=1:ntotal
        for k=1:ntotal
            for l=1:ntotal
                pb = 5;
                [H4(i,j,k,l),fracn] = ...
                    estimateH4(x(:,i),x(:,j),x(:,k),x(:,l),pb,q);
                while fracn < fraction;
                    pb = pb+1;
                    [H4(i,j,k,l),fracn] = ...
                        estimateH4(x(:,i),x(:,j),x(:,k),x(:,l),pb,q);
                end
            end
        end
    end
end
fprintf(1,'\n H4 calculated in %d seconds.',toc); 

ERT_time = ERT_time + toc;

% Entropy of one variable conditional to the entropy of a triple:
% H(Xm|(Xn,Xp,Xr)) = cond_entr4(m,n,p,r)
tic
cond_entr4 = zeros(ntotal,ntotal,ntotal,ntotal); 
rel_entr_reduc_4 = zeros(ntotal,ntotal,ntotal,ntotal); 
for m=1:ntotal
    for n=1:ntotal
        for p=1:ntotal
            for r=1:ntotal
                cond_entr4(m,n,p,r) = H4(n,p,r,m) - H3(n,p,r);
                rel_entr_reduc_4(m,n,p,r) = ...
                ( cond_entr3(m,n,p) - cond_entr4(m,n,p,r) ) / H1(m);
            end
        end
    end
end

% Find, for each variable, which of the remaining triplets of variables
% yields the smallest conditional entropy (the 3 most likely connections):
diag_ones_4 = zeros(ntotal,ntotal,ntotal,ntotal);
for i=1:ntotal
    diag_ones_4(i,i,:,:) = 1;
    diag_ones_4(i,:,i,:) = 1;
    diag_ones_4(i,:,:,i) = 1;
    diag_ones_4(:,i,i,:) = 1;
    diag_ones_4(:,i,:,i) = 1;
    diag_ones_4(:,:,i,i) = 1;
end
cond_entr4m = cond_entr4 + 1e3*diag_ones_4; % (add 1e3 to avoid minima)
min_diff_4 = zeros(ntotal,1);
spec_ind_4 = zeros(ntotal,1);
for i=1:ntotal
    [min_diff_4(i), spec_ind_4(i)] =...
        min(cond_entr4m(ERT_array_3(i,1),ERT_array_3(i,2),ERT_array_3(i,3),:));
end
ERT_array_4 = [(1:ntotal)', spec_ind_2, spec_ind_3, spec_ind_4];

ERT_time = ERT_time + toc;


%--------------------------------------------------------------------------
% ERT results messages:

fprintf(1,'\n 3 entropy reduction rounds have been carried out in %d seconds.',ERT_time);  

else
    fprintf(1,'\n The 1st and 2nd entropy reduction rounds have been carried out');
    fprintf(1,' in %d seconds.',ERT_time);
    fprintf(1,'\n To carry out a 3rd round, set ert_crit = 3.');
end

else
    fprintf(1,'\n The 1st entropy reduction round has been carried out');
    fprintf(1,' in %d seconds.',ERT_time); 
    fprintf(1,'\n To carry out a 2nd round, set ert_crit = 2 ');
    fprintf(1,'(alternatively, set ert_crit = 3 for 2nd & 3rd rounds).');
end

else
    fprintf(1,'\n Entropy reduction (ERT) has not been carried out.');
    fprintf(1,'\n The diagram will be plotted showing only the entropic distances.');
    fprintf(1,'\n Interactions can be estimated from the entropic distances.');
    fprintf(1,'\n To carry out ERT, set ert_crit >= 1.');
end


%--------------------------------------------------------------------------
% Transfer entropy:

tic

% We calculate the transfer entropy T as a sum of 4 terms, T1, T2, T3, T4:
T  = zeros(ntotal,ntotal);
T1 = zeros(ntotal,1);
T2 = zeros(ntotal,1);
T3 = zeros(ntotal,ntotal);
T4 = zeros(ntotal,ntotal);
for i=1:ntotal
    for j=1:ntotal
        pb = 5;
        [~,fracn,T1(j)] = ...
            estimateH2(x(1+taumin(i,j):end,j),x(1:end-taumin(i,j),j),pb,q); 
        while fracn < fraction; 
            pb = pb+1;
            [~,fracn,T1(j)] = ...
                estimateH2(x(1+taumin(i,j):end,j),x(1:end-taumin(i,j),j),pb,q);  
        end
        pb = 5;
        [~,fracn,T2(j)] = ...
            estimateH2(x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),j),pb,q); 
        while fracn < fraction; 
            pb = pb+1;
            [~,fracn,T2(j)] = ...
                estimateH2(x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),j),pb,q); 
        end
        pb = 5;
        [T3(i,j),fracn] = ...
            estimateH3(x(1+taumin(i,j):end,j),x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),i),pb,q); 
        while fracn < fraction; 
            pb = pb+1;
            [T3(i,j),fracn] = ...
                estimateH3(x(1+taumin(i,j):end,j),x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),i),pb,q);
        end
        pb = 5;
        [~,fracn,T4(i,j)] = ...
            estimateH2(x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),i),pb,q); 
        while fracn < fraction; 
            pb = pb+1;
            [~,fracn,T4(i,j)] = ...
                estimateH2(x(1:end-taumin(i,j),j),x(1:end-taumin(i,j),i),pb,q);   
        end
    end
end
for i=1:ntotal
    for j=1:ntotal
        T(i,j) = T1(j) - T2(j) - T3(j,i) + T4(j,i);
    end
end


%--------------------------------------------------------------------------
% Create 'con_array', the array of connections among entities:

con_array = zeros(ntotal,ntotal);

if exist('ERT_array_2','var') == 1
    for i=1:ntotal
        con_array(i,ERT_array_2(i,2)) = ...
        rel_entr_reduc_2(i,ERT_array_2(i,2));
    end
end

if exist('ERT_array_3','var') == 1
    for i=1:ntotal
        con_array(i,ERT_array_3(i,3)) = ...
        rel_entr_reduc_3(i,ERT_array_2(i,2),ERT_array_3(i,3));
    end
end

if exist('ERT_array_4','var') == 1
    for i=1:ntotal
        con_array(i,ERT_array_4(i,4)) = ...
        rel_entr_reduc_4(i,ERT_array_2(i,2),ERT_array_3(i,3),ERT_array_4(i,4));
    end
end    


%--------------------------------------------------------------------------
% Discard false positives in the connection array:

if threshold == 1
    thresholdOption = 'adapt';
    maxMI = max(max(con_array));
    if maxMI < 0.3, threshold = 0.0; 
    else
        if maxMI > 0.7, threshold = 0.2;
        else
            threshold = 0.5*(maxMI-0.3);
        end
    end
else
    thresholdOption = 'fixed';
end

for i=1:ntotal
    for j=1:ntotal
        if abs(con_array(i,j)) < threshold
            con_array(i,j) = 0;
        end
    end
end

fprintf(1,'\n Transfer entropies and connection array calculated in %d seconds.\n',toc); 


%--------------------------------------------------------------------------
% Return output structure:

% Connection array, needed for visualization:
Output.con_array = con_array;

% Additional results, may be used for further analysis:
Output.MI  = MI;
Output.MIl = MIl;
Output.MIm = MIm;
Output.MIs = MIs;
Output.H1  = H1;
Output.H2  = H2;
Output.T   = T;
if ert_crit >= 2
    Output.H3  = H3;
    Output.MI3 = MI3;
    if ert_crit >= 3
        Output.H4  = H4;
    end
end
Output.taumin  = taumin -1;
Output.dist = dist;
Output.cond_entr2 = cond_entr2;
if ert_crit >= 2, Output.cond_entr3 = cond_entr3; end
if ert_crit >= 3, Output.cond_entr4 = cond_entr4; end
if strcmpi(thresholdOption,'adapt')
    Output.adaptThres = threshold;
end
