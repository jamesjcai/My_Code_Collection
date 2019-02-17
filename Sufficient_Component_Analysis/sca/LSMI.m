function [score_cv,alphah,Phiy,Phix,u,sigmax_chosen]=LSMI(x,y,y_type,sigmax_list,lambda_list,b,fold)
%
% Least-Squares Mutual Information (with likelihood cross validation)
%
% Estimating a squared-loss variant of mutual information
%    \frac{1}{2}\int\int (\frac{p_{xy}(x,y)}{p_x(x)p_y(y)}-1)^2 p_x(x)p_y(y) dx dy
% from input-output samples
%    { (x_i,y_i) | x_i\in R^{dx}, y_i\in R^{dy} }_{i=1}^n
% drawn independently from a joint density p_{xy}(x,y).
% p_x(x) and p_y(y) are marginal densities of x and y, respectively.
%
% Usage:
%       [MIh,score_cv]=LSMI(x,y,y_type,sigmax_list,lambda_list,b)
%
% Input:
%    x          : dx by n input sample matrix
%    y          : dy by n output sample matrix
%    y_type     : if y_type=1, delta kernel is used for y;
%                 otherwise (or empty) Gaussian kernel is used.
%    sigmax_list : (candidates of) Gaussian width
%                 If sigmax_list is a vector, one of them is selected by cross validation.
%                 If sigmax_list is a scalar, this value is used without cross validation
%                 If sigmax_list is empty/undefined, Gaussian width is chosen from
%                 some default canditate list by cross validation
%    lambda_list: (OPTIONAL) regularization parameter
%                 If lambda_list is a vector, one of them is selected by cross validation.
%                 If lambda_list is a scalar, this value is used without cross validation
%                 If lambda_list is empty, Gaussian width is chosen from
%                 some default canditate list by cross validation
%    b          : number of Gaussian centers (if empty, b=200 is used);
%
% Output:
%    MIh        : estimated mutual information between x and y
%    score_cv   : cross validation score
%
% (c) Taiji Suzuki, Department of Mathematical Informatics, The University of Tokyo, Japan.
%     Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     s-taiji@stat.t.u-tokyo.ac.jp
%     sugi@cs.titech.ac.jp,

if nargin<2
    error('number of input arguments is not enough!!!')
end

n =size(x,2);
ny=size(y,2);
if n ~= ny
    error('x and y must have the same number of samples!!!')
end

%y_type
if nargin<4
    y_type=0;
end

if nargin < 5 || isempty(sigmax_list)
    medx = compmedDist(x');
    medy = compmedDist(y');
    
end

if nargin < 6 || isempty(lambda_list)
    lambda_list = [0.01 0.001];
end

if nargin<7 || isempty(b)
    b = 500;
end

if nargin<8 || isempty(fold)
    fold=5;
end

sigmax_list  = [0.25 0.5 0.75 1.0]*medx;
sigmay_list  = [0.25 0.5 0.75 1.0]*medy;
%sigmay = sigmay_list(sigmay_ind);

b=min(n,b);
%b = n;
rand('state',0);
randn('state',0);

rand_index=randperm(n);
rand_index_base = rand_index(1:b);
u=x(:,rand_index_base);
v=y(:,rand_index_base);

Phix_tmp=GaussBasis_sub(x,u)';
Phiy_tmp=GaussBasis_sub(y,v)';

if length(sigmax_list)==1 && length(lambda_list)==1
    sigma_chosen=sigmax_list;
    lambda_chosen=lambda_list;
    score_cv=-inf;
else
    %%%%%%%%%%%%%%%% Searching Gaussian kernel width `sigma_chosen'
    %%%%%%%%%%%%%%%% and regularization parameter `lambda_chosen'
    rand('state',0);
    randn('state',0);
    fold_index=[1:fold];
    cv_index=randperm(n);
    cv_split=floor([0:n-1]*fold./n)+1;
    
    
    %keyboard
    for sigmay_index = 1:length(sigmay_list)
        scores_cv=zeros(length(sigmax_list),length(lambda_list));
        sigmay = sigmay_list(sigmay_index);
        
        Phiy_sigma=Gaussian_kernel(Phiy_tmp,sigmay);
       
        
        HH_cv_Phiy = zeros(size(Phiy_sigma,1),size(Phiy_sigma,1),fold);
        for i=fold_index
            cv_index_tmp=cv_index(cv_split==i);
            HH_cv_Phiy(:,:,i)=Phiy_sigma(:,cv_index_tmp)*Phiy_sigma(:,cv_index_tmp)';
        end
        
        Tmpy = sum(HH_cv_Phiy,3);
        for sigmax_index=1:length(sigmax_list)
            sigma=sigmax_list(sigmax_index);
            Phix_sigma=Epanechnikov_kernel(Phix_tmp,sigma);
            
            Phi_sigma=Phix_sigma.*Phiy_sigma;
            
            HH_cv_Phix = zeros(size(Phi_sigma,1),size(Phi_sigma,1),fold);
            
            for i=fold_index
                cv_index_tmp=cv_index(cv_split==i);
                
                HH_cv_Phix(:,:,i)=Phix_sigma(:,cv_index_tmp)*Phix_sigma(:,cv_index_tmp)';
                hh_cv(:,:,i)=mymean(Phi_sigma(:,cv_index_tmp),2);
                n_cv(i)=length(cv_index_tmp);
            end
            
            Tmpx = sum(HH_cv_Phix,3);
            
            for i=fold_index
                cv_index_tr=fold_index(fold_index~=i);
                Hh_cv_tr = 1/(sum(n_cv(cv_index_tr))^2)*(Tmpx - HH_cv_Phix(:,:,i)).*(Tmpy - HH_cv_Phiy(:,:,i));
                Hh_cv_te=1/(n_cv(i)^2)*HH_cv_Phix(:,:,i).*HH_cv_Phiy(:,:,i);
                hh_cv_tr=mymean(hh_cv(:,:,cv_index_tr),3);
                hh_cv_te=hh_cv(:,:,i);
                for lambda_index=1:length(lambda_list)
                    lambda=lambda_list(lambda_index);
                    Tmp = Hh_cv_tr+lambda*eye(b);
                    alphah_cv=max(0,mylinsolve(Tmp,hh_cv_tr));
                    wh_cv=alphah_cv'*Hh_cv_te*alphah_cv/2-hh_cv_te'*alphah_cv;
                    scores_cv(sigmax_index,lambda_index)=scores_cv(sigmax_index,lambda_index)+wh_cv/fold;
                end % fold
            end % lambda
        end % sigma
        
        [scores_cv_tmp,lambda_chosen_index]=min(scores_cv,[],2);
        [scores_cv_all(sigmay_index),sigma_chosen_index]=min(scores_cv_tmp);
        lambda_chosen_tmp(sigmay_index)=lambda_list(lambda_chosen_index(sigma_chosen_index));
        sigmax_chosen_tmp(sigmay_index)=sigmax_list(sigma_chosen_index);
    end
end %length(sigmax_list)==1 && length(lambda_list)==1
%%%%%%%%%%%%%%%% Computing the final solution `MIh'

[score_cv, sigmay_chosen_index] = min(scores_cv_all);
sigmay_chosen = sigmay_list(sigmay_chosen_index);
lambda_chosen = lambda_chosen_tmp(sigmay_chosen_index);
sigmax_chosen = sigmax_chosen_tmp(sigmay_chosen_index);

%keyboard
Phix=Epanechnikov_kernel(Phix_tmp,sigmax_chosen);
Phiy=Gaussian_kernel(Phiy_tmp,sigmay_chosen);

Phi=Phix.*Phiy;
Hh=1/(n^2)*(Phix*Phix').*(Phiy*Phiy');
hh=mean(Phi,2);
%keyboard;
alphah=max(0,mylinsolve(Hh+lambda_chosen*eye(b),hh));
MIh=hh'*alphah/2-1/2;

%score_cv = MIh;
score_cv = -(score_cv + 0.5);
%keyboard
%sigma_chosen
