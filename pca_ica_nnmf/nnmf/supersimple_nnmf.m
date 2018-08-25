function [W,H]=supersimple_nnmf(X,k)
[r,c]=size(X); % c is # of samples, r is # of features
H=rand(k,c);
% H(H<eps)=0;
H=max(H,eps);
W=X/H;
% W(W<eps)=0;

W=max(W,eps);


            W=W.*((X*H')./(W*(H*H')));
%             W(W<eps)=0;
                W=max(W,eps);
            H=H.*((W'*X)./(W'*W*H));
%             H(H<eps)=0;
                H=max(H,eps)