function [A] = fpkm2tpm(A)

% B = A;
% for k=1:size(A,2)
%     c = sum(A);
%     B(:,k)= A(:,k)./c(k);
% end
% B = B*1e6;

% A=A./sum(A)*1e6;

[A]=normalize_libsize(A,1e6);