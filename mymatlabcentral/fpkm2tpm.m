function [B] = fpkm2tpm(A)

B = A;
for k=1:size(A,2)
    c = sum(A);
    B(:,k)= A(:,k)./c(k);
end
B = B*1e6;