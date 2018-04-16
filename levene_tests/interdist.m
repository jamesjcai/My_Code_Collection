function D = interdist(A)
  n=size(A,2);
  for k=1:n
    D(:,:,k)=A(:,k)-A(:,k)';
  end
  D=sqrt(sum(D.^2,3));
end
