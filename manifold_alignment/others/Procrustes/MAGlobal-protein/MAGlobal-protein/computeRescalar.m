%Compute the smallest distance between two distance matrices.
function scal= computeRescalar(X, Y)
  scal=trace(Y'*X)/trace(Y'*Y);
end