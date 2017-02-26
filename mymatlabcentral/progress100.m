function progress100(k,N,x)
if nargin<3
    x=100;
end
   if rem(k,x) == 0
        fprintf('%5d ...... %d\n', k, N);
   end