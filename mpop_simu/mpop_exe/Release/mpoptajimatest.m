d=[];
for k=1:100
    system('mpop -i ms.example -o new.example -N 400 -m 0.01 -r 0 -g 1600');
    [a]=snp_readmsoutfile('new.example');
    d(k)=tajima89d_test(a);
    k
end

mean(d)
