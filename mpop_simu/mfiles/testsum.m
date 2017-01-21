
t1=zeros(1,1000);
for k=1:1000
    k
t1(k)=tajima89d_test(ms_mex(400,4,0)+1);
end

t2=zeros(1,1000);
for j=1:1000
    j
t2(j)=tajima89d_test(simucode+1);
end

