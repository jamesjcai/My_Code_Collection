xtilde = zeros(n,n);
% generate sparse FFT data
for i=1:K;
    loopbreak = 0;
    while(~loopbreak)
        I(i) = 0+ceil(rand(1)*n/15);
        J(i) = 0+ceil(rand(1)*n/15);
        if(xtilde(I(i),J(i)) == 0)
            loopbreak = 1;
        end
    end
    IC(i) = randn();
    xtilde(I(i),J(i)) = IC(i);
    F(i) = sqrt(4*rand());
    damping(i) = -rand()*.1;
end
if(saveFLAG)
    save([filename,'_PARMS.mat']);
end