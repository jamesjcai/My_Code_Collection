


for k=1:1000
    load(sprintf('res/simures%d',k));
    seq1=segseq+1;
    % t=tajima89d_test(segseq+1);
    % t=faywu00h_test(seq,zeros(1,size(seq,2)));
    seq2=ms_mex(400,4,0,4)+1;
    
    t11(k)=tajima89d_test(seq1);
    t12(k)=tajima89d_test(seq2);
    
    t21(k)=faywu00h_test(seq1,zeros(1,size(seq1,2)));
    t22(k)=faywu00h_test(seq2,zeros(1,size(seq2,2)));
    
    t31(k)=thetapi(seq1)/size(seq1,2);
    t32(k)=thetapi(seq2)/size(seq2,2);    
    
    t41(k)=hapdiv_test(seq1);
    t42(k)=hapdiv_test(seq2);    
    k
end

figure;
subplot(2,2,1)
histfit(t11); hold on; histfit(t12);
h = get(gca,'Children');
set(h(2),'FaceColor',[.8 .8 1])
title('tajima d')

subplot(2,2,2)
histfit(t21); hold on; histfit(t22);
h = get(gca,'Children');
set(h(2),'FaceColor',[.8 .8 1])
title('fay&wu h')

subplot(2,2,3)
histfit(t31); hold on; histfit(t32);
h = get(gca,'Children');
set(h(2),'FaceColor',[.8 .8 1])
title('nucdiv')

subplot(2,2,4)
histfit(t41); hold on; histfit(t42);
h = get(gca,'Children');
set(h(2),'FaceColor',[.8 .8 1])
title('hapdiv')

% thetapi(seq3)/size(seq3,2)
% snp_hapfreqxdist(seq3)
% nucdiv(seq3)
