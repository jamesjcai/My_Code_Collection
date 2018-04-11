load testdata datax xy

figure;
plot(datax(:,1),datax(:,2),'o')
hold on
error_ellipse(cov(datax),mean(datax),'conf',0.95);
plot(xy(1),xy(2),'*r');



plot(0.8182,0.066,'*');

%%

% this should return FALSE
i_95pct_error_ellipse(datax,xy,0.95)

% this should return TRUE
i_95pct_error_ellipse(datax,mean(datax),0.95)

% this should return TRUE
i_95pct_error_ellipse(datax,[0.8182 0.066],0.95)

