idx=find(sce.g=="RPS28");
x=sce.X(idx,:);

[~,i]=sort(t_mono1);
X=sce.X(:,i);

x2=X(idx,:);
y=fft(x);
y2=fft(x2);
figure;
subplot(311)
plot(x);
xlim([0 length(x)])
subplot(312)
plot(abs(y(2:end)).^2);
xlim([0 length(x)])
subplot(313)
z=abs(y2(2:end)).^2;
z(z>1e5)=0;
plot(z);
xlim([0 length(x)])


