x=linspace(-30,30,2001);
y1=[];
for k=-12:2:21
    y1=[y1; 1/sqrt(2*pi)*exp(-(x-0.95*k).^2./16.0)];
end

y2=[];
for k=-12:2:21
    y2=[y2; 1/sqrt(2*pi)*exp(-(x-0.001*k).^2./160.0)];
end

y=[0.45*y2;y1];
y=y+0.009*randn(size(y));
%y(y<0.015)=0;

g=[zeros(1,size(y2,1)) ones(1,size(y1,1))];
i=randperm(size(y,1));
figure; i_joyplot3col(y(i,:),g(i),0.45,2);

%%

x=linspace(-30,30,2001);
y1=[];
for k=-12:2:21
    y1=[y1; 1/sqrt(2*pi)*exp(-(x-0.005*k).^2./16.0)];
end

y2=[];
for k=-12:2:21
    y2=[y2; 1/sqrt(2*pi)*exp(-(x-0.91*k).^2./160.0)];
end

y=[0.45*y2;y1];
y=y+0.009*randn(size(y));
%y(y<0.015)=0;

g=[zeros(1,size(y2,1)) ones(1,size(y1,1))];
i=randperm(size(y,1));
figure; i_joyplot3col(y(i,:),g(i),0.45,2);
