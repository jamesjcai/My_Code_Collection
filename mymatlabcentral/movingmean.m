function [y]=movingmean(data,windowSize)

methodid=3;

switch methodid
    case 1
y=zeros(1,length(data)-windowSize);
for n= windowSize:length(data),
       y(n)= mean(data(n-windowSize+1:n));
      %y[n] is the filtered signal
end
                                          
    case 2

y=filter(ones(1,windowSize)/windowSize,1,data);

    case 3     % moving var

        m=windowSize;
        x=data;
f=zeros(m,1)+1/m;
y=filter(f,1,x.^2)-filter(f,1,x).^2;
end