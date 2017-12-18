function [x,y]=quadraticcurveto(currentp,curvetop,controllp)

%The quadraticCurveTo method creates a line from the path's current point to the specified point, via a controlpoint.
%The quadraticCurveTo method takes four parameters, the x and y coordinates for the controlpoint, and the x and y coordinates for the line's destination.

%Current point moveTo(20,20)
%CurveTo point quadraticCurveTo(20,100,200,20)
%Controllpoint quadraticCurveTo(20,100,200,20)

o=curvetop(1);
i=curvetop(2);
d=currentp(1);
f=currentp(2);


x2=curvetop(1);
y2=curvetop(2);

x1=currentp(1);
y1=currentp(2);


if nargin<3
    
    x3=min([x1,x2])+abs(x1-x2)/2;
    y3=y1;
    controllp=[x3 y3];
    controllp=[(d+o)/2+(i-f)/4,(f+i)/2+(d-o)/4];
end
P=[[currentp';0],[controllp';0],[curvetop';0]];

n=3;
count = 1;

div = 50; %number of segments of the curve (Increase this value to obtain a
          %smoother curve

for u = 0:(1/div):1
    sum = [0 0 0]';
    for i = 1:n
        B = nchoosek(n,i-1)*(u^(i-1))*((1-u)^(n-i+1)); %B is the Bernstein polynomial value
        sum = sum + B*P(:,i);
    end
    B = nchoosek(n,n)*(u^(n));
    sum = sum + B*P(:,n);
    A(:,count) = sum; %the matrix containing the points of curve as column vectors. 
    count = count+1;  % count is the index of the points on the curve.
end

%plotting the curve
x = A(1,:);
y = A(2,:);

%plot(currentp,curvetop,'og-');

%{
for j = 1:n %plots the points
    plot(P(1,j),P(2,j),'*');
    hold on;
end
hold on
plot(x,y)
%}
