function plotnetcurve(A,xy)

[i,j] = find(A);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';

if isfloat(xy)
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
end
X=X(:); Y=Y(:);

for k=1:length(X)-1
    if ~(isnan(X(k))||isnan(X(k+1)))
    A(1)=X(k); B(1)=Y(k);
    A(2)=X(k+1); B(2)=Y(k+1);
    [a,b]=quadraticcurveto([A(1) B(1)],[A(2) B(2)]);
    plot(a,b);
    hold on
    end
    
end
plot(X,Y);

%{
n=size(A,1);
c=1;
for i=1:n-1
    for j=i+1:n
        if A(i,j)
            X{c}=[(xy(c,1),xy(c,2))];            
            
        c=c+1;
        end
    end
end
 
%}


%function definitions
function [] = lineplot(A,B)

x = [A(1) B(1)]; 
y = [A(2) B(2)]; 
plot(x,y,'-r'); %a dashed red  line 

 
