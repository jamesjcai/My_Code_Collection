X = linspace(0, 2*pi, 100);
X = X';
Y = sin(X) + randn(100,1);
foo = fit(X,Y, 'sin1')
scatter(X,Y);
hold on
plot(foo)

% https://www.mathworks.com/matlabcentral/answers/28237-curve-fitting-periodic-function
%%

figure;

 t = (1:50)';
 X = ones(50,3);
 X(:,2) = cos((2*pi)/50*t);
 X(:,3) = sin((2*pi)/50*t);
 y = 2*cos((2*pi)/50*t-pi/4)+randn(size(t));
 y = y(:);
 beta = X\y;
 yhat = beta(1)+beta(2)*cos((2*pi)/50*t)+beta(3)*sin((2*pi)/50*t);
 plot(t,y,'bo');
 hold on
 plot(t,yhat,'r','linewidth',2);

foo = fit(t,y, 'sin1')
plot(foo)