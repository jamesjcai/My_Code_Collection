function myboxplot(X,g)

boxplot(X,g,'colors','k');
hold on;
plot(g+1+0.025*randn(size(g)),X,'ko');
hold off;
end
