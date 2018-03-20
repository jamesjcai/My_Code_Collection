function newpdf = lognormpdf(data,mu,sigma)
newpdf = exp((-(log(data)-mu).^2)/(2*sigma^2))./(data*sigma*sqrt(2*pi));