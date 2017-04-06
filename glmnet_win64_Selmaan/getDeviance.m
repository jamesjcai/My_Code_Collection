function [frac, D_model, D_null] = getDeviance(y, yHat, mean_y_train, family)

if nargin < 4
    family = 'Poisson';
end

y = y(:);
mu = yHat(:);

switch family
    case 'Poisson'
        % The equations for this function were derived from first
        % principles (OneNote tag tataka), match online sources and the
        % glmnet output. The mu-y term is usually zero but might not be
        % if different data were used for fitting and testing so it is
        % left in. Sources:
        % http://thestatsgeek.com/2014/04/26/deviance-goodness-of-fit-test-for-poisson-regression/
        % http://stats.stackexchange.com/questions/15730/poisson-deviance-and-what-about-zero-observed-values
            
        D_model = 2 * sum(nanRep(y .* log(y ./ mu), 0) + mu - y);
        D_null = 2 * sum(nanRep(y .* log(y ./ mean_y_train), 0) + mean_y_train - y);
    case 'Gaussian'
        % This is simply the R-squared, https://en.wikipedia.org/wiki/Coefficient_of_determination
        D_model = sum((y-mu).^2);
        D_null = sum((y-mean_y_train).^2);
end

frac = 1 - D_model/D_null;