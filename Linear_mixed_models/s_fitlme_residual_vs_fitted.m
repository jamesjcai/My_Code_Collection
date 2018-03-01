% see https://www.mathworks.com/help/stats/linear-mixed-effects-models.html
% Fit a model and plot the residuals vs. fitted values, grouped byOrigin.
     load carsmall
     ds = dataset(MPG,Weight,Model_Year);
     lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)');
     r = residuals(lme);
     f = fitted(lme);
     gscatter(f,r,Origin)

     
%%   Example: Fit a model and make a probability plot of the residuals.
%            There are two outliers in the upper right of the plot, which
%            you can examine using the data cursor in the plot figure.
     load carsmall
     obsname = strcat(num2str((1:100)'),{' '},Model);
     ds = dataset(MPG,Weight,Model_Year,'ObsNames',obsname);
     lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
     plotResiduals(lme,'probability')   
     
%%
      load carsmall
      ds = table(MPG,Weight,Model_Year);
      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
      
      lme.randomEffects

      
%%

%   Example: Fit a model. Simulate new random values for the first
%            observation.
     load carsmall
     ds = dataset(MPG,Weight,Model_Year);
     lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)');
     newds = ds(ones(10000,1),:);
     r = random(lme,newds);
     hist(r)
%
%      % These random values share a common new random effect due to
%      % Model_Year, but their variance is comparable to the model value.
     v1 = var(r)
     v2 = lme.MSE

%%
%   Example: Model gas mileage as a function of car weight, with a random
%            effect due to model year. Compare a model having random
%            intercepts with one also having random slopes.
     load carsmall
     ds = dataset(MPG,Weight,Model_Year);
     lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')
     lme2 = fitlme(ds,'MPG ~ Weight + (1 | Model_Year) + (-1 + Weight|Model_Year)')
     compare(lme,lme2)
      
