
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
      
      

