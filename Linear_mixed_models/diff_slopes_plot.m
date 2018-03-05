     load carsmall
      ds = dataset(MPG,Weight,Model_Year);
      lme = fitlme(ds,'MPG ~ Weight + (1|Model_Year)')

%       Plot predicted values conditional on each year.
      gscatter(ds.Weight,ds.MPG,ds.Model_Year)
      ds2 = dataset;
      ds2.Weight = linspace(1500,5000)';
      ds2.Model_Year = repmat(70,100,1);
      line(ds2.Weight,predict(lme,ds2),'color','r');
      ds2.Model_Year(:) = 76;
      line(ds2.Weight,predict(lme,ds2),'color','g');
      ds2.Model_Year(:) = 82;
      line(ds2.Weight,predict(lme,ds2),'color','b');