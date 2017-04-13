load variancePartition_example_data.mat

% form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
% varPart <- fitExtractVarPartModel( geneExpr, form, info )
% > varPart
% Object of class "varPartResults"
%                Batch Individual      Tissue          Age  Residuals
% gene1   1.579636e-04  0.8903737 0.024688515 4.528933e-05 0.08473457
% gene2   0.000000e+00  0.8060315 0.101019178 3.336681e-04 0.09261569
% gene3   2.422515e-03  0.8901117 0.035618136 1.471720e-03 0.07037597

for k=1:1 % size(geneExpr,1)
    y=geneExpr(k,:)';
    
    tbl = table(y, info.Age, info.Batch, info.Height, info.Individual, info.Tissue,...
          'VariableNames',{'y','Age','Batch','Height','Individual','Tissue'});
    lme = fitlme(tbl,'y ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)');
    
end
