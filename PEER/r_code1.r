library(peer)
library("R.matlab")

# expr_Brain_Amygdala_72_67.mat                         
# expr_Brain_AnteriorCingulateCortex_BA24__84_77.mat    
# expr_Brain_Caudate_basalGanglia__117_101.mat          
# expr_Brain_CerebellarHemisphere_105_93.mat            
# expr_Brain_Cerebellum_125_110.mat                     
# expr_Brain_Cortex_114_97.mat                          
# expr_Brain_FrontalCortex_BA9__108_95.mat              
#expr_Brain_Hippocampus_94_82.mat                      
# expr_Brain_Hypothalamus_96_86.mat                     
# expr_Brain_NucleusAccumbens_basalGanglia__113_97.mat  
# expr_Brain_Putamen_basalGanglia__97_85.mat            
# expr_Brain_SpinalCord_cervicalC_1__71_64.mat          
#expr_Brain_SubstantiaNigra_63_54.mat 

data <- readMat("expr_Brain_SubstantiaNigra_63_54.mat")
expr=matrix(unlist(data), ncol = 54, byrow = FALSE)
covs = read.table('cov_Brain_SubstantiaNigra_63_54.txt',header=FALSE)
#covs2 = read.table('../peer_qn/cov_Brain_SubstantiaNigra_63_54.txt',header=FALSE)

model = PEER()
PEER_setPhenoMean(model,as.matrix(t(expr)))
PEER_setAdd_mean(model, TRUE)
PEER_setNk(model,15)
#PEER_setCovariates(model, as.matrix(covs2[,c(1,3,2)]))
PEER_setCovariates(model, as.matrix(covs))

PEER_update(model)
residuals = PEER_getResiduals(model)
factors = PEER_getX(model)
covsx = PEER_getCovariates(model)
weigthx = PEER_getW(model)
#PEER_plotModel(model)

writeMat("peer_Brain_SubstantiaNigra_63_54.mat",residuals=residuals,factors=factors,covsx=covsx,weigthx=weigthx)


