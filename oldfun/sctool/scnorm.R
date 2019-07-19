dtable<-read.csv('tmp.txt',header=FALSE)

Conditions = rep(c(1), each= 2339)
DataNorm <- SCnorm(Data = dtabe, Conditions=Conditions, PrintProgressPlots = TRUE, FilterCellNum = 10, PropToUse = .1, Thresh = .1, ditherCounts = TRUE)


dftable<-as.data.frame(dtable)

rownames(dftable)<-paste('gene',1:13535)

DataNorm <- SCnorm(Data = dftable, Conditions=Conditions, PrintProgressPlots = TRUE, FilterCellNum = 10, PropToUse = .1, Thresh = .1, ditherCounts = TRUE)
