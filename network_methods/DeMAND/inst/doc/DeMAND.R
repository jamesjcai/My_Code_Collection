### R code from vignette source 'DeMAND.Rnw'

###################################################
### code chunk number 1: DeMAND.Rnw:24-26 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("DeMAND")


###################################################
### code chunk number 2: DeMAND.Rnw:33-34
###################################################
library(DeMAND)


###################################################
### code chunk number 3: DeMAND.Rnw:49-51
###################################################
data(inputExample)
#ls()


###################################################
### code chunk number 4: DeMAND.Rnw:56-58
###################################################
dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
printDeMAND(dobj)


###################################################
### code chunk number 5: DeMAND.Rnw:65-67
###################################################
dobj <- runDeMAND(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
printDeMAND(dobj)


###################################################
### code chunk number 6: DeMAND.Rnw:99-100
###################################################
sessionInfo()


