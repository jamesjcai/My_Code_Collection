library(scde)

X<-read.csv('../X.txt',header = FALSE)
Y<-read.csv('../Y.txt',header = FALSE)
#Z1<-read.csv('Z1.txt',header = FALSE)
#Z2<-read.csv('Z2.txt',header = FALSE)
glist<-read.csv('../glist.txt',header = FALSE)
X<-as.data.frame(X)
rownames(X)<-t(glist)



cd <- clean.counts(X)
knn <- knn.error.models(cd)
varinfo <- pagoda.varnorm(knn, counts = cd)

#


# load example dataset
data(es.mef.small)

# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)  
table(sg)

# clean up the dataset
cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)


# EVALUATION NOT NEEDED
# calculate models
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)




