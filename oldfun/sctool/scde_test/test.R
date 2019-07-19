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

